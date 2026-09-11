// Analisis espejo de process_sum4.C para el scan End envuelto.
//
// Reusa la respuesta SUM4 + leading-edge disponible en
// congruent_sum4_timing.C. Esa cadena incluye el walk propio del leading-edge,
// pero NO implementa una correccion de time-walk (HOOK_WALK pendiente).
//
// Uso:
// root -l -b -q 'analysis/tb_mirror_sigma_vs_x.C("scan_resume_2_audit.csv","results/analysis_sigma_vs_x_2026-06-10")'

#include "congruent_sum4_timing.C"

#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLine.h>

#include <fstream>
#include <map>
#include <numeric>

namespace tbmirror {
constexpr double kSqrt2 = 1.4142135623730951;

struct ScanFile {
    double x = 0.0;
    std::string path;
};

struct MirrorFit {
    double meanNs = 0.0;
    double sigmaPs = 0.0;
    double sigmaErrPs = 0.0;
    double chi2Ndf = 0.0;
    double fwhmPs = 0.0;
    int n = 0;
};

struct LandauFit {
    double mpv = 0.0;
    double sigma = 0.0;
    double chi2Ndf = 0.0;
};

struct Result {
    double x = 0.0;
    MirrorFit raw;
    MirrorFit cut25;
    MirrorFit cut50;
    LandauFit landau;
    double npeLeft = 0.0;
    double npeRight = 0.0;
    double p25 = 0.0;
    double p50 = 0.0;
};

std::vector<std::string> Split(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(field);
    return fields;
}

std::vector<ScanFile> ReadAudit(const char* auditPath) {
    std::ifstream input(auditPath);
    std::string line;
    std::getline(input, line);
    std::map<double, std::string> unique;
    while (std::getline(input, line)) {
        const auto fields = Split(line);
        if (fields.size() < 7 || fields[0].empty() || fields[6] != "OK") continue;
        unique[std::stod(fields[0])] = fields[2];
    }
    std::vector<ScanFile> files;
    for (const auto& item : unique) files.push_back({item.first, item.second});
    return files;
}

double Quantile(std::vector<double> values, double probability) {
    if (values.empty()) return 0.0;
    const size_t index = static_cast<size_t>(probability * (values.size() - 1));
    std::nth_element(values.begin(), values.begin() + index, values.end());
    return values[index];
}

double Fwhm(const TH1D& hist) {
    const double half = hist.GetMaximum() / 2.0;
    int left = -1;
    int right = -1;
    for (int bin = 1; bin <= hist.GetNbinsX(); ++bin) {
        if (hist.GetBinContent(bin) >= half) {
            left = bin;
            break;
        }
    }
    for (int bin = hist.GetNbinsX(); bin >= 1; --bin) {
        if (hist.GetBinContent(bin) >= half) {
            right = bin;
            break;
        }
    }
    if (left < 0 || right <= left) return 0.0;
    return 1000.0 * (hist.GetBinCenter(right) - hist.GetBinCenter(left));
}

MirrorFit FitPeakSeeded(const std::vector<double>& values, const std::string& name) {
    MirrorFit result;
    result.n = static_cast<int>(values.size());
    if (values.size() < 20) return result;

    const double median = Median(values);
    const double robust = std::max(RobustSigma(values), 0.02);
    const double histLo = median - 8.0 * robust;
    const double histHi = median + 8.0 * robust;
    TH1D hist(name.c_str(), "", 200, histLo, histHi);
    hist.SetDirectory(nullptr);
    for (double value : values) hist.Fill(value);

    const int peakBin = hist.GetMaximumBin();
    const double peakX = hist.GetBinCenter(peakBin);
    const double peakY = hist.GetBinContent(peakBin);

    // Adaptacion de process_sum4.C: pico como semilla y ventana local fija.
    // El dato usa 2--6 ns porque su pico vive ahi; DeltaT_LR cambia con x,
    // por lo que trasladamos esa ventana de 4 ns alrededor del pico observado.
    const double fitLo = std::max(histLo, peakX - 2.0);
    const double fitHi = std::min(histHi, peakX + 2.0);
    TF1 gaus((name + "_gaus").c_str(), "gaus", fitLo, fitHi);
    gaus.SetParameters(peakY, peakX, 0.5);
    gaus.SetParLimits(0, 0.1, 10.0 * peakY);
    gaus.SetParLimits(1, fitLo, fitHi);
    gaus.SetParLimits(2, 0.02, 5.0);
    const int status = hist.Fit(&gaus, "QNR");
    if (status == 0) {
        result.meanNs = gaus.GetParameter(1);
        result.sigmaPs = 1000.0 * std::abs(gaus.GetParameter(2));
        result.sigmaErrPs = 1000.0 * gaus.GetParError(2);
        result.chi2Ndf = gaus.GetNDF() > 0 ? gaus.GetChisquare() / gaus.GetNDF() : 0.0;
    }
    result.fwhmPs = Fwhm(hist);
    return result;
}

LandauFit FitNpeLandau(const std::vector<double>& npe, const std::string& name) {
    LandauFit result;
    if (npe.size() < 20) return result;
    const double maxNpe = *std::max_element(npe.begin(), npe.end());
    TH1D hist(name.c_str(), "", 100, 0.0, maxNpe * 1.05);
    hist.SetDirectory(nullptr);
    for (double value : npe) hist.Fill(value);
    const int peakBin = hist.GetMaximumBin();
    const double peakX = hist.GetBinCenter(peakBin);
    const double peakY = hist.GetBinContent(peakBin);
    const double binWidth = hist.GetBinWidth(peakBin);
    const double fitLo = std::max(0.0, peakX - 10.0 * binWidth);
    const double fitHi = std::min(maxNpe * 1.05, peakX + 50.0 * binWidth);
    TF1 landau((name + "_landau").c_str(), "landau", fitLo, fitHi);
    landau.SetParameters(peakY, peakX, 7.0 * binWidth);
    landau.SetParLimits(0, 0.01 * peakY, 100.0 * peakY);
    landau.SetParLimits(1, std::max(fitLo, peakX - 16.0 * binWidth),
                       std::min(fitHi, peakX + 16.0 * binWidth));
    landau.SetParLimits(2, binWidth, 36.0 * binWidth);
    const int status = hist.Fit(&landau, "QNR");
    if (status == 0) {
        result.mpv = landau.GetParameter(1);
        result.sigma = std::abs(landau.GetParameter(2));
        result.chi2Ndf = landau.GetNDF() > 0 ? landau.GetChisquare() / landau.GetNDF() : 0.0;
    }
    return result;
}

void Style(TGraphErrors& graph, int color, int marker) {
    graph.SetLineColor(color);
    graph.SetMarkerColor(color);
    graph.SetLineWidth(2);
    graph.SetMarkerStyle(marker);
}

void SaveMaps(double x, const std::vector<double>& delta, const std::vector<double>& npe,
              const std::string& outputDir) {
    if (x != 0.0 && x != 350.0 && x != 690.0) return;
    const auto mm = static_cast<int>(x);
    const double loDt = Quantile(delta, 0.005);
    const double hiDt = Quantile(delta, 0.995);
    const double hiNpe = Quantile(npe, 0.995);
    TH2D map(("map_" + std::to_string(mm)).c_str(), "", 100, 0.0, hiNpe,
             120, loDt, hiDt);
    map.SetDirectory(nullptr);
    for (size_t i = 0; i < delta.size(); ++i) map.Fill(npe[i], delta[i]);
    TCanvas canvas(("map_canvas_" + std::to_string(mm)).c_str(), "", 950, 700);
    canvas.SetRightMargin(0.14);
    map.SetTitle(("Simulacion Geant4 - End readout, x=" + std::to_string(mm)
                  + " mm;NPE_{L}+NPE_{R};#DeltaT_{LR} [ns]").c_str());
    map.Draw("COLZ");
    canvas.SaveAs((outputDir + "/deltaT_vs_NPE_x" + std::to_string(mm) + "_End.png").c_str());
    canvas.SaveAs((outputDir + "/muon_" + std::to_string(mm) + "mm_End.png").c_str());
}

void SaveSummaryPlots(const std::vector<Result>& results, const std::string& outputDir) {
    TGraphErrors raw, cut25, cut50, left, right;
    for (size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        raw.SetPoint(i, r.x, r.raw.sigmaPs / kSqrt2);
        raw.SetPointError(i, 0.0, r.raw.sigmaErrPs / kSqrt2);
        cut25.SetPoint(i, r.x, r.cut25.sigmaPs / kSqrt2);
        cut25.SetPointError(i, 0.0, r.cut25.sigmaErrPs / kSqrt2);
        cut50.SetPoint(i, r.x, r.cut50.sigmaPs / kSqrt2);
        cut50.SetPointError(i, 0.0, r.cut50.sigmaErrPs / kSqrt2);
        left.SetPoint(i, r.x, r.npeLeft);
        right.SetPoint(i, r.x, r.npeRight);
    }
    Style(raw, kBlue + 1, 20);
    Style(cut25, kOrange + 7, 21);
    Style(cut50, kGreen + 2, 22);
    TCanvas sigmaCanvas("sigma_tb_mirror", "", 1100, 700);
    sigmaCanvas.SetGrid();
    raw.SetTitle("Simulacion Geant4 - End readout;Posicion x [mm];#sigma_{single} [ps]");
    raw.Draw("ALP");
    cut25.Draw("LP SAME");
    cut50.Draw("LP SAME");
    TLine reference(-690, 88, 690, 88);
    reference.SetLineColor(kRed + 1);
    reference.SetLineStyle(2);
    reference.SetLineWidth(2);
    reference.Draw();
    TLegend legend(0.15, 0.66, 0.49, 0.88);
    legend.AddEntry(&raw, "Sin corte NPE", "lp");
    legend.AddEntry(&cut25, "NPE #geq p25", "lp");
    legend.AddEntry(&cut50, "NPE #geq p50", "lp");
    legend.AddEntry(&reference, "Test beam (SUM4, reportado)", "l");
    legend.Draw();
    sigmaCanvas.SaveAs((outputDir + "/sigma_vs_x_End.png").c_str());
    sigmaCanvas.SaveAs((outputDir + "/muon_longitudinal_scan_End.png").c_str());

    Style(left, kBlue + 1, 20);
    Style(right, kRed + 1, 21);
    TCanvas npeCanvas("npe_scan", "", 1100, 700);
    npeCanvas.SetGrid();
    left.SetTitle("Simulacion Geant4 - End readout;Posicion x [mm];NPE medio por evento");
    left.Draw("ALP");
    right.Draw("LP SAME");
    TLegend npeLegend(0.40, 0.75, 0.62, 0.88);
    npeLegend.AddEntry(&left, "End-L", "lp");
    npeLegend.AddEntry(&right, "End-R", "lp");
    npeLegend.Draw();
    npeCanvas.SaveAs((outputDir + "/npe_vs_x_End.png").c_str());
}
}  // namespace tbmirror

void tb_mirror_sigma_vs_x(const char* auditPath = "scan_resume_2_audit.csv",
                          const char* outputDir = "results/analysis_sigma_vs_x_2026-06-10") {
    using namespace tbmirror;
    gStyle->SetOptStat(0);
    gSystem->mkdir(outputDir, true);
    const auto files = ReadAudit(auditPath);
    std::vector<Result> results;

    for (size_t fileIndex = 0; fileIndex < files.size(); ++fileIndex) {
        const auto& scan = files[fileIndex];
        TFile file(scan.path.c_str(), "READ");
        auto* tree = dynamic_cast<TTree*>(file.Get("sipm_hits"));
        if (file.IsZombie() || !tree) continue;

        std::vector<EventData> events(10000);
        std::vector<int> npeLeft(10000, 0), npeRight(10000, 0);
        TTreeReader reader(tree);
        TTreeReaderValue<int> eventId(reader, "event_id");
        TTreeReaderValue<int> globalId(reader, "global_id");
        TTreeReaderValue<double> time(reader, "time_ns");
        while (reader.Next()) {
            if (*eventId < 0 || *eventId >= 10000 || *globalId < 0 || *globalId > 15) continue;
            const int group = *globalId / 4;
            events[*eventId].groups[group].push_back(*time);
            if (*globalId < 8) ++npeLeft[*eventId];
            else ++npeRight[*eventId];
        }

        std::vector<double> delta, npe;
        delta.reserve(10000);
        npe.reserve(10000);
        long long totalLeft = 0;
        long long totalRight = 0;
        for (int event = 0; event < 10000; ++event) {
            totalLeft += npeLeft[event];
            totalRight += npeRight[event];
            const double l0 = LeadingEdgeTime(std::move(events[event].groups[0]));
            const double l1 = LeadingEdgeTime(std::move(events[event].groups[1]));
            const double r0 = LeadingEdgeTime(std::move(events[event].groups[2]));
            const double r1 = LeadingEdgeTime(std::move(events[event].groups[3]));
            const double left = Earliest(l0, l1);
            const double right = Earliest(r0, r1);
            if (!std::isfinite(left) || !std::isfinite(right)) continue;
            delta.push_back(left - right);
            npe.push_back(npeLeft[event] + npeRight[event]);
        }

        Result result;
        result.x = scan.x;
        result.npeLeft = static_cast<double>(totalLeft) / 10000.0;
        result.npeRight = static_cast<double>(totalRight) / 10000.0;
        result.p25 = Quantile(npe, 0.25);
        result.p50 = Quantile(npe, 0.50);
        result.raw = FitPeakSeeded(delta, "dt_raw_" + std::to_string(fileIndex));
        std::vector<double> dt25, dt50;
        for (size_t i = 0; i < delta.size(); ++i) {
            if (npe[i] >= result.p25) dt25.push_back(delta[i]);
            if (npe[i] >= result.p50) dt50.push_back(delta[i]);
        }
        result.cut25 = FitPeakSeeded(dt25, "dt_p25_" + std::to_string(fileIndex));
        result.cut50 = FitPeakSeeded(dt50, "dt_p50_" + std::to_string(fileIndex));
        result.landau = FitNpeLandau(npe, "npe_" + std::to_string(fileIndex));
        SaveMaps(scan.x, delta, npe, outputDir);
        results.push_back(result);
        std::cout << "x=" << scan.x << " mm: N_eff=" << result.raw.n
                  << ", sigma_single=" << result.raw.sigmaPs / tbmirror::kSqrt2 << " ps\n";
    }

    std::sort(results.begin(), results.end(), [](const Result& a, const Result& b) {
        return a.x < b.x;
    });
    std::ofstream csv(std::string(outputDir) + "/sigma_vs_x_End.csv");
    csv << "posicion_mm,sigma_LR_ps,err_sigma_LR_ps,sigma_single_ps,err_sigma_single_ps,"
           "fwhm_ps,sigma_single_cut_p25_ps,sigma_single_cut_p50_ps,NPE_L,NPE_R,MPV_NPE,"
           "chi2_ndf_gaus,chi2_ndf_landau\n";
    csv << std::fixed << std::setprecision(6);
    for (const auto& r : results) {
        csv << r.x << "," << r.raw.sigmaPs << "," << r.raw.sigmaErrPs << ","
            << r.raw.sigmaPs / tbmirror::kSqrt2 << ","
            << r.raw.sigmaErrPs / tbmirror::kSqrt2 << ","
            << r.raw.fwhmPs << "," << r.cut25.sigmaPs / tbmirror::kSqrt2 << ","
            << r.cut50.sigmaPs / tbmirror::kSqrt2 << "," << r.npeLeft << "," << r.npeRight << ","
            << r.landau.mpv << "," << r.raw.chi2Ndf << "," << r.landau.chi2Ndf << "\n";
    }
    SaveSummaryPlots(results, outputDir);
    std::cout << "TB-mirror analysis complete: " << results.size() << " positions\n";
}
