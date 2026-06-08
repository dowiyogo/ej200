// Analisis endurecido del scan longitudinal EJ-204 de alta estadistica.
//
// CFD@14% se implementa como proxy sobre arribos discretos: tiempo del
// k-esimo fotoelectron, k = round(0.14*N_pe). No se reconstruye waveform.
// End usa (t_L+t_R)/2 y exige >=10 PE en cada lado por evento. FPT se
// conserva como cross-check con el mismo corte End. El run configura
// /sipm/jitterSigma 0 ns: todas las resoluciones son INTRINSECAS.
//
// Eficiencia = fotones detectados / fotones de centelleo generados, usando
// los contadores master impresos por RunAction al final de cada beamOn.
//
// Uso:
//   root -l -b -q 'analysis/high_stats_position_scan.C("results/scan_hi_2026-06-08",10000)'

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace {
constexpr int kMinEndPe = 10;
constexpr double kCfdFraction = 0.14;

struct EventData {
    std::vector<double> left;
    std::vector<double> right;
    std::vector<double> top;
};

struct RunCounters {
    long long scint = 0;
    long long detected = 0;
};

struct FitResult {
    double sigmaPs = 0.0;
    double sigmaErrPs = 0.0;
    double chi2Ndf = 0.0;
    double rmsPs = 0.0;
    int n = 0;
};

struct Point {
    int run = 0;
    double x = 0.0;
    double yieldLeft = 0.0;
    double yieldRight = 0.0;
    double yieldTop = 0.0;
    double yieldLeftErr = 0.0;
    double yieldRightErr = 0.0;
    double yieldTopErr = 0.0;
    double efficiency = 0.0;
    bool endReliable = false;
    FitResult endCfd;
    FitResult topCfd;
    FitResult endFpt;
    FitResult topFpt;
    std::vector<double> endCfdTimes;
    std::vector<double> topCfdTimes;
    std::vector<double> endFptTimes;
    std::vector<double> topFptTimes;
};

struct FoldedPoint {
    double absX = 0.0;
    bool endReliable = false;
    FitResult endCfd;
    FitResult topCfd;
    FitResult endFpt;
    FitResult topFpt;
};

double Mean(const std::vector<double>& values) {
    if (values.empty()) return 0.0;
    double sum = 0.0;
    for (double value : values) sum += value;
    return sum / values.size();
}

double StdDev(const std::vector<double>& values) {
    if (values.size() < 2) return 0.0;
    const double mean = Mean(values);
    double sum = 0.0;
    for (double value : values) sum += (value - mean) * (value - mean);
    return std::sqrt(sum / (values.size() - 1));
}

double Median(std::vector<double> values) {
    if (values.empty()) return 0.0;
    const size_t middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    double median = values[middle];
    if (values.size() % 2 == 0) {
        const auto lower = std::max_element(values.begin(), values.begin() + middle);
        median = 0.5 * (median + *lower);
    }
    return median;
}

double RobustSigma(const std::vector<double>& values) {
    const double median = Median(values);
    std::vector<double> deviations;
    deviations.reserve(values.size());
    for (double value : values) deviations.push_back(std::abs(value - median));
    return 1.4826 * Median(std::move(deviations));
}

double CfdProxy(std::vector<double>& times) {
    const size_t k = std::max<size_t>(
        1, static_cast<size_t>(std::lround(kCfdFraction * times.size())));
    std::nth_element(times.begin(), times.begin() + k - 1, times.end());
    return times[k - 1];
}

FitResult FitCore(const std::vector<double>& times, const std::string& name) {
    FitResult result;
    result.n = static_cast<int>(times.size());
    result.rmsPs = 1000.0 * StdDev(times);
    if (times.size() < 20) return result;

    const double center = Median(times);
    const double robustSigma = std::max(RobustSigma(times), 1e-4);
    const double lo = center - 8.0 * robustSigma;
    const double hi = center + 8.0 * robustSigma;

    TH1D hist(name.c_str(), "", 100, lo, hi);
    hist.SetDirectory(nullptr);
    for (double time : times) hist.Fill(time);

    double mean = center;
    double sigma = robustSigma;
    TF1 gaus((name + "_gaus").c_str(), "gaus", lo, hi);
    int status = 1;
    for (int iteration = 0; iteration < 4; ++iteration) {
        const double fitLo = std::max(lo, mean - 2.0 * sigma);
        const double fitHi = std::min(hi, mean + 2.0 * sigma);
        gaus.SetRange(fitLo, fitHi);
        gaus.SetParameters(hist.GetMaximum(), mean, sigma);
        gaus.SetParLimits(2, 1e-6, hi - lo);
        status = hist.Fit(&gaus, "QNR");
        if (status != 0 || gaus.GetParameter(2) <= 0.0) break;
        mean = gaus.GetParameter(1);
        sigma = std::abs(gaus.GetParameter(2));
    }

    if (status == 0 && sigma > 0.0) {
        result.sigmaPs = 1000.0 * sigma;
        result.sigmaErrPs = 1000.0 * gaus.GetParError(2);
        result.chi2Ndf = gaus.GetNDF() > 0 ? gaus.GetChisquare() / gaus.GetNDF() : 0.0;
    } else {
        result.sigmaPs = result.rmsPs;
        result.sigmaErrPs = result.rmsPs / std::sqrt(2.0 * (times.size() - 1));
        result.chi2Ndf = -1.0;
    }
    return result;
}

std::map<int, RunCounters> ParseMasterCounters(const std::string& logPath, int nEvents) {
    std::ifstream input(logPath);
    std::map<int, RunCounters> counters;
    std::string line;
    int currentRun = -1;
    int events = -1;
    long long scint = 0;
    long long detected = 0;
    const std::regex integerPattern(R"(([-+]?[0-9]+))");
    auto lastInteger = [&](const std::string& text) {
        long long value = 0;
        for (std::sregex_iterator it(text.begin(), text.end(), integerPattern), end;
             it != end; ++it) {
            value = std::stoll((*it)[1]);
        }
        return value;
    };

    while (std::getline(input, line)) {
        if (line.find("Run ID") != std::string::npos) {
            currentRun = static_cast<int>(lastInteger(line));
        } else if (line.find("Events run") != std::string::npos) {
            events = static_cast<int>(lastInteger(line));
        } else if (line.find("Scint photons generated") != std::string::npos) {
            scint = lastInteger(line);
        } else if (line.find("Total photons detected") != std::string::npos) {
            detected = lastInteger(line);
            if (currentRun >= 0 && events == nEvents) counters[currentRun] = {scint, detected};
        }
    }
    return counters;
}

void Style(TGraphErrors& graph, int color, int marker) {
    graph.SetLineColor(color);
    graph.SetMarkerColor(color);
    graph.SetLineWidth(2);
    graph.SetMarkerStyle(marker);
    graph.SetMarkerSize(0.9);
}

void SaveYieldPlot(const std::vector<Point>& points, const std::string& base) {
    TGraphErrors left(points.size()), right(points.size()), top(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        left.SetPoint(i, points[i].x, points[i].yieldLeft);
        right.SetPoint(i, points[i].x, points[i].yieldRight);
        top.SetPoint(i, points[i].x, points[i].yieldTop);
        left.SetPointError(i, 0.0, points[i].yieldLeftErr);
        right.SetPointError(i, 0.0, points[i].yieldRightErr);
        top.SetPointError(i, 0.0, points[i].yieldTopErr);
    }
    Style(left, kBlue + 1, 20);
    Style(right, kRed + 1, 21);
    Style(top, kGreen + 2, 22);
    TCanvas canvas("yield_hi_canvas", "yield", 1100, 650);
    canvas.SetGrid();
    TMultiGraph multi;
    multi.SetTitle("EJ-204: yield detectado sin plegar;Posicion x [mm];Fotoelectrones detectados / evento");
    multi.Add(&left, "LP");
    multi.Add(&right, "LP");
    multi.Add(&top, "LP");
    multi.Draw("A");
    TLegend legend(0.68, 0.72, 0.90, 0.89);
    legend.AddEntry(&left, "End-left (IDs 0-7)", "lp");
    legend.AddEntry(&right, "End-right (IDs 8-15)", "lp");
    legend.AddEntry(&top, "Top (IDs 16-35)", "lp");
    legend.Draw();
    canvas.SaveAs((base + "/yield_vs_position.pdf").c_str());
    canvas.SaveAs((base + "/yield_vs_position.png").c_str());
}

void SaveEfficiencyPlot(const std::vector<Point>& points, const std::string& base) {
    TGraphErrors graph(points.size());
    for (size_t i = 0; i < points.size(); ++i) graph.SetPoint(i, points[i].x, points[i].efficiency);
    Style(graph, kOrange + 7, 20);
    TCanvas canvas("eff_hi_canvas", "efficiency", 1100, 650);
    canvas.SetGrid();
    graph.SetTitle("EJ-204: eficiencia de deteccion;Posicion x [mm];Fotones detectados / fotones de centelleo [%]");
    graph.Draw("ALP");
    graph.SetMinimum(0.0);
    canvas.SaveAs((base + "/efficiency_vs_position.pdf").c_str());
    canvas.SaveAs((base + "/efficiency_vs_position.png").c_str());
}

void SaveFoldedTimingPlot(const std::vector<FoldedPoint>& points, const std::string& base) {
    TGraphErrors end, top, unreliable;
    int iEnd = 0;
    int iTop = 0;
    int iBad = 0;
    for (const auto& point : points) {
        top.SetPoint(iTop, point.absX, point.topCfd.sigmaPs);
        top.SetPointError(iTop++, 0.0, point.topCfd.sigmaErrPs);
        if (point.endReliable) {
            end.SetPoint(iEnd, point.absX, point.endCfd.sigmaPs);
            end.SetPointError(iEnd++, 0.0, point.endCfd.sigmaErrPs);
        } else {
            unreliable.SetPoint(iBad, point.absX, point.endCfd.sigmaPs);
            unreliable.SetPointError(iBad++, 0.0, point.endCfd.sigmaErrPs);
        }
    }
    Style(end, kViolet + 1, 20);
    Style(top, kGreen + 2, 22);
    Style(unreliable, kRed + 1, 24);
    TCanvas canvas("timing_folded_canvas", "timing folded", 1100, 650);
    canvas.SetGrid();
    TMultiGraph multi;
    multi.SetTitle("EJ-204: resolucion temporal intrinseca CFD@14% plegada;|x| [mm];#sigma_{t} intrinseco [ps]");
    multi.Add(&end, "LP");
    multi.Add(&top, "LP");
    multi.Add(&unreliable, "P");
    multi.Draw("A");
    multi.SetMinimum(0.0);
    TLegend legend(0.56, 0.70, 0.90, 0.89);
    legend.AddEntry(&end, "End CFD@14% (fiable)", "lp");
    legend.AddEntry(&top, "Top CFD@14%", "lp");
    legend.AddEntry(&unreliable, "End: min yield < 10 PE", "p");
    legend.Draw();
    canvas.SaveAs((base + "/sigma_cfd_folded_vs_abs_position.pdf").c_str());
    canvas.SaveAs((base + "/sigma_cfd_folded_vs_abs_position.png").c_str());
}

void SaveUnfoldedTimingPlot(const std::vector<Point>& points, const std::string& base) {
    TGraphErrors end, top, unreliable;
    int iEnd = 0;
    int iTop = 0;
    int iBad = 0;
    for (const auto& point : points) {
        top.SetPoint(iTop, point.x, point.topCfd.sigmaPs);
        top.SetPointError(iTop++, 0.0, point.topCfd.sigmaErrPs);
        if (point.endReliable) {
            end.SetPoint(iEnd, point.x, point.endCfd.sigmaPs);
            end.SetPointError(iEnd++, 0.0, point.endCfd.sigmaErrPs);
        } else {
            unreliable.SetPoint(iBad, point.x, point.endCfd.sigmaPs);
            unreliable.SetPointError(iBad++, 0.0, point.endCfd.sigmaErrPs);
        }
    }
    Style(end, kViolet + 1, 20);
    Style(top, kGreen + 2, 22);
    Style(unreliable, kRed + 1, 24);
    TCanvas canvas("timing_unfolded_canvas", "timing unfolded", 1100, 650);
    canvas.SetGrid();
    TMultiGraph multi;
    multi.SetTitle("EJ-204: CFD@14% sin plegar (control de simetria);Posicion x [mm];#sigma_{t} intrinseco [ps]");
    multi.Add(&end, "LP");
    multi.Add(&top, "LP");
    multi.Add(&unreliable, "P");
    multi.Draw("A");
    multi.SetMinimum(0.0);
    TLegend legend(0.56, 0.70, 0.90, 0.89);
    legend.AddEntry(&end, "End CFD@14% (fiable)", "lp");
    legend.AddEntry(&top, "Top CFD@14%", "lp");
    legend.AddEntry(&unreliable, "End: min yield < 10 PE", "p");
    legend.Draw();
    canvas.SaveAs((base + "/sigma_cfd_unfolded_vs_position.pdf").c_str());
    canvas.SaveAs((base + "/sigma_cfd_unfolded_vs_position.png").c_str());
}

void SaveFptCfdPlot(const std::vector<FoldedPoint>& points, const std::string& base) {
    TGraphErrors endCfd, endFpt, topCfd, topFpt;
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& p = points[i];
        endCfd.SetPoint(i, p.absX, p.endCfd.sigmaPs);
        endCfd.SetPointError(i, 0.0, p.endCfd.sigmaErrPs);
        endFpt.SetPoint(i, p.absX, p.endFpt.sigmaPs);
        endFpt.SetPointError(i, 0.0, p.endFpt.sigmaErrPs);
        topCfd.SetPoint(i, p.absX, p.topCfd.sigmaPs);
        topCfd.SetPointError(i, 0.0, p.topCfd.sigmaErrPs);
        topFpt.SetPoint(i, p.absX, p.topFpt.sigmaPs);
        topFpt.SetPointError(i, 0.0, p.topFpt.sigmaErrPs);
    }
    Style(endCfd, kViolet + 1, 20);
    Style(endFpt, kMagenta + 1, 24);
    Style(topCfd, kGreen + 2, 22);
    Style(topFpt, kTeal + 2, 26);
    TCanvas canvas("fpt_cfd_canvas", "FPT vs CFD", 1100, 650);
    canvas.SetGrid();
    TMultiGraph multi;
    multi.SetTitle("EJ-204: cross-check FPT vs CFD@14% (plegado);|x| [mm];#sigma_{t} intrinseco [ps]");
    multi.Add(&endCfd, "LP");
    multi.Add(&endFpt, "LP");
    multi.Add(&topCfd, "LP");
    multi.Add(&topFpt, "LP");
    multi.Draw("A");
    multi.SetMinimum(0.0);
    TLegend legend(0.63, 0.65, 0.90, 0.89);
    legend.AddEntry(&endCfd, "End CFD@14%", "lp");
    legend.AddEntry(&endFpt, "End FPT", "lp");
    legend.AddEntry(&topCfd, "Top CFD@14%", "lp");
    legend.AddEntry(&topFpt, "Top FPT", "lp");
    legend.Draw();
    canvas.SaveAs((base + "/sigma_fpt_vs_cfd_folded.pdf").c_str());
    canvas.SaveAs((base + "/sigma_fpt_vs_cfd_folded.png").c_str());
}

std::vector<FoldedPoint> Fold(const std::vector<Point>& points) {
    std::map<int, std::vector<const Point*>> groups;
    for (const auto& point : points) groups[static_cast<int>(std::lround(std::abs(point.x)))].push_back(&point);
    std::vector<FoldedPoint> folded;
    for (const auto& [absX, members] : groups) {
        std::vector<double> endCfd, topCfd, endFpt, topFpt;
        bool reliable = true;
        for (const Point* point : members) {
            endCfd.insert(endCfd.end(), point->endCfdTimes.begin(), point->endCfdTimes.end());
            topCfd.insert(topCfd.end(), point->topCfdTimes.begin(), point->topCfdTimes.end());
            endFpt.insert(endFpt.end(), point->endFptTimes.begin(), point->endFptTimes.end());
            topFpt.insert(topFpt.end(), point->topFptTimes.begin(), point->topFptTimes.end());
            reliable = reliable && point->endReliable;
        }
        FoldedPoint point;
        point.absX = absX;
        point.endReliable = reliable;
        point.endCfd = FitCore(endCfd, "fold_end_cfd_" + std::to_string(absX));
        point.topCfd = FitCore(topCfd, "fold_top_cfd_" + std::to_string(absX));
        point.endFpt = FitCore(endFpt, "fold_end_fpt_" + std::to_string(absX));
        point.topFpt = FitCore(topFpt, "fold_top_fpt_" + std::to_string(absX));
        folded.push_back(point);
    }
    return folded;
}

void WriteSummary(const std::vector<Point>& points, const std::string& base, int nEvents) {
    std::ofstream csv(base + "/summary.csv");
    csv << "x_mm,n_events,yield_end_left,yield_end_left_err,yield_end_right,yield_end_right_err,"
           "yield_top,yield_top_err,efficiency_percent,end_reliable,"
           "n_eff_end_cfd,sigma_end_cfd_ps,sigma_end_cfd_err_ps,chi2_ndf_end_cfd,rms_end_cfd_ps,"
           "n_eff_top_cfd,sigma_top_cfd_ps,sigma_top_cfd_err_ps,chi2_ndf_top_cfd,rms_top_cfd_ps,"
           "n_eff_end_fpt,sigma_end_fpt_ps,sigma_end_fpt_err_ps,chi2_ndf_end_fpt,rms_end_fpt_ps,"
           "n_eff_top_fpt,sigma_top_fpt_ps,sigma_top_fpt_err_ps,chi2_ndf_top_fpt,rms_top_fpt_ps\n";
    csv << std::fixed << std::setprecision(5);
    for (const auto& p : points) {
        csv << p.x << "," << nEvents << "," << p.yieldLeft << "," << p.yieldLeftErr << ","
            << p.yieldRight << "," << p.yieldRightErr << "," << p.yieldTop << "," << p.yieldTopErr << ","
            << p.efficiency << "," << (p.endReliable ? "true" : "false") << ","
            << p.endCfd.n << "," << p.endCfd.sigmaPs << "," << p.endCfd.sigmaErrPs << "," << p.endCfd.chi2Ndf << "," << p.endCfd.rmsPs << ","
            << p.topCfd.n << "," << p.topCfd.sigmaPs << "," << p.topCfd.sigmaErrPs << "," << p.topCfd.chi2Ndf << "," << p.topCfd.rmsPs << ","
            << p.endFpt.n << "," << p.endFpt.sigmaPs << "," << p.endFpt.sigmaErrPs << "," << p.endFpt.chi2Ndf << "," << p.endFpt.rmsPs << ","
            << p.topFpt.n << "," << p.topFpt.sigmaPs << "," << p.topFpt.sigmaErrPs << "," << p.topFpt.chi2Ndf << "," << p.topFpt.rmsPs << "\n";
    }
}

void WriteFolded(const std::vector<FoldedPoint>& folded, const std::string& base) {
    std::ofstream csv(base + "/folded_summary.csv");
    csv << "abs_x_mm,end_reliable,n_eff_end_cfd,sigma_end_cfd_ps,sigma_end_cfd_err_ps,chi2_ndf_end_cfd,"
           "n_eff_top_cfd,sigma_top_cfd_ps,sigma_top_cfd_err_ps,chi2_ndf_top_cfd,"
           "sigma_end_fpt_ps,sigma_end_fpt_err_ps,sigma_top_fpt_ps,sigma_top_fpt_err_ps\n";
    csv << std::fixed << std::setprecision(5);
    for (const auto& p : folded) {
        csv << p.absX << "," << (p.endReliable ? "true" : "false") << ","
            << p.endCfd.n << "," << p.endCfd.sigmaPs << "," << p.endCfd.sigmaErrPs << "," << p.endCfd.chi2Ndf << ","
            << p.topCfd.n << "," << p.topCfd.sigmaPs << "," << p.topCfd.sigmaErrPs << "," << p.topCfd.chi2Ndf << ","
            << p.endFpt.sigmaPs << "," << p.endFpt.sigmaErrPs << ","
            << p.topFpt.sigmaPs << "," << p.topFpt.sigmaErrPs << "\n";
    }
}

void WriteAsymmetry(const std::vector<Point>& points, const std::string& base) {
    std::map<int, const Point*> negative;
    std::map<int, const Point*> positive;
    for (const auto& p : points) {
        const int absX = static_cast<int>(std::lround(std::abs(p.x)));
        if (p.x < 0) negative[absX] = &p;
        if (p.x > 0) positive[absX] = &p;
    }
    std::ofstream csv(base + "/asymmetry.csv");
    csv << "abs_x_mm,end_reliable,asym_sigma_end_cfd,asym_sigma_top_cfd\n";
    csv << std::fixed << std::setprecision(6);
    for (const auto& [absX, minus] : negative) {
        const auto found = positive.find(absX);
        if (found == positive.end()) continue;
        const Point* plus = found->second;
        auto asym = [](double positiveValue, double negativeValue) {
            const double average = 0.5 * (positiveValue + negativeValue);
            return average != 0.0 ? (positiveValue - negativeValue) / average : 0.0;
        };
        csv << absX << "," << (minus->endReliable && plus->endReliable ? "true" : "false") << ","
            << asym(plus->endCfd.sigmaPs, minus->endCfd.sigmaPs) << ","
            << asym(plus->topCfd.sigmaPs, minus->topCfd.sigmaPs) << "\n";
    }
}
} // namespace

void high_stats_position_scan(const char* resultDir = "results/scan_hi_2026-06-08",
                              int nEvents = 10000) {
    gStyle->SetOptStat(0);
    const std::string base(resultDir);
    const auto counters = ParseMasterCounters(base + "/run.log", nEvents);
    std::vector<Point> points;

    for (int run = 0; run < 100; ++run) {
        std::ostringstream fileName;
        fileName << base << "/photon_hits_run" << std::setw(3) << std::setfill('0') << run << ".root";
        if (gSystem->AccessPathName(fileName.str().c_str())) break;
        TFile file(fileName.str().c_str(), "READ");
        auto* tree = dynamic_cast<TTree*>(file.Get("sipm_hits"));
        if (file.IsZombie() || tree == nullptr) break;

        TTreeReader reader(tree);
        TTreeReaderValue<int> eventId(reader, "event_id");
        TTreeReaderValue<int> face(reader, "face_type");
        TTreeReaderValue<double> time(reader, "time_ns");
        TTreeReaderValue<double> gunX(reader, "gun_x_mm");
        std::vector<EventData> events(nEvents);
        double position = 0.0;
        while (reader.Next()) {
            if (*eventId < 0 || *eventId >= nEvents) continue;
            position = *gunX;
            auto& event = events[*eventId];
            if (*face == 0) event.left.push_back(*time);
            else if (*face == 1) event.right.push_back(*time);
            else if (*face == 2) event.top.push_back(*time);
        }

        Point point;
        point.run = run;
        point.x = position;
        std::vector<double> leftCounts, rightCounts, topCounts;
        leftCounts.reserve(nEvents);
        rightCounts.reserve(nEvents);
        topCounts.reserve(nEvents);
        for (auto& event : events) {
            leftCounts.push_back(event.left.size());
            rightCounts.push_back(event.right.size());
            topCounts.push_back(event.top.size());
            if (event.left.size() >= kMinEndPe && event.right.size() >= kMinEndPe) {
                const double leftFpt = *std::min_element(event.left.begin(), event.left.end());
                const double rightFpt = *std::min_element(event.right.begin(), event.right.end());
                point.endFptTimes.push_back(0.5 * (leftFpt + rightFpt));
                point.endCfdTimes.push_back(0.5 * (CfdProxy(event.left) + CfdProxy(event.right)));
            }
            if (!event.top.empty()) {
                point.topFptTimes.push_back(*std::min_element(event.top.begin(), event.top.end()));
                point.topCfdTimes.push_back(CfdProxy(event.top));
            }
        }

        point.yieldLeft = Mean(leftCounts);
        point.yieldRight = Mean(rightCounts);
        point.yieldTop = Mean(topCounts);
        point.yieldLeftErr = StdDev(leftCounts) / std::sqrt(nEvents);
        point.yieldRightErr = StdDev(rightCounts) / std::sqrt(nEvents);
        point.yieldTopErr = StdDev(topCounts) / std::sqrt(nEvents);
        point.endReliable = std::min(point.yieldLeft, point.yieldRight) >= kMinEndPe;
        point.endCfd = FitCore(point.endCfdTimes, "end_cfd_" + std::to_string(run));
        point.topCfd = FitCore(point.topCfdTimes, "top_cfd_" + std::to_string(run));
        point.endFpt = FitCore(point.endFptTimes, "end_fpt_" + std::to_string(run));
        point.topFpt = FitCore(point.topFptTimes, "top_fpt_" + std::to_string(run));
        const auto counter = counters.find(run);
        if (counter != counters.end() && counter->second.scint > 0) {
            point.efficiency = 100.0 * counter->second.detected / counter->second.scint;
        }
        std::cout << "Run " << run << ": x=" << point.x
                  << " mm, N_eff End/Top=" << point.endCfd.n << "/" << point.topCfd.n << "\n";
        points.push_back(std::move(point));
    }

    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) { return a.x < b.x; });
    const auto folded = Fold(points);
    WriteSummary(points, base, nEvents);
    WriteFolded(folded, base);
    WriteAsymmetry(points, base);
    SaveYieldPlot(points, base);
    SaveEfficiencyPlot(points, base);
    SaveFoldedTimingPlot(folded, base);
    SaveUnfoldedTimingPlot(points, base);
    SaveFptCfdPlot(folded, base);
    std::cout << "Analisis completado: " << points.size() << " posiciones; salida en " << base << "\n";
}
