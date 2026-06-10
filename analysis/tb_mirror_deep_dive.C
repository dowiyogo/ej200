// Evidencia intermedia y ablacion para CODEX_EXEC_09.
// Script nuevo: no modifica los analisis existentes.

#include "tb_mirror_sigma_vs_x.C"

#include <TPaveText.h>

namespace deepdive {
struct EventSet {
    std::vector<double> delta;
    std::vector<double> npe;
};

EventSet ReadEvents(const std::string& path) {
    TFile file(path.c_str(), "READ");
    auto* tree = dynamic_cast<TTree*>(file.Get("sipm_hits"));
    std::vector<EventData> events(10000);
    std::vector<int> npeLeft(10000, 0), npeRight(10000, 0);
    TTreeReader reader(tree);
    TTreeReaderValue<int> eventId(reader, "event_id");
    TTreeReaderValue<int> globalId(reader, "global_id");
    TTreeReaderValue<double> time(reader, "time_ns");
    while (reader.Next()) {
        if (*eventId < 0 || *eventId >= 10000 || *globalId < 0 || *globalId > 15) continue;
        events[*eventId].groups[*globalId / 4].push_back(*time);
        if (*globalId < 8) ++npeLeft[*eventId];
        else ++npeRight[*eventId];
    }
    EventSet result;
    for (int event = 0; event < 10000; ++event) {
        const double l0 = LeadingEdgeTime(std::move(events[event].groups[0]));
        const double l1 = LeadingEdgeTime(std::move(events[event].groups[1]));
        const double r0 = LeadingEdgeTime(std::move(events[event].groups[2]));
        const double r1 = LeadingEdgeTime(std::move(events[event].groups[3]));
        const double left = Earliest(l0, l1);
        const double right = Earliest(r0, r1);
        if (!std::isfinite(left) || !std::isfinite(right)) continue;
        result.delta.push_back(left - right);
        result.npe.push_back(npeLeft[event] + npeRight[event]);
    }
    return result;
}

tbmirror::MirrorFit DrawDelta(const EventSet& events, int x, const std::string& outputDir) {
    const double median = Median(events.delta);
    const double robust = std::max(RobustSigma(events.delta), 0.02);
    const double histLo = median - 8.0 * robust;
    const double histHi = median + 8.0 * robust;
    TH1D hist(("dt_" + std::to_string(x)).c_str(), "", 200, histLo, histHi);
    hist.SetDirectory(nullptr);
    for (double value : events.delta) hist.Fill(value);
    const double peakX = hist.GetBinCenter(hist.GetMaximumBin());
    const double peakY = hist.GetMaximum();
    const double fitLo = std::max(histLo, peakX - 2.0);
    const double fitHi = std::min(histHi, peakX + 2.0);
    TF1 gaus(("dt_gaus_" + std::to_string(x)).c_str(), "gaus", fitLo, fitHi);
    gaus.SetParameters(peakY, peakX, 0.5);
    gaus.SetParLimits(0, 0.1, 10.0 * peakY);
    gaus.SetParLimits(1, fitLo, fitHi);
    gaus.SetParLimits(2, 0.02, 5.0);
    hist.Fit(&gaus, "QNR");

    TCanvas canvas(("dt_canvas_" + std::to_string(x)).c_str(), "", 800, 650);
    hist.SetTitle(("Simulacion Geant4 - End readout, x=" + std::to_string(x)
                   + " mm;#DeltaT_{LR} [ns];Eventos").c_str());
    hist.SetLineWidth(2);
    hist.Draw("HIST");
    gaus.SetLineColor(kRed + 1);
    gaus.SetLineWidth(3);
    gaus.Draw("SAME");
    TPaveText box(0.56, 0.62, 0.89, 0.88, "NDC");
    box.SetFillColor(kWhite);
    box.SetTextAlign(12);
    box.AddText(Form("#mu = %.4f ns", gaus.GetParameter(1)));
    box.AddText(Form("#sigma = %.2f #pm %.2f ps", 1000.0 * std::abs(gaus.GetParameter(2)),
                     1000.0 * gaus.GetParError(2)));
    box.AddText(Form("#chi^{2}/ndf = %.2f", gaus.GetNDF() > 0
                     ? gaus.GetChisquare() / gaus.GetNDF() : 0.0));
    box.AddText(Form("N_{eff} = %zu", events.delta.size()));
    box.AddText(Form("fit = [%.2f, %.2f] ns", fitLo, fitHi));
    box.Draw();
    canvas.SaveAs((outputDir + "/deltaT_fit_x" + std::to_string(x) + "_End.png").c_str());

    auto fit = tbmirror::FitPeakSeeded(events.delta, "mirror_repeat_" + std::to_string(x));
    return fit;
}

tbmirror::LandauFit DrawNpe(const EventSet& events, int x, const std::string& outputDir) {
    const double maxNpe = *std::max_element(events.npe.begin(), events.npe.end());
    TH1D hist(("npe_" + std::to_string(x)).c_str(), "", 100, 0.0, maxNpe * 1.05);
    hist.SetDirectory(nullptr);
    for (double value : events.npe) hist.Fill(value);
    const int peakBin = hist.GetMaximumBin();
    const double peakX = hist.GetBinCenter(peakBin);
    const double peakY = hist.GetMaximum();
    const double width = hist.GetBinWidth(peakBin);
    const double fitLo = std::max(0.0, peakX - 10.0 * width);
    const double fitHi = std::min(maxNpe * 1.05, peakX + 50.0 * width);
    TF1 landau(("landau_" + std::to_string(x)).c_str(), "landau", fitLo, fitHi);
    landau.SetParameters(peakY, peakX, 7.0 * width);
    landau.SetParLimits(0, 0.01 * peakY, 100.0 * peakY);
    landau.SetParLimits(1, std::max(fitLo, peakX - 16.0 * width),
                       std::min(fitHi, peakX + 16.0 * width));
    landau.SetParLimits(2, width, 36.0 * width);
    hist.Fit(&landau, "QNR");

    TCanvas canvas(("npe_canvas_" + std::to_string(x)).c_str(), "", 800, 650);
    hist.SetTitle(("Simulacion Geant4 - End readout, x=" + std::to_string(x)
                   + " mm;NPE_{L}+NPE_{R};Eventos").c_str());
    hist.SetLineWidth(2);
    hist.Draw("HIST");
    landau.SetLineColor(kRed + 1);
    landau.SetLineWidth(3);
    landau.Draw("SAME");
    TPaveText box(0.56, 0.65, 0.89, 0.88, "NDC");
    box.SetFillColor(kWhite);
    box.SetTextAlign(12);
    box.AddText(Form("MPV = %.2f NPE", landau.GetParameter(1)));
    box.AddText(Form("#sigma_{Landau} = %.2f NPE", std::abs(landau.GetParameter(2))));
    box.AddText(Form("#chi^{2}/ndf = %.2f", landau.GetNDF() > 0
                     ? landau.GetChisquare() / landau.GetNDF() : 0.0));
    box.AddText(Form("fit = [%.1f, %.1f] NPE", fitLo, fitHi));
    box.Draw();
    canvas.SaveAs((outputDir + "/npe_landau_fit_x" + std::to_string(x) + "_End.png").c_str());
    return tbmirror::FitNpeLandau(events.npe, "landau_repeat_" + std::to_string(x));
}
}  // namespace deepdive

void tb_mirror_deep_dive(const char* auditPath = "scan_resume_2_audit.csv",
                         const char* outputDir =
                             "results/analysis_sigma_vs_x_2026-06-10/deep_dive") {
    using namespace deepdive;
    gStyle->SetOptStat(0);
    gSystem->mkdir(outputDir, true);
    const auto files = tbmirror::ReadAudit(auditPath);
    std::ofstream metrics(std::string(outputDir) + "/event_metrics.csv");
    metrics << "position_mm,mean_delta_ns,n_eff,fit_mean_ns,fit_sigma_lr_ps,"
               "fit_sigma_lr_err_ps,chi2_ndf_gaus,fwhm_ps,mpv_npe,sigma_landau,"
               "chi2_ndf_landau\n";
    metrics << std::fixed << std::setprecision(6);

    std::map<int, EventSet> representative;
    for (const auto& scan : files) {
        auto events = ReadEvents(scan.path);
        const auto fit = tbmirror::FitPeakSeeded(events.delta,
            "all_metric_" + std::to_string(static_cast<int>(scan.x)));
        const auto landau = tbmirror::FitNpeLandau(events.npe,
            "all_landau_" + std::to_string(static_cast<int>(scan.x)));
        metrics << scan.x << "," << Mean(events.delta) << "," << fit.n << ","
                << fit.meanNs << "," << fit.sigmaPs << "," << fit.sigmaErrPs << ","
                << fit.chi2Ndf << "," << fit.fwhmPs << "," << landau.mpv << ","
                << landau.sigma << "," << landau.chi2Ndf << "\n";
        const int x = static_cast<int>(scan.x);
        if (x == 0 || x == 350 || x == 690) representative[x] = std::move(events);
    }
    for (auto& item : representative) {
        DrawDelta(item.second, item.first, outputDir);
        DrawNpe(item.second, item.first, outputDir);
    }

    const std::string oldCenter = "results/scan_hi_2026-06-08/photon_hits_run015.root";
    const std::string newCenter =
        "results/scan_end_wrapped_2026-06-09/resume_missing/photon_hits_run003.root";
    const auto oldEvents = ReadEvents(oldCenter);
    const auto newEvents = ReadEvents(newCenter);
    const double p25 = tbmirror::Quantile(newEvents.npe, 0.25);
    const double p50 = tbmirror::Quantile(newEvents.npe, 0.50);
    std::vector<double> new25, new50;
    for (size_t i = 0; i < newEvents.delta.size(); ++i) {
        if (newEvents.npe[i] >= p25) new25.push_back(newEvents.delta[i]);
        if (newEvents.npe[i] >= p50) new50.push_back(newEvents.delta[i]);
    }
    std::ofstream ablation(std::string(outputDir) + "/ablation_results.csv");
    ablation << "variant,dataset,trigger_chain,fit_strategy,energy_cut,"
                "sigma_single_ps,sigma_single_err_ps,n_eff\n";
    auto write = [&](const std::string& variant, const std::string& dataset,
                     const std::string& fitName, const std::string& cut,
                     double sigma, double error, int n) {
        ablation << variant << "," << dataset
                 << ",SUM4_first_cluster_LE_4PE," << fitName << "," << cut << ","
                 << sigma / tbmirror::kSqrt2 << "," << error / tbmirror::kSqrt2
                 << "," << n << "\n";
    };
    const auto oldStage = FitCore(oldEvents.delta, "old_stage");
    const auto oldMirror = tbmirror::FitPeakSeeded(oldEvents.delta, "old_mirror");
    const auto newStage = FitCore(newEvents.delta, "new_stage");
    const auto newMirror = tbmirror::FitPeakSeeded(newEvents.delta, "new_mirror");
    const auto newMirror25 = tbmirror::FitPeakSeeded(new25, "new_mirror25");
    const auto newMirror50 = tbmirror::FitPeakSeeded(new50, "new_mirror50");
    write("historical_stage1_reproduced", "Top-open old scan", "iterative_MAD_pm2sigma",
          "none", oldStage.sigmaPs, oldStage.sigmaErrPs, oldStage.n);
    write("old_geometry_mirror_fit", "Top-open old scan", "peak_local_4ns", "none",
          oldMirror.sigmaPs, oldMirror.sigmaErrPs, oldMirror.n);
    write("wrapped_geometry_stage1_fit", "End-wrapped scan", "iterative_MAD_pm2sigma",
          "none", newStage.sigmaPs, newStage.sigmaErrPs, newStage.n);
    write("wrapped_geometry_mirror_fit", "End-wrapped scan", "peak_local_4ns", "none",
          newMirror.sigmaPs, newMirror.sigmaErrPs, newMirror.n);
    write("wrapped_mirror_NPE_p25", "End-wrapped scan", "peak_local_4ns", "NPE_ge_p25",
          newMirror25.sigmaPs, newMirror25.sigmaErrPs, newMirror25.n);
    write("wrapped_mirror_NPE_p50", "End-wrapped scan", "peak_local_4ns", "NPE_ge_p50",
          newMirror50.sigmaPs, newMirror50.sigmaErrPs, newMirror50.n);
}
