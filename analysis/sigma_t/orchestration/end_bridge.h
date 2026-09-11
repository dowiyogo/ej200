#include "../upstream/related/420addf/analysis/tb_mirror_sigma_vs_x.C"
#include <stdexcept>

// New orchestration only. Never call either imported plotting/main function.
namespace exec34 {
std::vector<double> analyze_end(const std::string& path, int nGenerated, TFile* output) {
    TFile input(path.c_str(), "READ");
    auto* tree = dynamic_cast<TTree*>(input.Get("sipm_hits"));
    if (input.IsZombie() || !tree) throw std::runtime_error("Unreadable sipm_hits");
    TTreeReader reader(tree);
    TTreeReaderValue<int> eventId(reader, "event_id");
    TTreeReaderValue<int> globalId(reader, "global_id");
    TTreeReaderValue<int> face(reader, "face_type");
    TTreeReaderValue<double> time(reader, "time_ns");
    std::vector<EventData> events(nGenerated);
    std::vector<int> countLeft(nGenerated, 0), countRight(nGenerated, 0);
    while (reader.Next()) {
        if (*eventId < 0 || *eventId >= nGenerated)
            throw std::runtime_error("event_id outside explicit generated range");
        if (!std::isfinite(*time)) throw std::runtime_error("Nonfinite hit time");
        if (*face == 2) continue;
        if (*face != 0 && *face != 1) throw std::runtime_error("Unknown face_type");
        const int group = GroupIndex(*globalId);
        if (group < 0 || group > 3 || (*face == 0) != (group < 2))
            throw std::runtime_error("END map/face mismatch");
        events[*eventId].groups[group].push_back(*time);
        if (*face == 0) ++countLeft[*eventId]; else ++countRight[*eventId];
    }
    output->cd();
    TTree eventTree("end_events", "All generated events, including zero hits");
    int event = 0, npeL = 0, npeR = 0;
    double tL, tR, delta;
    bool accepted;
    eventTree.Branch("event_id", &event);
    eventTree.Branch("npe_left", &npeL); eventTree.Branch("npe_right", &npeR);
    eventTree.Branch("t_left_ns", &tL); eventTree.Branch("t_right_ns", &tR);
    eventTree.Branch("delta_ns", &delta); eventTree.Branch("accepted", &accepted);
    std::vector<double> differences;
    int nonfiniteLeft = 0, nonfiniteRight = 0;
    for (event = 0; event < nGenerated; ++event) {
        const auto& g = events[event].groups;
        tL = Earliest(LeadingEdgeTime(g[0]), LeadingEdgeTime(g[1]));
        tR = Earliest(LeadingEdgeTime(g[2]), LeadingEdgeTime(g[3]));
        accepted = std::isfinite(tL) && std::isfinite(tR);
        delta = accepted ? tL-tR : std::numeric_limits<double>::quiet_NaN();
        npeL = countLeft[event]; npeR = countRight[event];
        if (!std::isfinite(tL)) ++nonfiniteLeft;
        if (!std::isfinite(tR)) ++nonfiniteRight;
        if (accepted) differences.push_back(delta);
        eventTree.Fill();
    }
    eventTree.Write();
    const auto core = FitCore(differences, "exec34_end_core");
    const auto peak = tbmirror::FitPeakSeeded(differences, "exec34_end_peak");
    // Persist precisely the histogram definitions used by the two primitives.
    if (!differences.empty()) {
        const double median = Median(differences);
        const double coreWidth = std::max(RobustSigma(differences), 1e-4);
        const double peakWidth = std::max(RobustSigma(differences), .02);
        TH1D hCore("end_core_hist", "", 100, median-8*coreWidth, median+8*coreWidth);
        TH1D hPeak("end_peak_hist", "", 200, median-8*peakWidth, median+8*peakWidth);
        for (double t: differences) {hCore.Fill(t); hPeak.Fill(t);}
        hCore.Write(); hPeak.Write();
    }
    return {core.sigmaPs, core.sigmaErrPs, core.chi2Ndf, double(core.n),
            double(core.usedFit), core.meanNs, core.rmsPs,
            peak.sigmaPs, peak.sigmaErrPs, peak.chi2Ndf, double(peak.n),
            peak.meanNs, peak.fwhmPs, double(nonfiniteLeft), double(nonfiniteRight)};
}
}
