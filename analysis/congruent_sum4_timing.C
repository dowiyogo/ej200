// Metricas de timing congruentes con el test beam: SUM4 + leading-edge.
//
// Etapa 1 intrinseca: usa los tiempos de PE ya simulados sin jitter y no aplica
// time-walk ni cortes ToT. Todos los parametros electronicos son hooks
// provisionales que deben reemplazarse con el repo/analisis de Constanza.
//
// End:
//   clusters SUM4 provisionales = {0..3}, {4..7}, {8..11}, {12..15}
//   t_L/R = primer cluster del extremo que cruza el umbral
//   deltaT = t_L - t_R; sigma_single = sigma(deltaT)/sqrt(2)
//
// Top (solo prediccion de simulacion; no estuvo en el test beam de abril):
//   grupos espejados A={16..50}, B={51..85}
//   deltaT = t_A - t_B; sigma_single = sigma(deltaT)/sqrt(2)
//
// Uso:
//   root -l -b -q 'analysis/congruent_sum4_timing.C("results/scan_hi_2026-06-08","results/congruent_metric_2026-06-08",10000)'

#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {
// HOOK_SPR: respuesta 1-PE provisional, diferencia de exponenciales normalizada.
constexpr double kSprRiseNs = 0.5;
constexpr double kSprFallNs = 5.0;
// HOOK_THR: umbral leading-edge provisional en amplitud equivalente a PE.
constexpr double kThresholdPe = 4.0;
// HOOK_LSB: metadata solamente; no se cuantizan tiempos en etapa 1.
constexpr double kTimeLsbPs = 24.0;
constexpr double kSqrt2 = 1.4142135623730951;

struct EventData {
    std::array<std::vector<double>, 6> groups;
};

struct FitResult {
    double meanNs = 0.0;
    double sigmaPs = 0.0;
    double sigmaErrPs = 0.0;
    double chi2Ndf = 0.0;
    double rmsPs = 0.0;
    int n = 0;
    bool usedFit = false;
};

struct Point {
    double x = 0.0;
    FitResult endDelta;
    FitResult topDelta;
    int endTriggered = 0;
    int topTriggered = 0;
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
        median = 0.5 * (median + *std::max_element(values.begin(), values.begin() + middle));
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
        result.meanNs = mean;
        result.sigmaPs = 1000.0 * sigma;
        result.sigmaErrPs = 1000.0 * gaus.GetParError(2);
        result.chi2Ndf = gaus.GetNDF() > 0 ? gaus.GetChisquare() / gaus.GetNDF() : 0.0;
        result.usedFit = true;
    } else {
        result.meanNs = Mean(times);
        result.sigmaPs = result.rmsPs;
        result.sigmaErrPs = result.rmsPs / std::sqrt(2.0 * (times.size() - 1));
        result.chi2Ndf = -1.0;
    }
    return result;
}

double SprPeakTime() {
    return kSprRiseNs * kSprFallNs / (kSprFallNs - kSprRiseNs)
           * std::log(kSprFallNs / kSprRiseNs);
}

double SprNorm() {
    const double peak = SprPeakTime();
    return 1.0 / (std::exp(-peak / kSprFallNs) - std::exp(-peak / kSprRiseNs));
}

double Pulse(double slowState, double fastState, double dt) {
    return SprNorm() * (slowState * std::exp(-dt / kSprFallNs)
                        - fastState * std::exp(-dt / kSprRiseNs));
}

double LeadingEdgeTime(std::vector<double> arrivals) {
    if (arrivals.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(arrivals.begin(), arrivals.end());
    double slowState = 0.0;
    double fastState = 0.0;
    size_t index = 0;
    double current = arrivals.front();

    while (index < arrivals.size()) {
        current = arrivals[index];
        size_t nextIndex = index;
        while (nextIndex < arrivals.size() && arrivals[nextIndex] == current) {
            slowState += 1.0;
            fastState += 1.0;
            ++nextIndex;
        }

        const double interval = nextIndex < arrivals.size()
            ? arrivals[nextIndex] - current
            : std::numeric_limits<double>::infinity();
        const double derivative0 = fastState / kSprRiseNs - slowState / kSprFallNs;
        if (derivative0 > 0.0) {
            const double peakDt = std::log((fastState * kSprFallNs)
                                           / (slowState * kSprRiseNs))
                                  / (1.0 / kSprRiseNs - 1.0 / kSprFallNs);
            const double reach = std::min(peakDt, interval);
            if (reach >= 0.0 && Pulse(slowState, fastState, reach) >= kThresholdPe) {
                double low = 0.0;
                double high = reach;
                for (int iteration = 0; iteration < 60; ++iteration) {
                    const double middle = 0.5 * (low + high);
                    if (Pulse(slowState, fastState, middle) >= kThresholdPe) high = middle;
                    else low = middle;
                }
                return current + high;
            }
        }

        if (nextIndex >= arrivals.size()) break;
        slowState *= std::exp(-interval / kSprFallNs);
        fastState *= std::exp(-interval / kSprRiseNs);
        index = nextIndex;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

int GroupIndex(int globalId) {
    // HOOK_MAP provisional.
    if (globalId >= 0 && globalId <= 3) return 0;
    if (globalId >= 4 && globalId <= 7) return 1;
    if (globalId >= 8 && globalId <= 11) return 2;
    if (globalId >= 12 && globalId <= 15) return 3;
    // Top is a simulation-only trigger-free split, not a test-beam mapping.
    if (globalId >= 16 && globalId <= 50) return 4;
    if (globalId >= 51 && globalId <= 85) return 5;
    return -1;
}

double Earliest(double first, double second) {
    // HOOK_ENDRED provisional: self-trigger = first cluster to cross.
    if (!std::isfinite(first)) return second;
    if (!std::isfinite(second)) return first;
    return std::min(first, second);
}

void Style(TGraphErrors& graph, int color, int marker) {
    graph.SetLineColor(color);
    graph.SetMarkerColor(color);
    graph.SetLineWidth(2);
    graph.SetMarkerStyle(marker);
    graph.SetMarkerSize(0.9);
}

void SaveResolutionPlot(const std::vector<Point>& points, const std::string& outputDir) {
    TGraphErrors end, top;
    int iEnd = 0;
    int iTop = 0;
    for (const auto& point : points) {
        if (point.endDelta.n >= 20) {
            end.SetPoint(iEnd, point.x, point.endDelta.sigmaPs / kSqrt2);
            end.SetPointError(iEnd++, 0.0, point.endDelta.sigmaErrPs / kSqrt2);
        }
        if (point.topDelta.n >= 20) {
            top.SetPoint(iTop, point.x, point.topDelta.sigmaPs / kSqrt2);
            top.SetPointError(iTop++, 0.0, point.topDelta.sigmaErrPs / kSqrt2);
        }
    }
    Style(end, kBlue + 1, 20);
    Style(top, kGreen + 2, 22);
    TCanvas canvas("sum4_resolution_canvas", "SUM4 resolution", 1100, 650);
    canvas.SetGrid();
    TMultiGraph multi;
    multi.SetTitle("Etapa 1 intrinseca: SUM4 + leading-edge;Posicion x [mm];#sigma_{single} = #sigma(#DeltaT)/#sqrt{2} [ps]");
    multi.Add(&end, "LP");
    multi.Add(&top, "LP");
    multi.Draw("A");
    multi.SetMinimum(0.0);
    TLine target(-690.0, 88.0, 690.0, 88.0);
    target.SetLineStyle(2);
    target.SetLineColor(kRed + 1);
    target.SetLineWidth(2);
    target.Draw();
    TLegend legend(0.52, 0.68, 0.90, 0.89);
    legend.AddEntry(&end, "End SUM4: #sigma(t_{L}-t_{R})/#sqrt{2}", "lp");
    legend.AddEntry(&top, "Top split A/B (sim only; no test-beam data)", "lp");
    legend.AddEntry(&target, "Gerardo: 88 ps single (dato; comparacion cualitativa)", "l");
    legend.Draw();
    canvas.SaveAs((outputDir + "/sigma_single_vs_position.pdf").c_str());
    canvas.SaveAs((outputDir + "/sigma_single_vs_position.png").c_str());
}

void SaveTriggerPlot(const std::vector<Point>& points, const std::string& outputDir, int nEvents) {
    TGraphErrors end(points.size()), top(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        end.SetPoint(i, points[i].x, 100.0 * points[i].endTriggered / nEvents);
        top.SetPoint(i, points[i].x, 100.0 * points[i].topTriggered / nEvents);
    }
    Style(end, kBlue + 1, 20);
    Style(top, kGreen + 2, 22);
    TCanvas canvas("sum4_eff_canvas", "SUM4 trigger efficiency", 1100, 650);
    canvas.SetGrid();
    TMultiGraph multi;
    multi.SetTitle("Eficiencia del observable #DeltaT con hooks provisionales;Posicion x [mm];Eventos con ambos lados/mitades [%]");
    multi.Add(&end, "LP");
    multi.Add(&top, "LP");
    multi.Draw("A");
    multi.SetMinimum(0.0);
    TLegend legend(0.57, 0.74, 0.90, 0.89);
    legend.AddEntry(&end, "End: ambos extremos disparan", "lp");
    legend.AddEntry(&top, "Top A/B (sim only)", "lp");
    legend.Draw();
    canvas.SaveAs((outputDir + "/trigger_efficiency_vs_position.pdf").c_str());
    canvas.SaveAs((outputDir + "/trigger_efficiency_vs_position.png").c_str());
}

void WriteHooks(const std::string& outputDir) {
    std::ofstream out(outputDir + "/hooks.txt");
    out << "HOOK_MAP=provisional: End SUM4 {0-3},{4-7},{8-11},{12-15}; "
           "Top simulation-only split {16-50},{51-85}\n"
        << "HOOK_SPR=provisional: normalized double exponential, tau_r="
        << kSprRiseNs << " ns, tau_f=" << kSprFallNs << " ns, peak amplitude=1 PE\n"
        << "HOOK_THR=provisional: absolute leading-edge threshold=" << kThresholdPe
        << " PE-equivalent amplitude\n"
        << "HOOK_ENDRED=provisional: first of two SUM4 clusters to cross per end\n"
        << "HOOK_WALK=not applied in stage 1\n"
        << "HOOK_TOTCUT=not applied in stage 1; ToT LSB unknown\n"
        << "HOOK_FIT=provisional: 100 bins over median +/-8*MAD-sigma; "
           "4 iterative Gaussian core fits in +/-2 sigma\n"
        << "HOOK_LSB=metadata only: approximately " << kTimeLsbPs
        << " ps/time-LSB; no quantization applied\n";
}
} // namespace

void congruent_sum4_timing(const char* inputDir = "results/scan_hi_2026-06-08",
                           const char* outputDir = "results/congruent_metric_2026-06-08",
                           int nEvents = 10000) {
    gStyle->SetOptStat(0);
    gSystem->mkdir(outputDir, true);
    std::vector<Point> points;

    for (int run = 0; run < 100; ++run) {
        std::ostringstream fileName;
        fileName << inputDir << "/photon_hits_run"
                 << std::setw(3) << std::setfill('0') << run << ".root";
        if (gSystem->AccessPathName(fileName.str().c_str())) break;
        TFile file(fileName.str().c_str(), "READ");
        auto* tree = dynamic_cast<TTree*>(file.Get("sipm_hits"));
        if (file.IsZombie() || tree == nullptr) break;

        TTreeReader reader(tree);
        TTreeReaderValue<int> eventId(reader, "event_id");
        TTreeReaderValue<int> globalId(reader, "global_id");
        TTreeReaderValue<double> time(reader, "time_ns");
        TTreeReaderValue<double> gunX(reader, "gun_x_mm");
        std::vector<EventData> events(nEvents);
        double position = 0.0;
        while (reader.Next()) {
            if (*eventId < 0 || *eventId >= nEvents) continue;
            position = *gunX;
            const int group = GroupIndex(*globalId);
            if (group >= 0) events[*eventId].groups[group].push_back(*time);
        }

        std::vector<double> endDelta, topDelta;
        endDelta.reserve(nEvents);
        topDelta.reserve(nEvents);
        for (auto& event : events) {
            std::array<double, 6> trigger;
            for (size_t group = 0; group < trigger.size(); ++group) {
                trigger[group] = LeadingEdgeTime(std::move(event.groups[group]));
            }
            const double left = Earliest(trigger[0], trigger[1]);
            const double right = Earliest(trigger[2], trigger[3]);
            if (std::isfinite(left) && std::isfinite(right)) endDelta.push_back(left - right);
            if (std::isfinite(trigger[4]) && std::isfinite(trigger[5])) {
                topDelta.push_back(trigger[4] - trigger[5]);
            }
        }

        Point point;
        point.x = position;
        point.endTriggered = static_cast<int>(endDelta.size());
        point.topTriggered = static_cast<int>(topDelta.size());
        point.endDelta = FitCore(endDelta, "sum4_end_delta_" + std::to_string(run));
        point.topDelta = FitCore(topDelta, "sum_top_delta_" + std::to_string(run));
        points.push_back(point);
        std::cout << "Run " << run << ": x=" << position
                  << " mm, N_eff End/Top=" << point.endTriggered << "/"
                  << point.topTriggered << ", sigma_single End/Top="
                  << point.endDelta.sigmaPs / kSqrt2 << "/"
                  << point.topDelta.sigmaPs / kSqrt2 << " ps\n";
    }

    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) { return a.x < b.x; });
    std::ofstream csv(std::string(outputDir) + "/summary.csv");
    csv << "x_mm,n_events,n_eff_end,n_eff_top_sim,trigger_eff_end_percent,trigger_eff_top_sim_percent,"
           "mean_delta_end_ns,sigma_lr_end_ps,sigma_lr_end_err_ps,sigma_single_end_ps,"
           "sigma_single_end_err_ps,chi2_ndf_end,rms_single_end_ps,fit_used_end,"
           "mean_delta_top_sim_ns,sigma_ab_top_sim_ps,sigma_ab_top_sim_err_ps,"
           "sigma_single_top_sim_ps,sigma_single_top_sim_err_ps,chi2_ndf_top_sim,"
           "rms_single_top_sim_ps,fit_used_top_sim\n";
    csv << std::fixed << std::setprecision(6);
    for (const auto& p : points) {
        csv << p.x << "," << nEvents << "," << p.endTriggered << "," << p.topTriggered << ","
            << 100.0 * p.endTriggered / nEvents << "," << 100.0 * p.topTriggered / nEvents << ","
            << p.endDelta.meanNs << "," << p.endDelta.sigmaPs << "," << p.endDelta.sigmaErrPs << ","
            << p.endDelta.sigmaPs / kSqrt2 << "," << p.endDelta.sigmaErrPs / kSqrt2 << ","
            << p.endDelta.chi2Ndf << "," << p.endDelta.rmsPs / kSqrt2 << ","
            << (p.endDelta.usedFit ? "true" : "false") << ","
            << p.topDelta.meanNs << "," << p.topDelta.sigmaPs << "," << p.topDelta.sigmaErrPs << ","
            << p.topDelta.sigmaPs / kSqrt2 << "," << p.topDelta.sigmaErrPs / kSqrt2 << ","
            << p.topDelta.chi2Ndf << "," << p.topDelta.rmsPs / kSqrt2 << ","
            << (p.topDelta.usedFit ? "true" : "false") << "\n";
    }

    SaveResolutionPlot(points, outputDir);
    SaveTriggerPlot(points, outputDir, nEvents);
    WriteHooks(outputDir);
    std::cout << "Analisis SUM4 congruente completado: " << points.size()
              << " posiciones; salida en " << outputDir << "\n";
}
