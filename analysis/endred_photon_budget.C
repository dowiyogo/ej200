// Diagnostico de HOOK_ENDRED usando los ROOT existentes del scan de 31 puntos.
//
// Compara por extremo:
//   first_sum4: primer cluster SUM4 que cruza el umbral (hook actual)
//   sum8:       los 8 SiPM se agregan en un solo pulso
//   mean_sum4:  promedio de los tiempos de los dos clusters SUM4
//
// Etapa intrinseca: no aplica jitter, walk ni cortes ToT. HOOK_SUM_WIDTH
// selecciona la curva "best" (4 o 8), sin cambiar las tres comparaciones.
//
// Uso:
// root -l -b -q 'analysis/endred_photon_budget.C("results/scan_hi_2026-06-08","results/diag_2026-06-09",10000,8)'

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
constexpr double kSprRiseNs = 0.5;   // HOOK_SPR provisional
constexpr double kSprFallNs = 5.0;
constexpr double kThresholdPe = 4.0; // HOOK_THR provisional
constexpr double kSqrt2 = 1.4142135623730951;

struct EventData {
    std::array<std::vector<double>, 16> channel;
};

struct FitResult {
    double sigmaPs = 0.0;
    double sigmaErrPs = 0.0;
    double chi2Ndf = 0.0;
    double rmsPs = 0.0;
    int n = 0;
};

struct ModeResult {
    FitResult delta;
    std::string name;
};

struct Point {
    double x = 0.0;
    double meanLeftPe = 0.0;
    double meanRightPe = 0.0;
    std::array<ModeResult, 3> modes;
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

FitResult FitCore(const std::vector<double>& values, const std::string& name) {
    FitResult result;
    result.n = static_cast<int>(values.size());
    result.rmsPs = 1000.0 * StdDev(values);
    if (values.size() < 20) return result;

    const double center = Median(values);
    const double robustSigma = std::max(RobustSigma(values), 1e-4);
    const double lo = center - 8.0 * robustSigma;
    const double hi = center + 8.0 * robustSigma;
    TH1D hist(name.c_str(), "", 100, lo, hi);
    hist.SetDirectory(nullptr);
    for (double value : values) hist.Fill(value);

    double mean = center;
    double sigma = robustSigma;
    TF1 gaus((name + "_gaus").c_str(), "gaus", lo, hi);
    int status = 1;
    for (int iteration = 0; iteration < 4; ++iteration) {
        gaus.SetRange(std::max(lo, mean - 2.0 * sigma), std::min(hi, mean + 2.0 * sigma));
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
        result.sigmaErrPs = result.rmsPs / std::sqrt(2.0 * (values.size() - 1));
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
    while (index < arrivals.size()) {
        const double current = arrivals[index];
        size_t next = index;
        while (next < arrivals.size() && arrivals[next] == current) {
            slowState += 1.0;
            fastState += 1.0;
            ++next;
        }
        const double interval = next < arrivals.size()
            ? arrivals[next] - current : std::numeric_limits<double>::infinity();
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
        if (next >= arrivals.size()) break;
        slowState *= std::exp(-interval / kSprFallNs);
        fastState *= std::exp(-interval / kSprRiseNs);
        index = next;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

std::vector<double> Merge(const EventData& event, int first, int last) {
    std::vector<double> merged;
    for (int id = first; id <= last; ++id) {
        merged.insert(merged.end(), event.channel[id].begin(), event.channel[id].end());
    }
    return merged;
}

double Earliest(double a, double b) {
    if (!std::isfinite(a)) return b;
    if (!std::isfinite(b)) return a;
    return std::min(a, b);
}

double AverageBoth(double a, double b) {
    if (!std::isfinite(a) || !std::isfinite(b)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return 0.5 * (a + b);
}

void Style(TGraphErrors& graph, int color, int marker) {
    graph.SetLineColor(color);
    graph.SetMarkerColor(color);
    graph.SetLineWidth(2);
    graph.SetMarkerStyle(marker);
    graph.SetMarkerSize(0.85);
}

void SavePlot(const std::vector<Point>& points, const std::string& outputDir, int bestSumWidth) {
    TGraphErrors first, sum8, combined, bestVsLight;
    for (size_t i = 0; i < points.size(); ++i) {
        const auto& point = points[i];
        const std::array<TGraphErrors*, 3> graphs = {&first, &sum8, &combined};
        for (size_t mode = 0; mode < graphs.size(); ++mode) {
            const auto& fit = point.modes[mode].delta;
            graphs[mode]->SetPoint(i, point.x, fit.sigmaPs / kSqrt2);
            graphs[mode]->SetPointError(i, 0.0, fit.sigmaErrPs / kSqrt2);
        }
        const size_t bestMode = bestSumWidth == 8 ? 1 : 0;
        const auto& best = point.modes[bestMode].delta;
        bestVsLight.SetPoint(i, point.meanLeftPe + point.meanRightPe, best.sigmaPs / kSqrt2);
        bestVsLight.SetPointError(i, 0.0, best.sigmaErrPs / kSqrt2);
    }
    Style(first, kBlue + 1, 20);
    Style(sum8, kRed + 1, 21);
    Style(combined, kGreen + 2, 22);
    TCanvas canvas("endred_canvas", "ENDRED comparison", 1150, 700);
    canvas.SetGrid();
    TMultiGraph multi;
    multi.SetTitle("Dependencia de HOOK_ENDRED;Posicion x [mm];#sigma_{single}(End) intrinseco [ps]");
    multi.Add(&first, "LP");
    multi.Add(&sum8, "LP");
    multi.Add(&combined, "LP");
    multi.Draw("A");
    multi.SetMinimum(0.0);
    TLegend legend(0.58, 0.70, 0.89, 0.89);
    legend.AddEntry(&first, "Primer SUM4 (actual)", "lp");
    legend.AddEntry(&sum8, "SUM8: un pulso/extremo", "lp");
    legend.AddEntry(&combined, "Promedio de dos tiempos SUM4", "lp");
    legend.Draw();
    canvas.SaveAs((outputDir + "/endred_comparison_vs_position.pdf").c_str());
    canvas.SaveAs((outputDir + "/endred_comparison_vs_position.png").c_str());

    Style(bestVsLight, kRed + 1, 21);
    TCanvas lightCanvas("endred_light_canvas", "Best ENDRED vs light", 1000, 650);
    lightCanvas.SetGrid();
    bestVsLight.SetTitle(("HOOK_SUM_WIDTH=" + std::to_string(bestSumWidth)
        + ";N_{PE,End-L}+N_{PE,End-R} por evento;#sigma_{single}(End) intrinseco [ps]").c_str());
    bestVsLight.Draw("AP");
    lightCanvas.SaveAs((outputDir + "/best_endred_vs_light.pdf").c_str());
    lightCanvas.SaveAs((outputDir + "/best_endred_vs_light.png").c_str());
}
} // namespace

void endred_photon_budget(const char* inputDir = "results/scan_hi_2026-06-08",
                          const char* outputDir = "results/diag_2026-06-09",
                          int nEvents = 10000,
                          int bestSumWidth = 8) {
    if (bestSumWidth != 4 && bestSumWidth != 8) {
        std::cerr << "HOOK_SUM_WIDTH must be 4 or 8\n";
        return;
    }
    gStyle->SetOptStat(0);
    gSystem->mkdir(outputDir, true);
    std::vector<Point> points;
    const std::array<std::string, 3> names = {"first_sum4", "sum8", "mean_sum4"};

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
        long long leftPe = 0;
        long long rightPe = 0;
        while (reader.Next()) {
            if (*eventId < 0 || *eventId >= nEvents || *globalId < 0 || *globalId > 15) continue;
            position = *gunX;
            events[*eventId].channel[*globalId].push_back(*time);
            if (*globalId < 8) ++leftPe;
            else ++rightPe;
        }

        std::array<std::vector<double>, 3> deltas;
        for (const auto& event : events) {
            const double l0 = LeadingEdgeTime(Merge(event, 0, 3));
            const double l1 = LeadingEdgeTime(Merge(event, 4, 7));
            const double r0 = LeadingEdgeTime(Merge(event, 8, 11));
            const double r1 = LeadingEdgeTime(Merge(event, 12, 15));
            const std::array<double, 3> left = {
                Earliest(l0, l1), LeadingEdgeTime(Merge(event, 0, 7)), AverageBoth(l0, l1)
            };
            const std::array<double, 3> right = {
                Earliest(r0, r1), LeadingEdgeTime(Merge(event, 8, 15)), AverageBoth(r0, r1)
            };
            for (size_t mode = 0; mode < deltas.size(); ++mode) {
                if (std::isfinite(left[mode]) && std::isfinite(right[mode])) {
                    deltas[mode].push_back(left[mode] - right[mode]);
                }
            }
        }

        Point point;
        point.x = position;
        point.meanLeftPe = static_cast<double>(leftPe) / nEvents;
        point.meanRightPe = static_cast<double>(rightPe) / nEvents;
        for (size_t mode = 0; mode < names.size(); ++mode) {
            point.modes[mode].name = names[mode];
            point.modes[mode].delta = FitCore(
                deltas[mode], names[mode] + "_run_" + std::to_string(run));
        }
        points.push_back(point);
        std::cout << "x=" << position << " mm, End PE L/R=" << point.meanLeftPe
                  << "/" << point.meanRightPe << ", sigma_single first/sum8/mean="
                  << point.modes[0].delta.sigmaPs / kSqrt2 << "/"
                  << point.modes[1].delta.sigmaPs / kSqrt2 << "/"
                  << point.modes[2].delta.sigmaPs / kSqrt2 << " ps\n";
    }

    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) { return a.x < b.x; });
    std::ofstream csv(std::string(outputDir) + "/endred_summary.csv");
    csv << "x_mm,mean_end_left_pe,mean_end_right_pe,mode,n_eff,sigma_lr_ps,"
           "sigma_lr_err_ps,sigma_single_ps,sigma_single_err_ps,chi2_ndf,rms_single_ps\n";
    csv << std::fixed << std::setprecision(6);
    for (const auto& point : points) {
        for (const auto& mode : point.modes) {
            const auto& fit = mode.delta;
            csv << point.x << "," << point.meanLeftPe << "," << point.meanRightPe << ","
                << mode.name << "," << fit.n << "," << fit.sigmaPs << "," << fit.sigmaErrPs
                << "," << fit.sigmaPs / kSqrt2 << "," << fit.sigmaErrPs / kSqrt2 << ","
                << fit.chi2Ndf << "," << fit.rmsPs / kSqrt2 << "\n";
        }
    }
    SavePlot(points, outputDir, bestSumWidth);
    std::ofstream hooks(std::string(outputDir) + "/endred_hooks.txt");
    hooks << "HOOK_SUM_WIDTH=" << bestSumWidth << " (accepted values: 4 or 8)\n"
          << "HOOK_MAP=provisional: End clusters {0-3},{4-7},{8-11},{12-15}\n"
          << "HOOK_SPR=provisional double exponential tau_r=" << kSprRiseNs
          << " ns, tau_f=" << kSprFallNs << " ns\n"
          << "HOOK_THR=provisional absolute leading-edge threshold=" << kThresholdPe << " PE\n"
          << "HOOK_ENDRED compared: first SUM4, one SUM8 pulse, mean of two SUM4 times\n"
          << "Final comparison with 88 ps requires Gerardo measurement position and SUM width.\n";
    std::cout << "ENDRED diagnostic complete: " << points.size() << " positions\n";
}
