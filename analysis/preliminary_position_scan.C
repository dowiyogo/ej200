// Analisis preliminar del scan longitudinal EJ-204.
//
// Estimador temporal:
//   - End: promedio de los primeros fotoelectrones de ambos extremos,
//          t_end = (t_L + t_R) / 2.
//   - Top: primer fotoelectron detectado en IDs 16-85.
// El run debe usar /sipm/jitterSigma 0 ns; por tanto sigma_t es INTRINSECO.
//
// Eficiencia = fotones detectados / fotones de centelleo generados, usando
// los contadores master impresos por RunAction al final de cada beamOn.
//
// Uso:
//   root -l -b -q 'analysis/preliminary_position_scan.C("results/scan_2026-06-08",200)'

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
struct EventData {
    int nLeft = 0;
    int nRight = 0;
    int nTop = 0;
    double firstLeft = std::numeric_limits<double>::infinity();
    double firstRight = std::numeric_limits<double>::infinity();
    double firstTop = std::numeric_limits<double>::infinity();
};

struct RunCounters {
    long long scint = 0;
    long long detected = 0;
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
    double sigmaEndPs = 0.0;
    double sigmaEndErrPs = 0.0;
    double sigmaTopPs = 0.0;
    double sigmaTopErrPs = 0.0;
    int nEndTiming = 0;
    int nTopTiming = 0;
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

std::pair<double, double> FitSigmaPs(const std::vector<double>& times,
                                    const std::string& name) {
    if (times.size() < 10) return {0.0, 0.0};
    auto sorted = times;
    std::sort(sorted.begin(), sorted.end());
    const double lo0 = sorted[static_cast<size_t>(0.01 * (sorted.size() - 1))];
    const double hi0 = sorted[static_cast<size_t>(0.99 * (sorted.size() - 1))];
    const double margin = std::max(0.5, 0.5 * (hi0 - lo0));
    const double lo = lo0 - margin;
    const double hi = hi0 + margin;

    TH1D hist(name.c_str(), "", 35, lo, hi);
    for (double time : times) hist.Fill(time);

    TF1 gaus((name + "_gaus").c_str(), "gaus", lo, hi);
    gaus.SetParameters(hist.GetMaximum(), Mean(times), std::max(StdDev(times), 0.01));
    gaus.SetParLimits(2, 1e-5, hi - lo);
    const int status = hist.Fit(&gaus, "QNR");
    if (status == 0 && gaus.GetParameter(2) > 0.0) {
        return {1000.0 * std::abs(gaus.GetParameter(2)),
                1000.0 * gaus.GetParError(2)};
    }
    const double sigma = StdDev(times);
    return {1000.0 * sigma, 1000.0 * sigma / std::sqrt(times.size())};
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
            if (currentRun >= 0 && events == nEvents) {
                counters[currentRun] = {scint, detected};
            }
        }
    }
    return counters;
}

void Style(TGraphErrors& graph, int color, int marker) {
    graph.SetLineColor(color);
    graph.SetMarkerColor(color);
    graph.SetLineWidth(2);
    graph.SetMarkerStyle(marker);
    graph.SetMarkerSize(1.1);
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

    TCanvas canvas("yield_canvas", "yield", 1000, 600);
    canvas.SetGrid();
    TMultiGraph multi;
    multi.SetTitle("EJ-204: yield detectado vs posicion;Posicion x [mm];Fotoelectrones detectados / evento");
    multi.Add(&left, "LP");
    multi.Add(&right, "LP");
    multi.Add(&top, "LP");
    multi.Draw("A");
    TLegend legend(0.68, 0.72, 0.90, 0.89);
    legend.AddEntry(&left, "End-left (IDs 0-7)", "lp");
    legend.AddEntry(&right, "End-right (IDs 8-15)", "lp");
    legend.AddEntry(&top, "Top (IDs 16-85)", "lp");
    legend.Draw();
    canvas.SaveAs((base + "/yield_vs_position.pdf").c_str());
    canvas.SaveAs((base + "/yield_vs_position.png").c_str());
}

void SaveTimingPlot(const std::vector<Point>& points, const std::string& base) {
    TGraphErrors end(points.size()), top(points.size());
    for (size_t i = 0; i < points.size(); ++i) {
        end.SetPoint(i, points[i].x, points[i].sigmaEndPs);
        top.SetPoint(i, points[i].x, points[i].sigmaTopPs);
        end.SetPointError(i, 0.0, points[i].sigmaEndErrPs);
        top.SetPointError(i, 0.0, points[i].sigmaTopErrPs);
    }
    Style(end, kViolet + 1, 20);
    Style(top, kGreen + 2, 22);

    TCanvas canvas("timing_canvas", "timing", 1000, 600);
    canvas.SetGrid();
    TMultiGraph multi;
    multi.SetTitle("EJ-204: resolucion temporal intrinseca (FPT);Posicion x [mm];#sigma_{t} intrinseco [ps]");
    multi.Add(&end, "LP");
    multi.Add(&top, "LP");
    multi.Draw("A");
    multi.SetMinimum(0.0);
    TLegend legend(0.55, 0.73, 0.90, 0.89);
    legend.AddEntry(&end, "End: (FPT_{L}+FPT_{R})/2", "lp");
    legend.AddEntry(&top, "Top: primer fotoelectron", "lp");
    legend.Draw();
    canvas.SaveAs((base + "/sigma_intrinsic_vs_position.pdf").c_str());
    canvas.SaveAs((base + "/sigma_intrinsic_vs_position.png").c_str());
}

void SaveEfficiencyPlot(const std::vector<Point>& points, const std::string& base) {
    TGraphErrors graph(points.size());
    for (size_t i = 0; i < points.size(); ++i) graph.SetPoint(i, points[i].x, points[i].efficiency);
    Style(graph, kOrange + 7, 20);
    TCanvas canvas("eff_canvas", "efficiency", 1000, 600);
    canvas.SetGrid();
    graph.SetTitle("EJ-204: eficiencia de deteccion;Posicion x [mm];Fotones detectados / fotones de centelleo [%]");
    graph.Draw("ALP");
    graph.SetMinimum(0.0);
    canvas.SaveAs((base + "/efficiency_vs_position.pdf").c_str());
    canvas.SaveAs((base + "/efficiency_vs_position.png").c_str());
}
} // namespace

void preliminary_position_scan(const char* resultDir = "results/scan_2026-06-08",
                               int nEvents = 200) {
    gStyle->SetOptStat(0);
    const std::string base(resultDir);
    const auto counters = ParseMasterCounters(base + "/run.log", nEvents);
    std::vector<Point> points;

    for (int run = 0; run < 100; ++run) {
        std::ostringstream fileName;
        fileName << base << "/photon_hits_run" << std::setw(3) << std::setfill('0') << run << ".root";
        if (gSystem->AccessPathName(fileName.str().c_str())) {
            if (run == 0) std::cerr << "No se encontro " << fileName.str() << "\n";
            break;
        }
        TFile file(fileName.str().c_str(), "READ");
        if (file.IsZombie()) {
            break;
        }
        auto* tree = dynamic_cast<TTree*>(file.Get("sipm_hits"));
        if (tree == nullptr) break;

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
            if (*face == 0) {
                ++event.nLeft;
                event.firstLeft = std::min(event.firstLeft, *time);
            } else if (*face == 1) {
                ++event.nRight;
                event.firstRight = std::min(event.firstRight, *time);
            } else if (*face == 2) {
                ++event.nTop;
                event.firstTop = std::min(event.firstTop, *time);
            }
        }

        std::vector<double> leftCounts, rightCounts, topCounts, endTimes, topTimes;
        for (const auto& event : events) {
            leftCounts.push_back(event.nLeft);
            rightCounts.push_back(event.nRight);
            topCounts.push_back(event.nTop);
            if (std::isfinite(event.firstLeft) && std::isfinite(event.firstRight)) {
                endTimes.push_back(0.5 * (event.firstLeft + event.firstRight));
            }
            if (std::isfinite(event.firstTop)) topTimes.push_back(event.firstTop);
        }

        Point point;
        point.run = run;
        point.x = position;
        point.yieldLeft = Mean(leftCounts);
        point.yieldRight = Mean(rightCounts);
        point.yieldTop = Mean(topCounts);
        point.yieldLeftErr = StdDev(leftCounts) / std::sqrt(nEvents);
        point.yieldRightErr = StdDev(rightCounts) / std::sqrt(nEvents);
        point.yieldTopErr = StdDev(topCounts) / std::sqrt(nEvents);
        const auto endFit = FitSigmaPs(endTimes, "end_" + std::to_string(run));
        const auto topFit = FitSigmaPs(topTimes, "top_" + std::to_string(run));
        point.sigmaEndPs = endFit.first;
        point.sigmaEndErrPs = endFit.second;
        point.sigmaTopPs = topFit.first;
        point.sigmaTopErrPs = topFit.second;
        point.nEndTiming = endTimes.size();
        point.nTopTiming = topTimes.size();
        const auto counter = counters.find(run);
        if (counter != counters.end() && counter->second.scint > 0) {
            point.efficiency = 100.0 * counter->second.detected / counter->second.scint;
        }
        points.push_back(point);
    }

    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) { return a.x < b.x; });
    std::ofstream csv(base + "/summary.csv");
    csv << "x_mm,n_events,yield_end_left,yield_end_right,yield_top,efficiency_percent,"
           "sigma_end_intrinsic_ps,sigma_end_err_ps,sigma_top_intrinsic_ps,sigma_top_err_ps,"
           "n_end_timing,n_top_timing\n";
    csv << std::fixed << std::setprecision(4);
    for (const auto& point : points) {
        csv << point.x << "," << nEvents << ","
            << point.yieldLeft << "," << point.yieldRight << "," << point.yieldTop << ","
            << point.efficiency << "," << point.sigmaEndPs << "," << point.sigmaEndErrPs << ","
            << point.sigmaTopPs << "," << point.sigmaTopErrPs << ","
            << point.nEndTiming << "," << point.nTopTiming << "\n";
    }

    SaveYieldPlot(points, base);
    SaveTimingPlot(points, base);
    SaveEfficiencyPlot(points, base);
    std::cout << "Analisis completado: " << points.size() << " posiciones; salida en " << base << "\n";
}
