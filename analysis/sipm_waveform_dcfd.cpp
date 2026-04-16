#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TSystem.h"

namespace fs = std::filesystem;

namespace {

struct Config {
    std::vector<std::string> inputs;
    std::string treeName = "sipm_hits";
    std::string timeBranch = "time_raw_ns";
    std::string face = "0";
    std::string histOut = "sipm_dcfd_resolution.png";
    std::string csvOut = "sipm_dcfd_times.csv";
    double fraction = 0.14;
    double dtPs = 10.0;
    double transitSigmaPs = 200.0;
    double tauRNs = 2.0;
    double tauFNs = 55.0;
    double pulseSigmaPe = 0.1;
    double electronicsSigmaPs = 0.0;
    double tailWindowFactor = 8.0;
    long long maxEvents = -1;
    unsigned int seed = 12345;
    bool savePlot = true;
};

struct HitRecord {
    int eventId = 0;
    int faceType = 0;
    double t0Ns = 0.0;
    double gunXmm = 0.0;
    int treeNumber = -1;
    std::string sourceFile;
};

struct EventResult {
    std::string eventUid;
    std::string sourceFile;
    int eventId = 0;
    int faceType = 0;
    int nPe = 0;
    double tDcfdNs = std::numeric_limits<double>::quiet_NaN();
    double amplitudePe = 0.0;
    double tPeakNs = 0.0;
    double gunXmm = 0.0;
};

struct FitResult {
    double meanNs = std::numeric_limits<double>::quiet_NaN();
    double sigmaNs = std::numeric_limits<double>::quiet_NaN();
    double sigmaErrNs = std::numeric_limits<double>::quiet_NaN();
};

void PrintUsage(const char* prog) {
    std::cout
        << "Uso:\n"
        << "  " << prog << " build/Data/*.root [opciones]\n\n"
        << "Opciones:\n"
        << "  --tree NAME                Nombre del TTree (default: sipm_hits)\n"
        << "  --time-branch NAME         Branch con t0 base (default: time_raw_ns)\n"
        << "  --face {0|1|2|all}         Cara a analizar (default: 0)\n"
        << "  --fraction X               Fraccion dCFD (default: 0.14)\n"
        << "  --dt-ps X                  Paso temporal del waveform en ps (default: 10)\n"
        << "  --transit-sigma-ps X       Sigma del transit time en ps (default: 200)\n"
        << "  --tau-r-ns X               Tau de subida en ns (default: 2)\n"
        << "  --tau-f-ns X               Tau de bajada en ns (default: 55)\n"
        << "  --pulse-sigma-pe X         Sigma de amplitud SPE en pe (default: 0.1)\n"
        << "  --electronics-sigma-ps X   Jitter electronico extra en ps (default: 0)\n"
        << "  --tail-window-factor X     Longitud del kernel en multiplos de tauF (default: 8)\n"
        << "  --max-events N             Limite de eventos reconstruidos\n"
        << "  --seed N                   Semilla RNG (default: 12345)\n"
        << "  --hist-out FILE            PNG de histograma (default: sipm_dcfd_resolution.png)\n"
        << "  --csv-out FILE             CSV de tiempos (default: sipm_dcfd_times.csv)\n"
        << "  --no-plot                  No guardar histograma PNG\n";
}

bool ParseDoubleArg(int& i, int argc, char** argv, double& target) {
    if (i + 1 >= argc) return false;
    target = std::stod(argv[++i]);
    return true;
}

bool ParseIntArg(int& i, int argc, char** argv, long long& target) {
    if (i + 1 >= argc) return false;
    target = std::stoll(argv[++i]);
    return true;
}

bool ParseUIntArg(int& i, int argc, char** argv, unsigned int& target) {
    if (i + 1 >= argc) return false;
    target = static_cast<unsigned int>(std::stoul(argv[++i]));
    return true;
}

bool ParseStringArg(int& i, int argc, char** argv, std::string& target) {
    if (i + 1 >= argc) return false;
    target = argv[++i];
    return true;
}

std::optional<Config> ParseArgs(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            PrintUsage(argv[0]);
            return std::nullopt;
        }
        if (arg == "--tree") {
            if (!ParseStringArg(i, argc, argv, cfg.treeName)) return std::nullopt;
        } else if (arg == "--time-branch") {
            if (!ParseStringArg(i, argc, argv, cfg.timeBranch)) return std::nullopt;
        } else if (arg == "--face") {
            if (!ParseStringArg(i, argc, argv, cfg.face)) return std::nullopt;
        } else if (arg == "--fraction") {
            if (!ParseDoubleArg(i, argc, argv, cfg.fraction)) return std::nullopt;
        } else if (arg == "--dt-ps") {
            if (!ParseDoubleArg(i, argc, argv, cfg.dtPs)) return std::nullopt;
        } else if (arg == "--transit-sigma-ps") {
            if (!ParseDoubleArg(i, argc, argv, cfg.transitSigmaPs)) return std::nullopt;
        } else if (arg == "--tau-r-ns") {
            if (!ParseDoubleArg(i, argc, argv, cfg.tauRNs)) return std::nullopt;
        } else if (arg == "--tau-f-ns") {
            if (!ParseDoubleArg(i, argc, argv, cfg.tauFNs)) return std::nullopt;
        } else if (arg == "--pulse-sigma-pe") {
            if (!ParseDoubleArg(i, argc, argv, cfg.pulseSigmaPe)) return std::nullopt;
        } else if (arg == "--electronics-sigma-ps") {
            if (!ParseDoubleArg(i, argc, argv, cfg.electronicsSigmaPs)) return std::nullopt;
        } else if (arg == "--tail-window-factor") {
            if (!ParseDoubleArg(i, argc, argv, cfg.tailWindowFactor)) return std::nullopt;
        } else if (arg == "--max-events") {
            if (!ParseIntArg(i, argc, argv, cfg.maxEvents)) return std::nullopt;
        } else if (arg == "--seed") {
            if (!ParseUIntArg(i, argc, argv, cfg.seed)) return std::nullopt;
        } else if (arg == "--hist-out") {
            if (!ParseStringArg(i, argc, argv, cfg.histOut)) return std::nullopt;
        } else if (arg == "--csv-out") {
            if (!ParseStringArg(i, argc, argv, cfg.csvOut)) return std::nullopt;
        } else if (arg == "--no-plot") {
            cfg.savePlot = false;
        } else if (!arg.empty() && arg.rfind("--", 0) == 0) {
            std::cerr << "Opcion no reconocida: " << arg << "\n";
            return std::nullopt;
        } else {
            cfg.inputs.push_back(arg);
        }
    }
    if (cfg.inputs.empty()) {
        cfg.inputs.push_back("build/Data/*.root");
    }
    return cfg;
}

std::vector<double> BuildKernel(double dtNs, double tauRNs, double tauFNs, double tailWindowFactor) {
    const double kernelSpanNs = std::max(tailWindowFactor * tauFNs, 5.0 * dtNs);
    const int nBins = static_cast<int>(std::ceil(kernelSpanNs / dtNs)) + 1;
    std::vector<double> kernel;
    kernel.reserve(nBins);
    for (int i = 0; i < nBins; ++i) {
        const double t = i * dtNs;
        kernel.push_back((1.0 - std::exp(-t / tauRNs)) * std::exp(-t / tauFNs));
    }
    return kernel;
}

bool FaceMatches(const std::string& requested, int faceType) {
    if (requested == "all") return true;
    return faceType == std::stoi(requested);
}

double ComputeDcfdTime(const std::vector<double>& timeNs, const std::vector<double>& waveform, double fraction) {
    if (waveform.empty()) return std::numeric_limits<double>::quiet_NaN();

    const auto maxIt = std::max_element(waveform.begin(), waveform.end());
    const int peakIdx = static_cast<int>(std::distance(waveform.begin(), maxIt));
    const double peakAmp = *maxIt;
    if (!(peakAmp > 0.0)) return std::numeric_limits<double>::quiet_NaN();

    const double threshold = fraction * peakAmp;
    for (int i = 0; i <= peakIdx; ++i) {
        if (waveform[i] >= threshold) {
            if (i == 0) return timeNs.front();
            const double x0 = timeNs[i - 1];
            const double x1 = timeNs[i];
            const double y0 = waveform[i - 1];
            const double y1 = waveform[i];
            if (std::abs(y1 - y0) < 1e-15) return x1;
            return x0 + (threshold - y0) * (x1 - x0) / (y1 - y0);
        }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

EventResult ProcessEvent(
    const std::vector<HitRecord>& hits,
    const std::vector<double>& kernel,
    const Config& cfg,
    std::mt19937& rng
) {
    EventResult result;
    if (hits.empty()) return result;

    const double dtNs = cfg.dtPs * 1e-3;
    const double transitSigmaNs = cfg.transitSigmaPs * 1e-3;
    const double electronicsSigmaNs = cfg.electronicsSigmaPs * 1e-3;

    std::normal_distribution<double> transitDist(0.0, transitSigmaNs);
    std::normal_distribution<double> amplitudeDist(1.0, cfg.pulseSigmaPe);
    std::normal_distribution<double> electronicsDist(0.0, electronicsSigmaNs);

    std::vector<double> shiftedTimes;
    shiftedTimes.reserve(hits.size());
    std::vector<double> amplitudes;
    amplitudes.reserve(hits.size());

    for (const auto& hit : hits) {
        shiftedTimes.push_back(hit.t0Ns + transitDist(rng));
        amplitudes.push_back(std::max(0.0, amplitudeDist(rng)));
    }

    const auto [minIt, maxIt] = std::minmax_element(shiftedTimes.begin(), shiftedTimes.end());
    const double preNs = std::max(5.0 * transitSigmaNs, dtNs);
    const double startNs = *minIt - preNs;
    const double stopNs = *maxIt + dtNs * static_cast<double>(kernel.size());
    const int nBins = static_cast<int>(std::ceil((stopNs - startNs) / dtNs)) + 1;

    std::vector<double> waveform(static_cast<size_t>(nBins), 0.0);
    std::vector<double> timeAxis(static_cast<size_t>(nBins), 0.0);
    for (int i = 0; i < nBins; ++i) {
        timeAxis[static_cast<size_t>(i)] = startNs + dtNs * static_cast<double>(i);
    }

    for (size_t ih = 0; ih < shiftedTimes.size(); ++ih) {
        int hitBin = static_cast<int>(std::llround((shiftedTimes[ih] - startNs) / dtNs));
        hitBin = std::max(0, std::min(hitBin, nBins - 1));
        for (size_t k = 0; k < kernel.size(); ++k) {
            const int idx = hitBin + static_cast<int>(k);
            if (idx >= nBins) break;
            waveform[static_cast<size_t>(idx)] += amplitudes[ih] * kernel[k];
        }
    }

    if (electronicsSigmaNs > 0.0) {
        for (double& sample : waveform) {
            sample += electronicsDist(rng);
        }
    }

    const double tDcfdNs = ComputeDcfdTime(timeAxis, waveform, cfg.fraction);
    const auto peakIt = std::max_element(waveform.begin(), waveform.end());
    const int peakIdx = static_cast<int>(std::distance(waveform.begin(), peakIt));

    result.sourceFile = hits.front().sourceFile;
    result.eventId = hits.front().eventId;
    result.faceType = hits.front().faceType;
    result.nPe = static_cast<int>(hits.size());
    result.tDcfdNs = tDcfdNs;
    result.amplitudePe = (peakIt != waveform.end()) ? *peakIt : 0.0;
    result.tPeakNs = timeAxis.empty() ? 0.0 : timeAxis[static_cast<size_t>(peakIdx)];
    result.gunXmm = hits.front().gunXmm;

    std::ostringstream uid;
    uid << hits.front().sourceFile << ":" << hits.front().treeNumber << ":"
        << hits.front().eventId << ":" << hits.front().faceType;
    result.eventUid = uid.str();
    return result;
}

void WriteCsv(const std::string& path, const std::vector<EventResult>& results) {
    std::ofstream out(path);
    out << "event_uid,source_file,event_id,face_type,n_pe,t_dcfd_ns,amplitude_pe,t_peak_ns,gun_x_mm\n";
    out << std::fixed << std::setprecision(8);
    for (const auto& row : results) {
        out << row.eventUid << ','
            << row.sourceFile << ','
            << row.eventId << ','
            << row.faceType << ','
            << row.nPe << ','
            << row.tDcfdNs << ','
            << row.amplitudePe << ','
            << row.tPeakNs << ','
            << row.gunXmm << '\n';
    }
}

FitResult FitResolution(const std::vector<EventResult>& results, const Config& cfg) {
    FitResult fit;
    if (results.empty()) return fit;

    std::vector<double> times;
    times.reserve(results.size());
    for (const auto& row : results) {
        if (std::isfinite(row.tDcfdNs)) times.push_back(row.tDcfdNs);
    }
    if (times.empty()) return fit;

    const auto [minIt, maxIt] = std::minmax_element(times.begin(), times.end());
    double lo = *minIt;
    double hi = *maxIt;
    if (!(hi > lo)) {
        lo -= 0.5;
        hi += 0.5;
    } else {
        const double margin = 0.15 * (hi - lo);
        lo -= margin;
        hi += margin;
    }

    TH1D hist(
        "hDcfdTime",
        Form("Tiempo dCFD %.0f%%;t_{dCFD} [ns];Eventos / bin", cfg.fraction * 100.0),
        std::min(80, std::max(25, static_cast<int>(2.0 * std::sqrt(static_cast<double>(times.size()))))),
        lo, hi
    );
    hist.SetDirectory(nullptr);
    for (double t : times) hist.Fill(t);

    fit.meanNs = hist.GetMean();
    fit.sigmaNs = hist.GetRMS();
    if (hist.GetEntries() > 1) {
        fit.sigmaErrNs = fit.sigmaNs / std::sqrt(2.0 * (hist.GetEntries() - 1.0));
    }

    if (hist.GetEntries() >= 5 && fit.sigmaNs > 0.0) {
        TF1 gaus("gausFit", "gaus", fit.meanNs - 2.0 * fit.sigmaNs, fit.meanNs + 2.0 * fit.sigmaNs);
        gaus.SetParameters(hist.GetMaximum(), fit.meanNs, fit.sigmaNs);
        if (hist.Fit(&gaus, "QNR") == 0) {
            fit.meanNs = gaus.GetParameter(1);
            fit.sigmaNs = std::abs(gaus.GetParameter(2));
            fit.sigmaErrNs = gaus.GetParError(2);
        }
    }

    if (cfg.savePlot) {
        gROOT->SetBatch(kTRUE);
        TCanvas canvas("cDcfd", "dCFD timing", 900, 600);
        hist.SetLineColor(kAzure + 2);
        hist.SetFillColorAlpha(kAzure - 9, 0.65);
        hist.Draw("hist");
        if (hist.GetEntries() >= 5 && std::isfinite(fit.sigmaNs) && fit.sigmaNs > 0.0) {
            TF1 drawGaus("drawGaus", "gaus", hist.GetXaxis()->GetXmin(), hist.GetXaxis()->GetXmax());
            drawGaus.SetParameters(hist.GetMaximum(), fit.meanNs, fit.sigmaNs);
            drawGaus.SetLineColor(kRed + 1);
            drawGaus.SetLineWidth(2);
            drawGaus.Draw("same");
        }
        canvas.SaveAs(cfg.histOut.c_str());
    }

    return fit;
}

std::vector<EventResult> AnalyzeChain(TChain& chain, const Config& cfg) {
    int eventId = 0;
    int faceType = 0;
    double timeNs = 0.0;
    double gunXmm = 0.0;

    bool hasTimeBranch = (chain.GetBranch(cfg.timeBranch.c_str()) != nullptr);
    const char* activeTimeBranch = hasTimeBranch ? cfg.timeBranch.c_str() : "time_ns";
    const bool hasGunX = (chain.GetBranch("gun_x_mm") != nullptr);

    chain.SetBranchStatus("*", 0);
    chain.SetBranchStatus("event_id", 1);
    chain.SetBranchStatus("face_type", 1);
    chain.SetBranchStatus(activeTimeBranch, 1);
    if (hasGunX) chain.SetBranchStatus("gun_x_mm", 1);

    chain.SetBranchAddress("event_id", &eventId);
    chain.SetBranchAddress("face_type", &faceType);
    chain.SetBranchAddress(activeTimeBranch, &timeNs);
    if (hasGunX) chain.SetBranchAddress("gun_x_mm", &gunXmm);

    std::vector<EventResult> results;
    const std::vector<double> kernel = BuildKernel(cfg.dtPs * 1e-3, cfg.tauRNs, cfg.tauFNs, cfg.tailWindowFactor);
    std::mt19937 rng(cfg.seed);

    std::vector<HitRecord> currentHits;
    int currentEvent = std::numeric_limits<int>::min();
    int currentFace = std::numeric_limits<int>::min();
    int currentTree = -999999;

    const Long64_t nEntries = chain.GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        chain.GetEntry(i);
        const int treeNumber = chain.GetTreeNumber();
        if (!FaceMatches(cfg.face, faceType)) continue;

        if (!currentHits.empty() &&
            (eventId != currentEvent || faceType != currentFace || treeNumber != currentTree)) {
            results.push_back(ProcessEvent(currentHits, kernel, cfg, rng));
            currentHits.clear();
            if (cfg.maxEvents > 0 && static_cast<long long>(results.size()) >= cfg.maxEvents) break;
        }

        if (cfg.maxEvents > 0 && static_cast<long long>(results.size()) >= cfg.maxEvents) break;

        HitRecord hit;
        hit.eventId = eventId;
        hit.faceType = faceType;
        hit.t0Ns = timeNs;
        hit.gunXmm = hasGunX ? gunXmm : 0.0;
        hit.treeNumber = treeNumber;
        hit.sourceFile = chain.GetCurrentFile() ? fs::path(chain.GetCurrentFile()->GetName()).filename().string() : "unknown";

        currentEvent = eventId;
        currentFace = faceType;
        currentTree = treeNumber;
        currentHits.push_back(hit);
    }

    if (!currentHits.empty() &&
        (cfg.maxEvents < 0 || static_cast<long long>(results.size()) < cfg.maxEvents)) {
        results.push_back(ProcessEvent(currentHits, kernel, cfg, rng));
    }

    return results;
}

}  // namespace

int main(int argc, char** argv) {
    const auto cfgOpt = ParseArgs(argc, argv);
    if (!cfgOpt.has_value()) {
        return (argc > 1 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) ? 0 : 1;
    }
    const Config cfg = *cfgOpt;

    TChain chain(cfg.treeName.c_str());
    int addedFiles = 0;
    for (const auto& input : cfg.inputs) {
        const int before = chain.GetNtrees();
        chain.Add(input.c_str());
        if (chain.GetNtrees() > before) {
            addedFiles += chain.GetNtrees() - before;
        }
    }

    if (chain.GetEntries() == 0 || addedFiles == 0) {
        std::cerr << "Error: no se encontraron entradas en el TTree '" << cfg.treeName << "'.\n";
        return 2;
    }

    const auto results = AnalyzeChain(chain, cfg);
    if (results.empty()) {
        std::cerr << "Error: no se reconstruyeron eventos para la seleccion pedida.\n";
        return 3;
    }

    WriteCsv(cfg.csvOut, results);
    const FitResult fit = FitResolution(results, cfg);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Archivos procesados: " << addedFiles << "\n";
    std::cout << "Entradas leidas: " << chain.GetEntries() << "\n";
    std::cout << "Seleccion de cara: " << cfg.face << "\n";
    std::cout << "Eventos reconstruidos: " << results.size() << "\n";
    std::cout << "Pulso Broadcom: tau_r = " << cfg.tauRNs << " ns, tau_f = " << cfg.tauFNs << " ns\n";
    std::cout << "dCFD: " << cfg.fraction * 100.0 << "%\n";
    std::cout << "Paso temporal: " << cfg.dtPs << " ps\n";
    std::cout << "Sigma transit-time: " << cfg.transitSigmaPs << " ps\n";
    std::cout << "Jitter electronica: " << cfg.electronicsSigmaPs << " ps\n";
    std::cout << "Tiempo medio reconstruido: " << fit.meanNs << " ns\n";
    std::cout << "Resolucion temporal sigma_t: " << fit.sigmaNs * 1e3 << " ps\n";
    std::cout << "Error estadistico del sigma: " << fit.sigmaErrNs * 1e3 << " ps\n";
    std::cout << "CSV: " << cfg.csvOut << "\n";
    if (cfg.savePlot) {
        std::cout << "Histograma: " << cfg.histOut << "\n";
    }

    return 0;
}
