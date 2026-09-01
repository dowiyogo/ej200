// SiPMRankingScan_v2.C
// Análisis de rendimiento de SiPMs individuales en scan longitudinal.
//
// Modos de operación (auto-detectados según los datos):
//   - SINGLE-EVENT: un muón por posición → rendimiento fotónico + t_first
//   - MULTI-EVENT : N muones por posición → distribución de time-marks → σ_t
//
// Uso:
//   root -l -b -q 'analysis/SiPMRankingScan_v2.C()'
//   root -l -b -q 'analysis/SiPMRankingScan_v2.C(kMeanFirstN, 3, 30)'
//
// Argumentos:
//   mode      : kFirstPhoton | kNthPhoton | kMeanFirstN
//   nParam    : N para kNthPhoton / kMeanFirstN (ignorado en kFirstPhoton)
//   minEvents : mínimo de eventos por (sensor,x) para sigma válida [multi-event]

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "TCanvas.h"
#include "TChain.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TString.h"
#include "TSystem.h"

// ── Enums & structs ──────────────────────────────────────────────────────────

enum TimeMarkMode {
    kFirstPhoton = 0,
    kNthPhoton   = 1,
    kMeanFirstN  = 2
};

struct SensorPos {
    int global_id;
    int face_type;
    int local_id;
};

struct SensorXData {
    // Siempre disponible (1 o más eventos)
    int    n_events   = 0;
    int    n_photons  = 0;        // total fotones detectados
    double t_first_ns = std::numeric_limits<double>::quiet_NaN();  // primer fotón
    double t_mean_ns  = std::numeric_limits<double>::quiet_NaN();  // media de todos los fotones
    double t_rms_ns   = std::numeric_limits<double>::quiet_NaN();  // RMS de todos los fotones

    // Solo con n_events > 1 (distribución de time-marks)
    double sigma_ns     = std::numeric_limits<double>::quiet_NaN();
    double sigma_err_ns = std::numeric_limits<double>::quiet_NaN();
    double mean_mark_ns = std::numeric_limits<double>::quiet_NaN();
};

// ── Utilidades ───────────────────────────────────────────────────────────────

namespace {

const char* ModeName(int mode) {
    switch (mode) {
        case kFirstPhoton: return "first-photon";
        case kNthPhoton:   return "nth-photon";
        case kMeanFirstN:  return "mean-first-n";
        default:           return "unknown";
    }
}

const char* FaceName(int face) {
    if (face == 0) return "end-left";
    if (face == 1) return "end-right";
    if (face == 2) return "top";
    return "unknown";
}

Double_t ComputeTimeMark(const std::vector<Double_t>& times, int mode, size_t nParam) {
    if (times.empty()) return -1.0;
    std::vector<Double_t> s = times;
    std::sort(s.begin(), s.end());
    if (mode == kFirstPhoton) return s.front();
    if (mode == kNthPhoton) {
        return (nParam == 0 || s.size() < nParam) ? -1.0 : s[nParam - 1];
    }
    if (mode == kMeanFirstN) {
        const size_t n = std::min(s.size(), nParam);
        if (n == 0) return -1.0;
        double sum = 0.0;
        for (size_t i = 0; i < n; ++i) sum += s[i];
        return sum / static_cast<double>(n);
    }
    return -1.0;
}

double VecMean(const std::vector<double>& v) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    double s = 0.0;
    for (double x : v) s += x;
    return s / static_cast<double>(v.size());
}

double VecRMS(const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    const double mu = VecMean(v);
    double s2 = 0.0;
    for (double x : v) { const double d = x - mu; s2 += d * d; }
    return std::sqrt(s2 / static_cast<double>(v.size() - 1));
}

// Ajuste gaussiano al núcleo de la distribución. Devuelve sigma [ns].
struct GausResult { double mean, sigma, sigma_err, rms; int n; bool fitted; };

GausResult FitGaussCore(const std::vector<double>& marks, const TString& name) {
    GausResult r{};
    r.n = static_cast<int>(marks.size());
    r.fitted = false;
    if (marks.empty()) return r;
    r.mean  = VecMean(marks);
    r.rms   = VecRMS(marks);
    r.sigma = r.rms;
    r.sigma_err = (r.n > 1) ? r.rms / std::sqrt(2.0 * (r.n - 1)) : 0.0;

    if (r.n < 20 || !(r.rms > 0.0)) return r;

    TH1D h(name, name, 80, r.mean - 4.0 * r.rms, r.mean + 4.0 * r.rms);
    h.SetDirectory(nullptr);
    for (double t : marks) h.Fill(t);

    TF1 f(name + "_g", "gaus", r.mean - 1.5 * r.rms, r.mean + 1.5 * r.rms);
    f.SetParameters(std::max(1.0, h.GetMaximum()), r.mean, std::max(0.02, r.rms));
    f.SetParLimits(2, 1e-4, 5.0 * r.rms);
    if (h.Fit(&f, "QNR") == 0 && std::isfinite(f.GetParameter(2)) && f.GetParameter(2) > 0.0) {
        r.mean      = f.GetParameter(1);
        r.sigma     = std::abs(f.GetParameter(2));
        r.sigma_err = f.GetParError(2);
        r.fitted    = true;
    }
    return r;
}

bool AddScanFiles(TChain& chain) {
    bool any = false;
    const std::vector<TString> dirs = {"", "build/", "../build/"};
    for (int i = 0; i <= 20; ++i) {
        bool got = false;
        for (const auto& d : dirs) {
            TString name = Form("%sphoton_hits_run%03d.root", d.Data(), i);
            if (!gSystem->AccessPathName(name)) {
                chain.Add(name);
                any = got = true;
                break;
            }
        }
        if (!got)
            std::cerr << "  [aviso] no se encontro photon_hits_run"
                      << Form("%03d", i) << ".root\n";
    }
    return any;
}

// Crea TH2D con ejes sensor × posición X
TH2D* MakeMap(const char* name, const char* title,
              const std::vector<int>& xPos, const std::vector<int>& sensors) {
    const int nx = static_cast<int>(xPos.size());
    const int ns = static_cast<int>(sensors.size());
    const double dx = (nx > 1) ? (xPos[1] - xPos[0]) / 2.0 : 32.5;
    auto* h = new TH2D(name, title,
                       nx, xPos.front() - dx, xPos.back() + dx,
                       ns, sensors.front() - 0.5, sensors.back() + 0.5);
    h->GetXaxis()->SetTitle("Posicion X del muon [mm]");
    h->GetYaxis()->SetTitle("global_id del SiPM");
    return h;
}

} // namespace

// ── Función principal ────────────────────────────────────────────────────────

void SiPMRankingScan_v2(Int_t mode = kFirstPhoton, Int_t nParam = 1, Int_t minEvents = 30) {

    // ── 1. Cargar datos ──────────────────────────────────────────────────────
    TChain chain("sipm_hits");
    std::cerr << "Buscando archivos de scan...\n";
    if (!AddScanFiles(chain)) {
        std::cerr << "Error: no se encontraron archivos del scan.\n";
        return;
    }
    const Long64_t totalEntries = chain.GetEntries();
    std::cout << "\nSiPMRankingScan_v2: modo=" << ModeName(mode)
              << "  nParam=" << nParam
              << "  minEvents=" << minEvents << "\n"
              << "  Entradas totales en el chain: " << totalEntries << "\n\n";

    Int_t event_id = -1, face_type = -1, global_id = -1, local_id = -1;
    Double_t time_ns = 0.0, gun_x_mm = 0.0;
    chain.SetBranchAddress("event_id",  &event_id);
    chain.SetBranchAddress("face_type", &face_type);
    chain.SetBranchAddress("global_id", &global_id);
    chain.SetBranchAddress("local_id",  &local_id);
    chain.SetBranchAddress("time_ns",   &time_ns);
    chain.SetBranchAddress("gun_x_mm",  &gun_x_mm);

    // ── 2. Acumular hits: sensor → {x,evento} → tiempos ────────────────────
    // También guardamos todos los tiempos por (sensor,x) para estadística bruta
    using EventKey = std::pair<int,int>;   // (x_rounded, event_id)
    std::map<int, std::map<EventKey, std::vector<double>>> hitsByEvent;   // per evento
    std::map<int, std::map<int, std::vector<double>>>      hitsByX;       // todos los fotones
    std::map<int, SensorPos> sensorInfo;
    std::set<int> xSet;
    std::set<int> eventSet; // event_ids globales para detectar modo

    for (Long64_t i = 0; i < totalEntries; ++i) {
        chain.GetEntry(i);
        const int x = static_cast<int>(std::lround(gun_x_mm));
        hitsByEvent[global_id][{x, event_id}].push_back(time_ns);
        hitsByX[global_id][x].push_back(time_ns);
        sensorInfo[global_id] = {global_id, face_type, local_id};
        xSet.insert(x);
        eventSet.insert(event_id);
    }

    const std::vector<int> xPositions(xSet.begin(), xSet.end());
    std::vector<int> sensors;
    for (const auto& kv : sensorInfo) sensors.push_back(kv.first);
    std::sort(sensors.begin(), sensors.end());

    if (xPositions.empty()) { std::cerr << "Error: sin posiciones X.\n"; return; }

    // ── 3. Diagnostico de modo ──────────────────────────────────────────────
    // Contar eventos únicos por posición (tomamos el máximo)
    std::map<int, std::set<int>> eventsPerX;
    for (Long64_t i = 0; i < totalEntries; ++i) {
        chain.GetEntry(i);
        eventsPerX[static_cast<int>(std::lround(gun_x_mm))].insert(event_id);
    }
    int maxEventsPerX = 0;
    for (const auto& kv : eventsPerX)
        maxEventsPerX = std::max(maxEventsPerX, (int)kv.second.size());

    const bool singleEventMode = (maxEventsPerX <= 1);

    std::cout << "── Diagnostico de datos ──────────────────────────────\n";
    std::cout << "  Posiciones X: " << xPositions.size() << "\n";
    std::cout << "  Sensores:     " << sensors.size() << "\n";
    std::cout << "  Eventos max por posicion: " << maxEventsPerX << "\n";
    if (singleEventMode) {
        std::cout << "\n  *** MODO SINGLE-EVENT ***\n"
                  << "  Solo hay 1 muon por posicion.\n"
                  << "  No es posible calcular sigma_t (requiere distribucion estadistica).\n"
                  << "  Se mostrara: rendimiento fotonico y tiempo de primer foton.\n"
                  << "  Para activar modo multi-event, regenerar con /run/beamOn >= 200\n\n";
    } else {
        std::cout << "  Modo multi-event activo.\n\n";
    }

    // ── 4. Calcular estadisticas por (sensor, x) ────────────────────────────
    // data[sensor][x] = SensorXData
    std::map<int, std::map<int, SensorXData>> data;

    for (int sid : sensors) {
        const auto& byX = hitsByX[sid];
        const auto& byEv = hitsByEvent[sid];

        for (int x : xPositions) {
            SensorXData d;

            // Fotones brutos
            auto itAll = byX.find(x);
            if (itAll != byX.end()) {
                const auto& allT = itAll->second;
                d.n_photons = static_cast<int>(allT.size());
                auto sorted = allT;
                std::sort(sorted.begin(), sorted.end());
                d.t_first_ns = sorted.front();
                d.t_mean_ns  = VecMean(allT);
                d.t_rms_ns   = VecRMS(allT);
            }

            // Eventos disponibles en esta (sensor,x)
            std::vector<double> marks;
            for (const auto& ev : byEv) {
                if (ev.first.first != x) continue;
                const double tm = ComputeTimeMark(ev.second, mode, static_cast<size_t>(nParam));
                if (tm >= 0.0) marks.push_back(tm);
            }
            d.n_events = static_cast<int>(marks.size());

            if (!singleEventMode && d.n_events >= minEvents) {
                const GausResult g = FitGaussCore(marks,
                    Form("h_%d_%d", sid, x));
                d.sigma_ns     = g.sigma;
                d.sigma_err_ns = g.sigma_err;
                d.mean_mark_ns = g.mean;
            }

            data[sid][x] = d;
        }
    }

    // ── 5. Histogramas 2D ────────────────────────────────────────────────────
    TH2D* hYield = MakeMap("hYield",
        "Fotones detectados por SiPM; X [mm]; global_id",
        xPositions, sensors);
    TH2D* hTfirst = MakeMap("hTfirst",
        "Tiempo de primer foton [ns]; X [mm]; global_id",
        xPositions, sensors);
    TH2D* hTmean = MakeMap("hTmean",
        "Tiempo medio de fotones [ns]; X [mm]; global_id",
        xPositions, sensors);

    TH2D* hSigma = nullptr;
    TH2D* hNev   = nullptr;
    if (!singleEventMode) {
        hSigma = MakeMap("hSigma",
            Form("#sigma time-mark [ps] (%s, n=%d); X [mm]; global_id",
                 ModeName(mode), nParam),
            xPositions, sensors);
        hNev = MakeMap("hNev",
            "Eventos por (sensor, X); X [mm]; global_id",
            xPositions, sensors);
    }

    for (int sid : sensors) {
        for (int x : xPositions) {
            const auto& d = data[sid][x];
            if (d.n_photons == 0) continue;
            hYield->SetBinContent(hYield->FindBin(x, sid), d.n_photons);
            if (std::isfinite(d.t_first_ns))
                hTfirst->SetBinContent(hTfirst->FindBin(x, sid), d.t_first_ns);
            if (std::isfinite(d.t_mean_ns))
                hTmean->SetBinContent(hTmean->FindBin(x, sid), d.t_mean_ns);
            if (!singleEventMode) {
                hNev->SetBinContent(hNev->FindBin(x, sid), d.n_events);
                if (std::isfinite(d.sigma_ns))
                    hSigma->SetBinContent(hSigma->FindBin(x, sid), d.sigma_ns * 1000.0);
            }
        }
    }

    // ── 6. Ranking por cara ─────────────────────────────────────────────────
    std::cout << "── Ranking por posicion X ────────────────────────────\n";

    // Variables para gráficas de tendencia
    struct FaceVec { std::vector<double> x, y; };
    std::map<int, FaceVec> faceYield, faceTfirst;
    std::map<int, FaceVec> faceSigma;

    for (int x : xPositions) {
        // Recopilar datos de todos los sensores en esta posición
        struct SEntry {
            int sid, face, local;
            double yield, tFirst, sigma;
        };
        std::vector<SEntry> entries;
        for (int sid : sensors) {
            const auto& d = data[sid][x];
            if (d.n_photons == 0) continue;
            SEntry e;
            e.sid   = sid;
            e.face  = sensorInfo[sid].face_type;
            e.local = sensorInfo[sid].local_id;
            e.yield  = d.n_photons;
            e.tFirst = std::isfinite(d.t_first_ns) ? d.t_first_ns : 1e9;
            e.sigma  = std::isfinite(d.sigma_ns) ? d.sigma_ns * 1000.0 : 1e9;
            entries.push_back(e);
        }
        if (entries.empty()) continue;

        // Mejor yield y t_first por cara
        for (int face = 0; face <= 2; ++face) {
            double bestYield = 0, bestTfirst = 1e9, bestSigma = 1e9;
            bool hasData = false;
            for (const auto& e : entries) {
                if (e.face != face) continue;
                if (e.yield  > bestYield)  bestYield  = e.yield;
                if (e.tFirst < bestTfirst) bestTfirst = e.tFirst;
                if (e.sigma  < bestSigma)  bestSigma  = e.sigma;
                hasData = true;
            }
            if (!hasData) continue;
            faceYield[face].x.push_back(x);
            faceYield[face].y.push_back(bestYield);
            faceTfirst[face].x.push_back(x);
            faceTfirst[face].y.push_back(bestTfirst);
            if (!singleEventMode && bestSigma < 1e8) {
                faceSigma[face].x.push_back(x);
                faceSigma[face].y.push_back(bestSigma);
            }
        }

        // Print ranking
        auto byYield = entries;
        std::sort(byYield.begin(), byYield.end(),
            [](const SEntry& a, const SEntry& b){ return a.yield > b.yield; });
        auto byTfirst = entries;
        std::sort(byTfirst.begin(), byTfirst.end(),
            [](const SEntry& a, const SEntry& b){ return a.tFirst < b.tFirst; });

        std::cout << "\n  x = " << std::setw(5) << x << " mm\n";
        std::cout << "    Top 3 por rendimiento fotonico:\n";
        for (size_t i = 0; i < std::min<size_t>(3, byYield.size()); ++i) {
            const auto& e = byYield[i];
            std::cout << "      " << i+1 << ". gid=" << e.sid
                      << " (" << FaceName(e.face) << " local=" << e.local << ")"
                      << "  fotones=" << (int)e.yield
                      << "  t_first=" << std::fixed << std::setprecision(3) << e.tFirst << " ns\n";
        }
        std::cout << "    Top 3 por t_first (mas rapido primero):\n";
        for (size_t i = 0; i < std::min<size_t>(3, byTfirst.size()); ++i) {
            const auto& e = byTfirst[i];
            std::cout << "      " << i+1 << ". gid=" << e.sid
                      << " (" << FaceName(e.face) << " local=" << e.local << ")"
                      << "  t_first=" << std::fixed << std::setprecision(3) << e.tFirst << " ns"
                      << "  fotones=" << (int)e.yield << "\n";
        }
        if (!singleEventMode) {
            auto bySigma = entries;
            std::sort(bySigma.begin(), bySigma.end(),
                [](const SEntry& a, const SEntry& b){ return a.sigma < b.sigma; });
            std::cout << "    Top 3 por sigma_t (mejor resolucion):\n";
            for (size_t i = 0; i < std::min<size_t>(3, bySigma.size()); ++i) {
                const auto& e = bySigma[i];
                if (e.sigma >= 1e8) break;
                std::cout << "      " << i+1 << ". gid=" << e.sid
                          << " (" << FaceName(e.face) << " local=" << e.local << ")"
                          << "  sigma=" << std::setprecision(0) << e.sigma << " ps"
                          << "  fotones=" << (int)e.yield << "\n";
            }
        }
    }

    // ── 7. Gráficas ─────────────────────────────────────────────────────────
    const int kColors[3] = {kBlue+1, kRed+1, kGreen+2};
    const int kMarkers[3] = {20, 21, 22};
    const char* kFaceLabels[3] = {"end-left", "end-right", "top"};

    // Canvas 1: mapas 2D
    {
        TCanvas* c = new TCanvas("cMaps", "Mapas 2D por SiPM", 1800, singleEventMode ? 400 : 600);
        const int ncols = singleEventMode ? 3 : 4;
        c->Divide(ncols, 1);

        c->cd(1); hYield->Draw("COLZ");
        gPad->SetRightMargin(0.14);
        gPad->SetBottomMargin(0.13);

        c->cd(2); hTfirst->Draw("COLZ");
        gPad->SetRightMargin(0.14);

        c->cd(3); hTmean->Draw("COLZ");
        gPad->SetRightMargin(0.14);

        if (!singleEventMode) {
            c->cd(4); hSigma->Draw("COLZ");
            gPad->SetRightMargin(0.14);
        }
        c->Update();
    }

    // Canvas 2: mejor rendimiento por cara vs x
    {
        TCanvas* c = new TCanvas("cYieldVsX", "Rendimiento fotonico por cara", 950, 520);
        TGraph* first = nullptr;
        TLegend* leg = new TLegend(0.62, 0.68, 0.90, 0.88);
        for (int face = 0; face <= 2; ++face) {
            auto& fv = faceYield[face];
            if (fv.x.empty()) continue;
            auto* g = new TGraph(static_cast<int>(fv.x.size()), fv.x.data(), fv.y.data());
            g->SetTitle("Mejor rendimiento fotonico por cara; X [mm]; n_{fotos} del mejor SiPM");
            g->SetLineColor(kColors[face]);
            g->SetMarkerColor(kColors[face]);
            g->SetMarkerStyle(kMarkers[face]);
            g->SetLineWidth(2);
            g->SetMinimum(0.0);
            if (!first) { g->Draw("ALP"); first = g; }
            else          g->Draw("LP SAME");
            leg->AddEntry(g, kFaceLabels[face], "lp");
        }
        if (first) leg->Draw();
        c->Update();
    }

    // Canvas 3: t_first del mejor SiPM por cara vs x
    {
        TCanvas* c = new TCanvas("cTfirstVsX", "Tiempo primer foton por cara", 950, 520);
        TGraph* first = nullptr;
        TLegend* leg = new TLegend(0.62, 0.68, 0.90, 0.88);
        for (int face = 0; face <= 2; ++face) {
            auto& fv = faceTfirst[face];
            if (fv.x.empty()) continue;
            auto* g = new TGraph(static_cast<int>(fv.x.size()), fv.x.data(), fv.y.data());
            g->SetTitle("Tiempo de primer foton (mejor SiPM por cara); X [mm]; t_{first} [ns]");
            g->SetLineColor(kColors[face]);
            g->SetMarkerColor(kColors[face]);
            g->SetMarkerStyle(kMarkers[face]);
            g->SetLineWidth(2);
            if (!first) { g->Draw("ALP"); first = g; }
            else          g->Draw("LP SAME");
            leg->AddEntry(g, kFaceLabels[face], "lp");
        }
        if (first) leg->Draw();
        c->Update();
    }

    // Canvas 4: sigma por cara (solo multi-event)
    if (!singleEventMode && !faceSigma.empty()) {
        TCanvas* c = new TCanvas("cSigmaVsX", "Resolucion temporal por cara", 950, 520);
        TGraph* first = nullptr;
        TLegend* leg = new TLegend(0.62, 0.68, 0.90, 0.88);
        for (int face = 0; face <= 2; ++face) {
            auto& fv = faceSigma[face];
            if (fv.x.empty()) continue;
            auto* g = new TGraph(static_cast<int>(fv.x.size()), fv.x.data(), fv.y.data());
            g->SetTitle(Form("Mejor #sigma_{core} por cara (%s, n=%d); X [mm]; #sigma [ps]",
                             ModeName(mode), nParam));
            g->SetLineColor(kColors[face]);
            g->SetMarkerColor(kColors[face]);
            g->SetMarkerStyle(kMarkers[face]);
            g->SetLineWidth(2);
            g->SetMinimum(0.0);
            if (!first) { g->Draw("ALP"); first = g; }
            else          g->Draw("LP SAME");
            leg->AddEntry(g, kFaceLabels[face], "lp");
        }
        if (first) leg->Draw();
        c->Update();
    }

    // ── 8. Resumen global de sensores ───────────────────────────────────────
    std::cout << "\n── Resumen global de sensores ────────────────────────\n";
    struct GlobalSensor {
        int sid, face, local;
        double totalPhotons;
        double avgTfirst;
        double avgSigmaPs;
        int    nxWithData;
    };
    std::vector<GlobalSensor> globalRank;
    for (int sid : sensors) {
        GlobalSensor gs;
        gs.sid   = sid;
        gs.face  = sensorInfo[sid].face_type;
        gs.local = sensorInfo[sid].local_id;
        gs.totalPhotons = 0;
        double sumTfirst = 0; int nTfirst = 0;
        double sumSigma  = 0; int nSigma  = 0;
        gs.nxWithData = 0;
        for (int x : xPositions) {
            const auto& d = data[sid][x];
            if (d.n_photons == 0) continue;
            gs.nxWithData++;
            gs.totalPhotons += d.n_photons;
            if (std::isfinite(d.t_first_ns)) { sumTfirst += d.t_first_ns; nTfirst++; }
            if (std::isfinite(d.sigma_ns))   { sumSigma  += d.sigma_ns * 1000.0; nSigma++; }
        }
        gs.avgTfirst   = nTfirst > 0 ? sumTfirst / nTfirst : 1e9;
        gs.avgSigmaPs  = nSigma  > 0 ? sumSigma  / nSigma  : 1e9;
        if (gs.nxWithData > 0) globalRank.push_back(gs);
    }

    // Ordenar por total de fotones
    auto byPhotons = globalRank;
    std::sort(byPhotons.begin(), byPhotons.end(),
        [](const GlobalSensor& a, const GlobalSensor& b){
            return a.totalPhotons > b.totalPhotons;
        });
    std::cout << "  Top 5 por total de fotones detectados (todas las posiciones):\n";
    for (size_t i = 0; i < std::min<size_t>(5, byPhotons.size()); ++i) {
        const auto& g = byPhotons[i];
        std::cout << "    " << i+1 << ". gid=" << g.sid
                  << " (" << FaceName(g.face) << " local=" << g.local << ")"
                  << "  total=" << (int)g.totalPhotons
                  << "  n_posiciones=" << g.nxWithData << "\n";
    }

    if (!singleEventMode) {
        auto bySigma = globalRank;
        std::sort(bySigma.begin(), bySigma.end(),
            [](const GlobalSensor& a, const GlobalSensor& b){
                return a.avgSigmaPs < b.avgSigmaPs;
            });
        std::cout << "  Top 5 por sigma promedio [ps]:\n";
        for (size_t i = 0; i < std::min<size_t>(5, bySigma.size()); ++i) {
            const auto& g = bySigma[i];
            if (g.avgSigmaPs >= 1e8) break;
            std::cout << "    " << i+1 << ". gid=" << g.sid
                      << " (" << FaceName(g.face) << " local=" << g.local << ")"
                      << "  avg_sigma=" << std::setprecision(0) << g.avgSigmaPs << " ps\n";
        }
    }

    if (singleEventMode) {
        std::cout << "\n╔════════════════════════════════════════════════════╗\n"
                  << "║ ADVERTENCIA: datos con 1 evento por posicion.      ║\n"
                  << "║ Para resolucion temporal σ_t regenerar con:        ║\n"
                  << "║   /run/beamOn 500  (o mas)  en scan_step.mac       ║\n"
                  << "║ El scan actual ya tiene beamOn 2000 configurado.   ║\n"
                  << "╚════════════════════════════════════════════════════╝\n\n";
    }
}
