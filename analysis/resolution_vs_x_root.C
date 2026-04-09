#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TLine.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TMultiGraph.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr double BAR_HALF_X = 700.0;

struct EventAccum {
    int x_mm = 0;
    int event_id = -1;

    int nLeft = 0;
    int nRight = 0;
    int nTop = 0;

    bool hasLeft = false;
    bool hasRight = false;
    bool hasTop = false;
    bool hasEnd = false;

    double minLeft = std::numeric_limits<double>::infinity();
    double minRight = std::numeric_limits<double>::infinity();
    double minTop = std::numeric_limits<double>::infinity();
    double minEnd = std::numeric_limits<double>::infinity();
};

struct FitResult {
    double mu_ns = std::numeric_limits<double>::quiet_NaN();
    double sigma_ns = std::numeric_limits<double>::quiet_NaN();
    double sigma_err_ns = std::numeric_limits<double>::quiet_NaN();
    int nEvents = 0;
};

struct ProfilePoint {
    double x_mm = 0.0;
    double y = 0.0;
    double ey = 0.0;
    int n = 0;
};

struct SummaryByX {
    std::vector<double> fptLeft;
    std::vector<double> fptRight;
    std::vector<double> fptEnd;
    std::vector<double> fptTop;
    std::vector<double> asym;
    std::vector<double> nLeft;
    std::vector<double> nRight;
    std::vector<double> nTop;
};

uint64_t make_key(int x_mm, int event_id) {
    const uint64_t ux = static_cast<uint32_t>(x_mm + 100000);
    const uint64_t ue = static_cast<uint32_t>(event_id);
    return (ux << 32) | ue;
}

double vec_mean(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    double s = 0.0;
    for (double x : v) s += x;
    return s / static_cast<double>(v.size());
}

double vec_std(const std::vector<double>& v) {
    if (v.size() < 2) return 0.0;
    const double m = vec_mean(v);
    double s2 = 0.0;
    for (double x : v) {
        const double d = x - m;
        s2 += d * d;
    }
    return std::sqrt(s2 / static_cast<double>(v.size() - 1));
}

double percentile(std::vector<double> v, double q) {
    if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
    std::sort(v.begin(), v.end());
    const double pos = q * (v.size() - 1);
    const size_t lo = static_cast<size_t>(std::floor(pos));
    const size_t hi = static_cast<size_t>(std::ceil(pos));
    const double frac = pos - lo;
    if (lo == hi) return v[lo];
    return v[lo] * (1.0 - frac) + v[hi] * frac;
}

FitResult fit_fpt_distribution(const std::vector<double>& times_ns, const std::string& hname) {
    FitResult out;
    out.nEvents = static_cast<int>(times_ns.size());
    if (times_ns.size() < 10) return out;

    const double q1 = percentile(times_ns, 0.01);
    const double q99 = percentile(times_ns, 0.99);
    double lo = q1;
    double hi = q99;
    const double margin = std::max(0.5 * (hi - lo), 0.5);
    lo -= margin;
    hi += margin;
    if (!(hi > lo)) {
        lo = *std::min_element(times_ns.begin(), times_ns.end()) - 0.5;
        hi = *std::max_element(times_ns.begin(), times_ns.end()) + 0.5;
    }

    TH1D h(hname.c_str(), "", 40, lo, hi);
    for (double t : times_ns) h.Fill(t);

    const double mu0 = vec_mean(times_ns);
    double sigma0 = vec_std(times_ns);
    if (sigma0 < 1e-6) sigma0 = 0.1;
    const double A0 = std::max(1.0, h.GetMaximum());

    TF1 gaus((hname + "_gaus").c_str(), "gaus", lo, hi);
    gaus.SetParameters(A0, mu0, sigma0);
    gaus.SetParLimits(2, 1e-4, hi - lo);

    const int fit_status = h.Fit(&gaus, "QNR");
    if (fit_status == 0) {
        out.mu_ns = gaus.GetParameter(1);
        out.sigma_ns = std::abs(gaus.GetParameter(2));
        out.sigma_err_ns = gaus.GetParError(2);
    } else {
        out.mu_ns = mu0;
        out.sigma_ns = sigma0;
        out.sigma_err_ns = sigma0 / std::sqrt(static_cast<double>(times_ns.size()));
    }
    return out;
}

std::vector<ProfilePoint> build_profile_mean_sem(const std::map<int, std::vector<double>>& byX) {
    std::vector<ProfilePoint> out;
    for (const auto& kv : byX) {
        const auto& v = kv.second;
        if (v.empty()) continue;
        const double m = vec_mean(v);
        const double s = vec_std(v);
        const double sem = (v.size() > 1) ? s / std::sqrt(static_cast<double>(v.size())) : 0.0;
        out.push_back(ProfilePoint{static_cast<double>(kv.first), m, sem, static_cast<int>(v.size())});
    }
    return out;
}

TGraphErrors* make_graph(const std::vector<ProfilePoint>& pts, const std::string& name) {
    auto* g = new TGraphErrors(static_cast<int>(pts.size()));
    g->SetName(name.c_str());
    for (int i = 0; i < static_cast<int>(pts.size()); ++i) {
        g->SetPoint(i, pts[i].x_mm, pts[i].y);
        g->SetPointError(i, 0.0, pts[i].ey);
    }
    return g;
}

void style_graph(TGraphErrors* g, int color, int marker, int lineStyle = 1, double alpha = 1.0) {
    g->SetLineColorAlpha(color, alpha);
    g->SetMarkerColorAlpha(color, alpha);
    g->SetLineWidth(2);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.0);
    g->SetLineStyle(lineStyle);
}

void save_resolution_plot(const std::map<int, SummaryByX>& byX, const std::string& outPdf) {
    std::vector<ProfilePoint> endLeft, endRight, endBoth, top;
    for (const auto& kv : byX) {
        const int x = kv.first;
        auto add_fit = [&](const std::vector<double>& times, std::vector<ProfilePoint>& dest, const std::string& base) {
            FitResult fr = fit_fpt_distribution(times, base + std::to_string(x));
            if (std::isfinite(fr.sigma_ns) && fr.sigma_ns > 0) {
                dest.push_back(ProfilePoint{static_cast<double>(x), fr.sigma_ns * 1e3, fr.sigma_err_ns * 1e3, fr.nEvents});
            }
        };
        add_fit(kv.second.fptLeft, endLeft, "h_fpt_left_");
        add_fit(kv.second.fptRight, endRight, "h_fpt_right_");
        add_fit(kv.second.fptEnd, endBoth, "h_fpt_end_");
        add_fit(kv.second.fptTop, top, "h_fpt_top_");
    }

    auto* gLeft = make_graph(endLeft, "gLeft");
    auto* gRight = make_graph(endRight, "gRight");
    auto* gEnd = make_graph(endBoth, "gEnd");
    auto* gTop = make_graph(top, "gTop");

    style_graph(gLeft, kBlue + 1, 22, 2, 0.6);
    style_graph(gRight, kRed + 1, 23, 2, 0.6);
    style_graph(gEnd, kViolet + 1, 20, 1, 1.0);
    style_graph(gTop, kGreen + 2, 21, 1, 1.0);

    TCanvas c("c_resolution", "resolution_vs_x", 1100, 550);
    c.SetGrid();

    TMultiGraph mg;
    mg.Add(gEnd, "LP");
    mg.Add(gTop, "LP");
    mg.Add(gLeft, "LP");
    mg.Add(gRight, "LP");
    mg.Draw("A");
    mg.GetXaxis()->SetTitle("Posicion longitudinal del muon x [mm]");
    mg.GetYaxis()->SetTitle("Resolucion temporal #sigma_{t} [ps]");
    mg.GetXaxis()->SetLimits(-BAR_HALF_X - 30.0, BAR_HALF_X + 30.0);
    mg.SetMinimum(0.0);

    TLine target(-BAR_HALF_X - 30.0, 100.0, BAR_HALF_X + 30.0, 100.0);
    target.SetLineStyle(3);
    target.SetLineWidth(2);
    target.Draw();

    TLegend leg(0.60, 0.62, 0.90, 0.88);
    leg.AddEntry(gEnd, "End SiPMs (L+R combinados)", "lp");
    leg.AddEntry(gTop, "Top SiPMs", "lp");
    leg.AddEntry(gLeft, "End izquierdo (-X)", "lp");
    leg.AddEntry(gRight, "End derecho (+X)", "lp");
    leg.AddEntry(&target, "Objetivo TD SHiP (100 ps)", "l");
    leg.Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.12, 0.94, "Resolucion temporal vs posicion longitudinal - EJ-200");

    c.SaveAs(outPdf.c_str());
}

void save_asymmetry_plot(const std::map<int, SummaryByX>& byX, const std::string& outPdf) {
    std::map<int, std::vector<double>> asymByX;
    for (const auto& kv : byX) asymByX[kv.first] = kv.second.asym;
    auto pts = build_profile_mean_sem(asymByX);
    auto* g = make_graph(pts, "gAsym");
    style_graph(g, kViolet + 1, 20);

    TCanvas c("c_asym", "asymmetry_vs_x", 1000, 450);
    c.SetGrid();
    g->Draw("AP");
    g->GetXaxis()->SetTitle("Posicion del muon x [mm]");
    g->GetYaxis()->SetTitle("Asimetria (N_{L}-N_{R})/(N_{L}+N_{R})");
    g->GetXaxis()->SetLimits(-BAR_HALF_X - 30.0, BAR_HALF_X + 30.0);
    g->SetMinimum(-1.05);
    g->SetMaximum(1.05);

    TLine zero(-BAR_HALF_X - 30.0, 0.0, BAR_HALF_X + 30.0, 0.0);
    zero.SetLineStyle(3);
    zero.Draw();

    if (pts.size() >= 2) {
        double sx=0, sy=0, sxx=0, sxy=0;
        const double n = static_cast<double>(pts.size());
        for (const auto& p : pts) {
            sx += p.x_mm;
            sy += p.y;
            sxx += p.x_mm * p.x_mm;
            sxy += p.x_mm * p.y;
        }
        const double denom = n * sxx - sx * sx;
        if (std::abs(denom) > 1e-12) {
            const double a = (n * sxy - sx * sy) / denom;
            const double b = (sy - a * sx) / n;
            TF1 fline("fline", "[0]*x+[1]", -BAR_HALF_X, BAR_HALF_X);
            fline.SetParameters(a, b);
            fline.SetLineStyle(2);
            fline.SetLineColor(kBlack);
            fline.Draw("same");
            TLegend leg(0.58, 0.75, 0.90, 0.90);
            leg.AddEntry(g, "<A> #pm SEM", "lp");
            leg.AddEntry(&fline, Form("Ajuste lineal: slope = %.2f m^{-1}", a * 1000.0), "l");
            leg.Draw();
        }
    }

    c.SaveAs(outPdf.c_str());
}

void save_nphotons_plot(const std::map<int, SummaryByX>& byX, const std::string& outPdf) {
    std::map<int, std::vector<double>> nLeftByX, nRightByX, nTopByX;
    for (const auto& kv : byX) {
        nLeftByX[kv.first] = kv.second.nLeft;
        nRightByX[kv.first] = kv.second.nRight;
        nTopByX[kv.first] = kv.second.nTop;
    }

    auto* gLeft = make_graph(build_profile_mean_sem(nLeftByX), "gNLeft");
    auto* gRight = make_graph(build_profile_mean_sem(nRightByX), "gNRight");
    auto* gTop = make_graph(build_profile_mean_sem(nTopByX), "gNTop");
    style_graph(gLeft, kBlue + 1, 20);
    style_graph(gRight, kRed + 1, 20);
    style_graph(gTop, kGreen + 2, 21);

    TCanvas c("c_nphot", "n_photons_vs_x", 1100, 450);
    c.SetGrid();
    TMultiGraph mg;
    mg.Add(gLeft, "LP");
    mg.Add(gRight, "LP");
    mg.Add(gTop, "LP");
    mg.Draw("A");
    mg.GetXaxis()->SetTitle("Posicion del muon x [mm]");
    mg.GetYaxis()->SetTitle("Fotones detectados por evento");
    mg.GetXaxis()->SetLimits(-BAR_HALF_X - 30.0, BAR_HALF_X + 30.0);

    TLegend leg(0.70, 0.73, 0.90, 0.88);
    leg.AddEntry(gLeft, "End izq. (-X)", "lp");
    leg.AddEntry(gRight, "End der. (+X)", "lp");
    leg.AddEntry(gTop, "Top (+Y)", "lp");
    leg.Draw();

    c.SaveAs(outPdf.c_str());
}

std::vector<int> choose_example_positions(const std::map<int, SummaryByX>& byX) {
    std::vector<int> xvals;
    for (const auto& kv : byX) xvals.push_back(kv.first);
    if (xvals.empty()) return {};
    std::vector<double> targets = {-0.9 * BAR_HALF_X, 0.0, 0.9 * BAR_HALF_X};
    std::vector<int> chosen;
    for (double t : targets) {
        int best = xvals.front();
        double bestd = std::abs(best - t);
        for (int x : xvals) {
            const double d = std::abs(x - t);
            if (d < bestd) {
                bestd = d;
                best = x;
            }
        }
        if (std::find(chosen.begin(), chosen.end(), best) == chosen.end()) chosen.push_back(best);
    }
    std::sort(chosen.begin(), chosen.end());
    return chosen;
}

void save_fpt_examples(const std::map<int, SummaryByX>& byX, bool useTop, const std::string& outPdf) {
    auto chosen = choose_example_positions(byX);
    if (chosen.empty()) return;

    const int n = static_cast<int>(chosen.size());
    TCanvas c(useTop ? "c_fpt_top" : "c_fpt_end", "fpt_examples", 450 * n, 450);
    c.Divide(n, 1);

    for (int i = 0; i < n; ++i) {
        const int x = chosen[i];
        const auto& vec = useTop ? byX.at(x).fptTop : byX.at(x).fptEnd;
        if (vec.empty()) continue;

        const double q1 = percentile(vec, 0.01);
        const double q99 = percentile(vec, 0.99);
        double lo = q1;
        double hi = q99;
        const double margin = std::max(0.5 * (hi - lo), 0.5);
        lo -= margin;
        hi += margin;
        if (!(hi > lo)) {
            lo = *std::min_element(vec.begin(), vec.end()) - 0.5;
            hi = *std::max_element(vec.begin(), vec.end()) + 0.5;
        }

        c.cd(i + 1);
        gPad->SetGrid();
        TH1D h(Form("h_%s_%d", useTop ? "top" : "end", x), "", 30, lo, hi);
        for (double t : vec) h.Fill(t);
        h.SetFillColorAlpha(kGray + 1, 0.6);
        h.SetLineColor(kBlack);
        h.GetXaxis()->SetTitle("FPT [ns]");
        h.GetYaxis()->SetTitle(i == 0 ? "Eventos / bin" : "");
        h.SetTitle(Form("%s; x = %d mm", useTop ? "Top SiPMs" : "End SiPMs", x));
        h.Draw("hist");

        FitResult fr = fit_fpt_distribution(vec, Form("hfit_%s_%d", useTop ? "top" : "end", x));
        if (std::isfinite(fr.sigma_ns)) {
            TF1 gaus(Form("g_%s_%d", useTop ? "top" : "end", x), "gaus", lo, hi);
            gaus.SetParameters(h.GetMaximum(), fr.mu_ns, fr.sigma_ns);
            gaus.SetLineColor(kRed + 1);
            gaus.SetLineWidth(2);
            gaus.Draw("same");

            TLatex t;
            t.SetNDC();
            t.SetTextSize(0.045);
            t.DrawLatex(0.55, 0.84, Form("#mu = %.2f ns", fr.mu_ns));
            t.DrawLatex(0.55, 0.76, Form("#sigma = %.0f #pm %.0f ps", fr.sigma_ns * 1e3, fr.sigma_err_ns * 1e3));
        }
    }

    c.SaveAs(outPdf.c_str());
}

void save_summary_csv(const std::map<int, SummaryByX>& byX, const std::string& outCsv) {
    std::ofstream out(outCsv);
    out << "x_mm,sigma_end_left_ps,sigma_end_right_ps,sigma_end_both_ps,sigma_top_ps,n_left,n_right,n_top\n";
    for (const auto& kv : byX) {
        const int x = kv.first;
        const auto frL = fit_fpt_distribution(kv.second.fptLeft, "tmpL" + std::to_string(x));
        const auto frR = fit_fpt_distribution(kv.second.fptRight, "tmpR" + std::to_string(x));
        const auto frE = fit_fpt_distribution(kv.second.fptEnd, "tmpE" + std::to_string(x));
        const auto frT = fit_fpt_distribution(kv.second.fptTop, "tmpT" + std::to_string(x));
        out << x << ","
            << frL.sigma_ns * 1e3 << ","
            << frR.sigma_ns * 1e3 << ","
            << frE.sigma_ns * 1e3 << ","
            << frT.sigma_ns * 1e3 << ","
            << kv.second.nLeft.size() << ","
            << kv.second.nRight.size() << ","
            << kv.second.nTop.size() << "\n";
    }
}

} // namespace

void resolution_vs_x_root(const char* rootFile = "photon_hits.root") {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TFile* f = TFile::Open(rootFile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: no pude abrir '" << rootFile << "'\n";
        return;
    }

    TTree* tree = dynamic_cast<TTree*>(f->Get("sipm_hits"));
    if (!tree) {
        std::cerr << "ERROR: el archivo no contiene el TTree 'sipm_hits'\n";
        f->Close();
        return;
    }

    TTreeReader reader(tree);
    TTreeReaderValue<int> event_id(reader, "event_id");
    TTreeReaderValue<int> face_type(reader, "face_type");
    TTreeReaderValue<double> time_ns(reader, "time_ns");
    TTreeReaderValue<double> gun_x_mm(reader, "gun_x_mm");

    std::unordered_map<uint64_t, EventAccum> events;
    events.reserve(static_cast<size_t>(tree->GetEntries() / 10 + 1000));

    long long nEntries = 0;
    while (reader.Next()) {
        ++nEntries;
        const int x = static_cast<int>(std::llround(*gun_x_mm));
        const int evt = *event_id;
        const int face = *face_type;
        const double t = *time_ns;

        const uint64_t key = make_key(x, evt);
        auto& acc = events[key];
        acc.x_mm = x;
        acc.event_id = evt;

        if (face == 0) {
            ++acc.nLeft;
            acc.hasLeft = true;
            acc.hasEnd = true;
            if (t < acc.minLeft) acc.minLeft = t;
            if (t < acc.minEnd) acc.minEnd = t;
        } else if (face == 1) {
            ++acc.nRight;
            acc.hasRight = true;
            acc.hasEnd = true;
            if (t < acc.minRight) acc.minRight = t;
            if (t < acc.minEnd) acc.minEnd = t;
        } else if (face == 2) {
            ++acc.nTop;
            acc.hasTop = true;
            if (t < acc.minTop) acc.minTop = t;
        }
    }

    std::cout << "Entradas leidas del TTree: " << nEntries << "\n";
    std::cout << "Eventos unicos (x,event_id): " << events.size() << "\n";

    std::map<int, SummaryByX> byX;
    for (const auto& kv : events) {
        const auto& e = kv.second;
        auto& sx = byX[e.x_mm];

        if (e.hasLeft) sx.fptLeft.push_back(e.minLeft);
        if (e.hasRight) sx.fptRight.push_back(e.minRight);
        if (e.hasEnd) sx.fptEnd.push_back(e.minEnd);
        if (e.hasTop) sx.fptTop.push_back(e.minTop);

        if (e.nLeft + e.nRight > 0) {
            const double asym = static_cast<double>(e.nLeft - e.nRight) / static_cast<double>(e.nLeft + e.nRight);
            sx.asym.push_back(asym);
        }

        sx.nLeft.push_back(static_cast<double>(e.nLeft));
        sx.nRight.push_back(static_cast<double>(e.nRight));
        sx.nTop.push_back(static_cast<double>(e.nTop));
    }

    std::cout << "Posiciones x encontradas: " << byX.size();
    if (!byX.empty()) {
        std::cout << " (" << byX.begin()->first << " a " << byX.rbegin()->first << " mm)";
    }
    std::cout << "\n";

    save_resolution_plot(byX, "resolution_vs_x.pdf");
    save_fpt_examples(byX, false, "fpt_dist_end.pdf");
    save_fpt_examples(byX, true, "fpt_dist_top.pdf");
    save_asymmetry_plot(byX, "asymmetry_vs_x.pdf");
    save_nphotons_plot(byX, "n_photons_vs_x.pdf");
    save_summary_csv(byX, "resolution_vs_x_summary.csv");

    std::cout << "Listo. Archivos generados:\n"
              << "  resolution_vs_x.pdf\n"
              << "  fpt_dist_end.pdf\n"
              << "  fpt_dist_top.pdf\n"
              << "  asymmetry_vs_x.pdf\n"
              << "  n_photons_vs_x.pdf\n"
              << "  resolution_vs_x_summary.csv\n";

    f->Close();
}
