// fpt_vs_n_profile.C
//
// Para cada SiPM construye un TProfile con <t_n> vs n promediando sobre
// todos los eventos del archivo.  Un canvas por cara (end-left, end-right,
// top).  En el titulo se muestra gid + posicion fisica relativa a la cara,
// y detras del grafico un watermark con la posicion del haz de muones.
//
// Uso:
//   root -l 'fpt_vs_n_profile.C("photon_hits_run000.root")'
//   root -l 'fpt_vs_n_profile.C("photon_hits_run000.root", 100)'
//   root -l 'fpt_vs_n_profile.C("photon_hits_run000.root", 150, 0)'  // filtra gun_x=0
//
// TProfile con opcion "s": la barra de error es la RMS de t_n sobre
// eventos (dispersion fisica, no SEM).

void fpt_vs_n_profile(const char* fname = "photon_hits_run000.root",
                      int nMax = 150,
                      double gunXSelect = -9999.0) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleFontSize(0.09);   // tamano del titulo del pad (clave
                                      // para que se lea en grids densos)

    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) { printf("ERROR: no puedo abrir %s\n", fname); return; }
    TTree* t = (TTree*)f->Get("sipm_hits");
    if (!t) { printf("ERROR: sin TTree sipm_hits\n"); return; }

    Int_t    ev, gid, face, lid;
    Double_t tns, xmm, ymm, gunX;
    t->SetBranchAddress("event_id",  &ev);
    t->SetBranchAddress("global_id", &gid);
    t->SetBranchAddress("local_id",  &lid);
    t->SetBranchAddress("face_type", &face);
    t->SetBranchAddress("time_ns",   &tns);
    t->SetBranchAddress("x_mm",      &xmm);
    t->SetBranchAddress("y_mm",      &ymm);
    t->SetBranchAddress("gun_x_mm",  &gunX);

    std::map<std::pair<int,int>, std::vector<double>> hits;
    std::map<int,int>    faceByGid, localByGid;
    std::map<int,double> sumX, sumY;
    std::map<int,long>   nHitsPerGid;
    std::set<int>        gunXSet;

    const Long64_t nE = t->GetEntries();
    const bool useFilter = (gunXSelect > -9998.0);
    Long64_t nKept = 0;
    for (Long64_t i = 0; i < nE; ++i) {
        t->GetEntry(i);
        if (useFilter && std::abs(gunX - gunXSelect) > 0.5) continue;
        hits[{ev, gid}].push_back(tns);
        faceByGid[gid]  = face;
        localByGid[gid] = lid;
        sumX[gid] += xmm;
        sumY[gid] += ymm;
        nHitsPerGid[gid]++;
        gunXSet.insert((int)std::round(gunX));
        ++nKept;
    }
    printf("Leidas %lld entradas, %lld aceptadas\n", nE, nKept);
    printf("Posiciones gun_x: ");
    int cnt = 0;
    for (int x : gunXSet) {
        printf("%d ", x);
        if (++cnt > 10) { printf("..."); break; }
    }
    printf("(%zu total)\n", gunXSet.size());
    if (hits.empty()) { printf("Sin datos. Abortando.\n"); return; }

    // Posicion media por SiPM (los pre-step points caen sobre su cara)
    std::map<int,double> meanX, meanY;
    for (auto& kv : nHitsPerGid) {
        int g = kv.first;
        double n = (double)kv.second;
        meanX[g] = sumX[g] / n;
        meanY[g] = sumY[g] / n;
    }

    std::vector<int> gidLeft, gidRight, gidTop;
    for (auto& kv : faceByGid) {
        if      (kv.second == 0) gidLeft.push_back(kv.first);
        else if (kv.second == 1) gidRight.push_back(kv.first);
        else if (kv.second == 2) gidTop.push_back(kv.first);
    }
    std::sort(gidLeft.begin(),  gidLeft.end());
    std::sort(gidRight.begin(), gidRight.end());
    std::sort(gidTop.begin(),   gidTop.end());

    std::map<int,TProfile*> prof;
    for (auto& kv : faceByGid) {
        const int g   = kv.first;
        const int fc  = kv.second;
        const int loc = localByGid[g];
        double      relPos   = (fc == 2) ? meanX[g] : meanY[g];
        const char* axisName = (fc == 2) ? "x"      : "y";
        TString name  = TString::Format("prof_gid%d", g);
        TString title = TString::Format(
            "gid=%d (L=%d)  %s=%+.0f mm;n-esimo foton;#LTt_{n}#GT [ns]",
            g, loc, axisName, relPos);
        prof[g] = new TProfile(name, title, nMax, 0.5, nMax + 0.5, "s");
    }

    for (auto& kv : hits) {
        const int g = kv.first.second;
        auto& times = kv.second;
        if (!prof.count(g)) continue;
        std::sort(times.begin(), times.end());
        const int nUse = std::min((int)times.size(), nMax);
        for (int n = 0; n < nUse; ++n) prof[g]->Fill(n + 1, times[n]);
    }

    TString wmText;
    if (useFilter) {
        wmText = TString::Format("#mu @ x = %+.0f mm", gunXSelect);
    } else if (gunXSet.size() == 1) {
        wmText = TString::Format("#mu @ x = %+d mm", *gunXSet.begin());
    } else {
        wmText = TString::Format("#mu scan (%zu pos.)", gunXSet.size());
    }

    auto draw = [&](const std::vector<int>& gids, const char* cname,
                    const char* ctitle, int color) {
        if (gids.empty()) { printf("  %s: sin SiPMs\n", cname); return; }
        const int n = (int)gids.size();
        const int ncols = (int)std::ceil(std::sqrt((double)n));
        const int nrows = (int)std::ceil((double)n / ncols);
        auto* c = new TCanvas(cname, ctitle, 290 * ncols, 240 * nrows);
        c->Divide(ncols, nrows, 0.002, 0.002);
        for (int i = 0; i < n; ++i) {
            c->cd(i + 1);
            gPad->SetGrid();
            gPad->SetLeftMargin(0.18);
            gPad->SetBottomMargin(0.17);
            gPad->SetTopMargin(0.14);  // espacio para el titulo

            auto* p = prof[gids[i]];
            p->SetMarkerStyle(20);
            p->SetMarkerSize(0.5);
            p->SetMarkerColor(color);
            p->SetLineColor(color);
            p->GetXaxis()->SetTitleSize(0.055);
            p->GetYaxis()->SetTitleSize(0.055);
            p->GetXaxis()->SetLabelSize(0.05);
            p->GetYaxis()->SetLabelSize(0.05);
            p->GetXaxis()->SetTitleOffset(1.1);
            p->GetYaxis()->SetTitleOffset(1.5);
            p->Draw("E1");             // <- dibuja ejes + titulo + datos

            // Watermark encima (visualmente "detras" por alpha bajo)
            TLatex wm;
            wm.SetNDC();
            wm.SetTextAlign(22);
            wm.SetTextSize(0.08);
            wm.SetTextAngle(20);
            wm.SetTextColorAlpha(kGray, 0.12);
            wm.DrawLatex(0.55, 0.50, wmText);
        }
        c->Update();
        printf("  %s: %d SiPMs dibujados\n", cname, n);
    };

    printf("Dibujando canvas separados por cara...\n");
    draw(gidLeft,  "cLeft",  "End-left SiPMs: <t_n> vs n",  kBlue+1);
    draw(gidRight, "cRight", "End-right SiPMs: <t_n> vs n", kRed+1);
    draw(gidTop,   "cTop",   "Top SiPMs: <t_n> vs n",       kGreen+2);

    printf("Listo: %zu end-left, %zu end-right, %zu top\n",
           gidLeft.size(), gidRight.size(), gidTop.size());
}
