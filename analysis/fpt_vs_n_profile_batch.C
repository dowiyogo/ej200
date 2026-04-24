// fpt_vs_n_profile_batch.C
//
// Extension batch de fpt_vs_n_profile.C: itera sobre los archivos
// photon_hits_run000.root ... photon_hits_run020.root, llama a la funcion
// fpt_vs_n_profile por cada uno, y guarda los tres canvas como PNG con
// nombres del tipo:
//
//   muon_-650mm_left.png
//   muon_-650mm_right.png
//   muon_-650mm_top.png
//
// El formato visual (titulo con gid + posicion fisica) lo impone
// fpt_vs_n_profile.C y se conserva tal cual.
//
// Watermark:
//   El watermark de fpt_vs_n_profile.C usa SetTextColorAlpha(kGray, 0.12)
//   que en modo batch se pierde al convertir a PNG porque el canal alpha
//   no se preserva correctamente en el backend raster. Para garantizar
//   que el watermark quede visible en los PNG, el batch dibuja aqui un
//   watermark adicional con color gris solido (sin alpha) antes de salvar.
//
// Uso:
//   root -l -b -q fpt_vs_n_profile_batch.C
//   root -l -b -q 'fpt_vs_n_profile_batch.C("photon_hits_run", "fpt_png", 0, 20, 20)'

// Lee gun_x_mm del primer hit del archivo.  Devuelve -9999 si falla.
double fpt_get_gun_x(const char* fname) {
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) { if (f) delete f; return -9999.0; }
    TTree* t = (TTree*)f->Get("sipm_hits");
    if (!t || t->GetEntries() == 0) { delete f; return -9999.0; }
    Double_t gx = 0.0;
    t->SetBranchAddress("gun_x_mm", &gx);
    t->GetEntry(0);
    delete f;
    return gx;
}

// Dibuja un watermark robusto para PNG en cada subpad del canvas.
// Usa color gris indice 17 (kGray), sin alpha, que renderiza bien en
// raster export.  kGray-1 es aun mas tenue; probar kGray o kGray+1 si
// queda muy tenue.
void fpt_watermark_pads(TCanvas* c, const TString& txt) {
    if (!c) return;
    TIter next(c->GetListOfPrimitives());
    TObject* obj;
    while ((obj = next())) {
        if (!obj->InheritsFrom("TPad")) continue;
        TPad* pad = (TPad*)obj;
        pad->cd();
        TLatex* wm = new TLatex();
        wm->SetNDC();
        wm->SetTextAlign(22);
        wm->SetTextSize(0.09);
        wm->SetTextAngle(20);
        wm->SetTextColor(kGray);    // solido, sin alpha -> visible en PNG
        wm->DrawLatex(0.55, 0.50, txt);
        pad->Modified();
    }
    c->Modified();
    c->Update();
}

void fpt_vs_n_profile_batch(const char* prefix   = "photon_hits_run",
                            const char* outDir   = "fpt_png",
                            int         runStart = 0,
                            int         runEnd   = 20,
                            int         nMax     = 20) {
    if (gROOT->LoadMacro("fpt_vs_n_profile.C") < 0) {
        printf("ERROR: no pude cargar fpt_vs_n_profile.C "
               "(debe estar en el directorio actual).\n");
        return;
    }

    gSystem->mkdir(outDir, kTRUE);
    printf("Guardando PNGs en: %s/  (nMax = %d fotones)\n", outDir, nMax);

    const char* faceCanvas[3] = {"cLeft", "cRight", "cTop"};
    const char* faceLabel [3] = {"left",  "right",  "top"};

    int nOK = 0, nFail = 0;
    for (int run = runStart; run <= runEnd; ++run) {
        TString fname = TString::Format("%s%03d.root", prefix, run);
        if (gSystem->AccessPathName(fname)) {
            printf("\n[run %03d] %s no existe, salteando\n", run, fname.Data());
            ++nFail;
            continue;
        }

        double gx = fpt_get_gun_x(fname);
        if (gx < -9000.0) {
            printf("\n[run %03d] no pude leer gun_x_mm, salteando\n", run);
            ++nFail;
            continue;
        }
        int xPos = (int)std::round(gx);
        printf("\n=== Run %03d: %s, x = %+d mm ===\n", run, fname.Data(), xPos);

        TString cmd = TString::Format(
            "fpt_vs_n_profile(\"%s\", %d)", fname.Data(), nMax);
        gROOT->ProcessLine(cmd);

        // Watermark robusto + SaveAs para cada canvas
        TString wmText = TString::Format("#mu @ x = %+d mm", xPos);
        for (int i = 0; i < 3; ++i) {
            TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(faceCanvas[i]);
            if (!c) { printf("  [!] canvas %s ausente\n", faceCanvas[i]); continue; }

            fpt_watermark_pads(c, wmText);

            TString outName = TString::Format("%s/muon_%+dmm_%s.png",
                outDir, xPos, faceLabel[i]);
            c->SaveAs(outName);
        }

        // Limpieza entre iteraciones
        gROOT->GetListOfCanvases()->Delete();
        TList* dirlist = gDirectory->GetList();
        if (dirlist) {
            TIter next(dirlist);
            TObject* obj;
            std::vector<TObject*> toDelete;
            while ((obj = next())) {
                TString n = obj->GetName();
                if (n.BeginsWith("prof_gid")) toDelete.push_back(obj);
            }
            for (auto* o : toDelete) { gDirectory->Remove(o); delete o; }
        }

        ++nOK;
    }

    printf("\n==============================\n");
    printf(" Resumen: %d runs OK, %d saltados\n", nOK, nFail);
    printf(" Salida : %s/\n", outDir);
    printf("==============================\n");
}
