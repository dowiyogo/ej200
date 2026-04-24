// fpt_vs_n_profile_batch_slides.C
//
// Version batch que:
//   1) itera sobre photon_hits_runNNN.root,
//   2) llama a fpt_vs_n_profile.C por cada run,
//   3) guarda los canvas como PNG,
//   4) escribe un manifiesto TSV con la metadata real del scan,
//   5) genera un Beamer .tex usando fpt_manifest_to_beamer.py.
//
// La presentacion no infiere captions ni orden desde nombres de archivo:
// usa run, gun_x_mm y cara conocidos durante el loop del scan.
//
// Uso:
//   root -l -b -q fpt_vs_n_profile_batch_slides.C
//   root -l -b -q 'fpt_vs_n_profile_batch_slides.C("photon_hits_run", "fpt_png", 0, 20, 20, "presentacion_timing_detector_scan.tex")'
//
// Regenerar solo el Beamer desde el manifiesto:
//   python3 fpt_manifest_to_beamer.py fpt_png/fpt_beamer_manifest.tsv presentacion_timing_detector_scan.tex

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

struct FptSlideRow {
    int run;
    int xPos;
    std::string rootFile;
    std::string face;
    std::string faceLabel;
    std::string sipmIds;
    std::string imagePath;
    std::string caption;
};

TString fpt_slides_join_path(const char* dir, const char* name) {
    TString d(dir);
    if (d.EndsWith("/")) return d + name;
    return d + "/" + name;
}

std::string fpt_slides_tsv_escape(const std::string& text) {
    std::string out;
    out.reserve(text.size());
    for (char ch : text) {
        if (ch == '\t' || ch == '\n' || ch == '\r') out.push_back(' ');
        else out.push_back(ch);
    }
    return out;
}

// Lee gun_x_mm del primer hit del archivo. Devuelve -9999 si falla.
double fpt_slides_get_gun_x(const char* fname) {
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
void fpt_slides_watermark_pads(TCanvas* c, const TString& txt) {
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
        wm->SetTextColor(kGray);
        wm->DrawLatex(0.55, 0.50, txt);
        pad->Modified();
    }
    c->Modified();
    c->Update();
}

void fpt_slides_cleanup_profiles() {
    gROOT->GetListOfCanvases()->Delete();
    TList* dirlist = gDirectory->GetList();
    if (!dirlist) return;

    TIter next(dirlist);
    TObject* obj;
    std::vector<TObject*> toDelete;
    while ((obj = next())) {
        TString n = obj->GetName();
        if (n.BeginsWith("prof_gid")) toDelete.push_back(obj);
    }
    for (auto* o : toDelete) {
        gDirectory->Remove(o);
        delete o;
    }
}

bool fpt_slides_write_manifest(const char* manifestPath,
                               const std::vector<FptSlideRow>& rows) {
    std::ofstream out(manifestPath);
    if (!out) return false;

    out << "slide_index\trun\troot_file\tgun_x_mm\tface\tface_label\t"
        << "sipm_ids\timage_path\tcaption\n";
    for (size_t i = 0; i < rows.size(); ++i) {
        const auto& r = rows[i];
        out << (i + 1) << '\t'
            << r.run << '\t'
            << fpt_slides_tsv_escape(r.rootFile) << '\t'
            << r.xPos << '\t'
            << fpt_slides_tsv_escape(r.face) << '\t'
            << fpt_slides_tsv_escape(r.faceLabel) << '\t'
            << fpt_slides_tsv_escape(r.sipmIds) << '\t'
            << fpt_slides_tsv_escape(r.imagePath) << '\t'
            << fpt_slides_tsv_escape(r.caption) << '\n';
    }
    return true;
}

void fpt_vs_n_profile_batch_slides(const char* prefix      = "photon_hits_run",
                                   const char* outDir      = "fpt_png",
                                   int         runStart    = 0,
                                   int         runEnd      = 20,
                                   int         nMax        = 20,
                                   const char* outputTex   = "presentacion_timing_detector_scan.tex",
                                   const char* manifestTsv = "fpt_beamer_manifest.tsv") {
    if (gROOT->LoadMacro("fpt_vs_n_profile.C") < 0) {
        printf("ERROR: no pude cargar fpt_vs_n_profile.C "
               "(debe estar en el directorio actual).\n");
        return;
    }

    gSystem->mkdir(outDir, kTRUE);
    TString manifestPath = fpt_slides_join_path(outDir, manifestTsv);

    printf("Guardando PNGs en: %s/  (nMax = %d fotones)\n", outDir, nMax);
    printf("Manifest TSV : %s\n", manifestPath.Data());
    printf("Beamer TEX   : %s\n", outputTex);

    const char* faceCanvas[3] = {"cLeft", "cRight", "cTop"};
    const char* faceName  [3] = {"End-left", "End-right", "Top"};
    const char* faceLabel [3] = {"left", "right", "top"};
    const char* faceIds   [3] = {"0-7", "8-15", "16-35"};

    int nOK = 0, nFail = 0;
    std::vector<FptSlideRow> slides;

    for (int run = runStart; run <= runEnd; ++run) {
        TString fname = TString::Format("%s%03d.root", prefix, run);
        if (gSystem->AccessPathName(fname)) {
            printf("\n[run %03d] %s no existe, salteando\n", run, fname.Data());
            ++nFail;
            continue;
        }

        double gx = fpt_slides_get_gun_x(fname);
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

        TString wmText = TString::Format("#mu @ x = %+d mm", xPos);
        for (int i = 0; i < 3; ++i) {
            TCanvas* c = (TCanvas*)gROOT->GetListOfCanvases()->FindObject(faceCanvas[i]);
            if (!c) {
                printf("  [!] canvas %s ausente\n", faceCanvas[i]);
                continue;
            }

            fpt_slides_watermark_pads(c, wmText);

            TString pngName = TString::Format("muon_%+dmm_%s.png", xPos, faceLabel[i]);
            TString outName = fpt_slides_join_path(outDir, pngName);
            c->SaveAs(outName);

            std::ostringstream caption;
            caption << "Muon @ x = " << (xPos >= 0 ? "+" : "") << xPos
                    << " mm -- " << faceName[i] << " SiPMs (IDs "
                    << faceIds[i] << ")";

            FptSlideRow row;
            row.run = run;
            row.xPos = xPos;
            row.rootFile = fname.Data();
            row.face = faceName[i];
            row.faceLabel = faceLabel[i];
            row.sipmIds = faceIds[i];
            row.imagePath = outName.Data();
            row.caption = caption.str();
            slides.push_back(row);
        }

        fpt_slides_cleanup_profiles();
        ++nOK;
    }

    if (!fpt_slides_write_manifest(manifestPath.Data(), slides)) {
        printf("\nERROR: no pude escribir manifest %s\n", manifestPath.Data());
        return;
    }

    TString pyCmd = TString::Format(
        "python3 fpt_manifest_to_beamer.py \"%s\" \"%s\"",
        manifestPath.Data(), outputTex);
    int pyStatus = gSystem->Exec(pyCmd);
    if (pyStatus != 0) {
        printf("\nERROR: fallo generando Beamer con:\n  %s\n", pyCmd.Data());
        printf("Puedes regenerarlo manualmente con el mismo comando.\n");
        return;
    }

    printf("\n==============================\n");
    printf(" Resumen : %d runs OK, %d saltados\n", nOK, nFail);
    printf(" PNGs    : %s/\n", outDir);
    printf(" Manifest: %s\n", manifestPath.Data());
    printf(" Beamer  : %s\n", outputTex);
    printf(" Slides  : %zu\n", slides.size());
    printf("==============================\n");
}
