// Genera el audit CSV de tb_mirror_sigma_vs_x.C desde la corrida limpia de t0 (31 posiciones)
void make_audit_clean(const char* dir = "/home/rrios/ej200/results/scan_end_wrapped_2026-06-09",
                      const char* out = "/tmp/sigma_vs_x_clean_20261005/audit_clean.csv") {
  std::ofstream csv(out);
  csv << "x_mm,run,path,entries,gun_x_min,gun_x_max,status\n";
  for (int run = 0; run < 31; ++run) {
    TString path = TString::Format("%s/photon_hits_run%03d.root", dir, run);
    TFile f(path, "READ");
    auto* t = dynamic_cast<TTree*>(f.Get("sipm_hits"));
    if (f.IsZombie() || !t) { csv << ",," << path << ",0,,,FALTA\n"; continue; }
    // posición del cañón leída del propio archivo, no asumida por número de run
    const double xmin = t->GetMinimum("gun_x_mm"), xmax = t->GetMaximum("gun_x_mm");
    csv << xmin << "," << run << "," << path << "," << t->GetEntries() << ","
        << xmin << "," << xmax << "," << (xmin == xmax ? "OK" : "X_NO_UNICO") << "\n";
  }
}
