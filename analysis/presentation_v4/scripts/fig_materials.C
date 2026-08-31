// fig_materials.C — Material comparison END-only Vikuiti scan
// Reads scan_end_vikuiti for EJ-200, EJ-204, EJ-230
// Generates sigma_T vs x figure with external legend
{
    gROOT->SetBatch(kTRUE);
    
    TStyle *st = new TStyle("ship","");
    st->SetOptStat(0); st->SetOptFit(0);
    st->SetCanvasColor(0); st->SetPadColor(0);
    st->SetFrameFillColor(0); st->SetCanvasBorderMode(0);
    st->SetPadBorderMode(0); st->SetFrameBorderMode(0);
    for (auto ax : {"x","y","z"}) {
        st->SetTitleFont(42,ax); st->SetLabelFont(42,ax);
        st->SetTitleSize(0.055,ax); st->SetLabelSize(0.048,ax);
        st->SetTitleOffset(1.0,ax);
    }
    st->SetTextFont(42); st->SetTextSize(0.05);
    st->SetPadTopMargin(0.07); st->SetPadRightMargin(0.22);
    st->SetPadBottomMargin(0.14); st->SetPadLeftMargin(0.14);
    st->SetTitleOffset(1.15,"y");
    st->SetHistLineWidth(2);
    st->SetPadTickX(1); st->SetPadTickY(1);
    st->SetLegendBorderSize(1); st->SetLegendFont(42); st->SetLegendFillColor(0);
    st->SetMarkerSize(1.3);
    gROOT->SetStyle("ship"); gROOT->ForceStyle();

    // EJ-230 optimal numbers from CSV (phase_sparse_top_optimal.csv)
    // These are from the sigma_vs_x scan in the END-only geometry
    // We'll read the phase_sparse_top_optimal.csv for all materials at N=20
    // Actually use verification.json for the hybrid scan
    // For material comparison, use the END-only scan: scan_end_vikuiti
    
    // Read phase_ab_optimal.csv which has END-only results per x for EJ200/204/230
    const char *csvfile = "/home/rrios/ej200/analysis/optim/phase_ab_optimal.csv";
    ifstream fin(csvfile);
    if (!fin.is_open()) {
        printf("ERROR: cannot open %s\n", csvfile);
        return;
    }
    
    // Parse CSV: scint, x, k_opt, sig_top_opt, m_opt, sig_end_opt, sig_comb
    string line;
    getline(fin, line); // header
    printf("Header: %s\n", line.c_str());
    
    map<string, vector<double>> x_vals, sig_end, sig_top;
    while (getline(fin, line)) {
        if (line.empty()) continue;
        stringstream ss(line);
        string sc; double x, k, st, m, se, sc2;
        char comma;
        getline(ss, sc, ',');
        ss >> x >> comma >> k >> comma >> st >> comma >> m >> comma >> se;
        x_vals[sc].push_back(x);
        sig_end[sc].push_back(se);
        sig_top[sc].push_back(st);
    }
    fin.close();
    
    for (auto &[sc, xv] : x_vals) {
        printf("Scint %s: %zu points\n", sc.c_str(), xv.size());
        for (size_t i=0; i<xv.size(); i++)
            printf("  x=%.0f sig_end=%.2f sig_top=%.2f\n", xv[i], sig_end[sc][i], sig_top[sc][i]);
    }
}
