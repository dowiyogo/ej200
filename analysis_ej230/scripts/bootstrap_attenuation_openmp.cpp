/**
 * bootstrap_attenuation_openmp.cpp
 *
 * C++/OpenMP pairs bootstrap for attenuation-length fits on
 * End-only + Mylar EJ-230 data (ported from EJ-204 mirror analysis).
 *
 * Models implemented (matching analysis_endonly_mylar.py):
 *   M1  N0 * exp(-d/lambda)
 *   M2  As*exp(-d/ls) + Al*exp(-d/ll)   with  ll = ls + exp(log_delta) > ls
 *   M3  A*exp(-d/lambda) + C
 *   M4  N0*exp(-d/lambda)               fitted on d > D_TAIL_CM only
 *
 * Compile:
 *   g++ -O2 -fopenmp -std=c++17 bootstrap_attenuation_openmp.cpp \
 *       -o bootstrap_attenuation -lm
 *
 * Run modes:
 *   ./bootstrap_attenuation --quick     (n_boot=50,  seed=42)
 *   ./bootstrap_attenuation --standard  (n_boot=200, seed=42)   [default]
 *   ./bootstrap_attenuation --full      (n_boot=1000,seed=42)
 *
 * Input: attenuation_curve.csv   (columns: x_mm, npe_left, npe_right, ...)
 * Output: bootstrap_results_openmp.csv  +  bootstrap_summary_openmp.csv
 *
 * Uses Levenberg-Marquardt via an internal minimal implementation.
 * No external ROOT/GSL dependency.
 */

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
static constexpr double D_TAIL_CM   = 40.0;   // M4: fit only d > this value
static constexpr int    MAX_ITER_LM = 2000;   // Levenberg-Marquardt iterations
static constexpr double TOL_LM      = 1e-10;  // convergence tolerance

// ---------------------------------------------------------------------------
// Data structures
// ---------------------------------------------------------------------------
struct DataPoint {
    double d_cm;   // distance from near end (cm)
    double npe;    // mean NPE
    double err;    // SEM on NPE (>0)
};

struct FitResult {
    // M1
    double m1_lam;    double m1_chi2ndf;   bool m1_ok;
    // M2
    double m2_ls;     double m2_ll;
    double m2_As;     double m2_Al;
    double m2_chi2ndf; bool m2_ok;
    // M3
    double m3_lam;    double m3_C;
    double m3_chi2ndf; bool m3_ok;
    // M4 (tail)
    double m4_lam;    double m4_chi2ndf;   bool m4_ok;
};

// ---------------------------------------------------------------------------
// Minimal Levenberg-Marquardt (generic, npar <= 5)
// ---------------------------------------------------------------------------
// Model: functor f(d, p[]) -> double
// Minimises chi2 = sum_i ((y_i - f(d_i, p)) / err_i)^2
// ---------------------------------------------------------------------------
template<typename Model>
static bool lm_fit(
    const std::vector<DataPoint>& data,
    std::vector<double>& p,          // in: initial guess; out: fitted values
    double& chi2ndf,
    const Model& model,
    const std::vector<DataPoint>* subset = nullptr  // nullptr = use all data
) {
    const auto& pts = subset ? *subset : data;
    int n = (int)pts.size();
    int k = (int)p.size();
    if (n <= k) { chi2ndf = std::numeric_limits<double>::quiet_NaN(); return false; }

    auto chi2_fn = [&](const std::vector<double>& par) {
        double s = 0;
        for (const auto& pt : pts) {
            double r = (pt.npe - model(pt.d_cm, par)) / pt.err;
            s += r * r;
        }
        return s;
    };

    double lam = 1e-3;
    double chi2_cur = chi2_fn(p);

    for (int iter = 0; iter < MAX_ITER_LM; ++iter) {
        // Numerical Jacobian J (n × k), eps per parameter
        std::vector<std::vector<double>> J(n, std::vector<double>(k));
        for (int j = 0; j < k; ++j) {
            double eps = std::max(1e-6 * std::abs(p[j]), 1e-8);
            std::vector<double> pp = p;
            pp[j] += eps;
            for (int i = 0; i < n; ++i)
                J[i][j] = (model(pts[i].d_cm, pp) - model(pts[i].d_cm, p))
                           / (eps * pts[i].err);
        }
        // Residuals
        std::vector<double> r(n);
        for (int i = 0; i < n; ++i)
            r[i] = (pts[i].npe - model(pts[i].d_cm, p)) / pts[i].err;

        // JtJ + lambda*diag, Jt*r
        std::vector<std::vector<double>> A(k, std::vector<double>(k, 0));
        std::vector<double> b(k, 0);
        for (int j = 0; j < k; ++j) {
            for (int l = 0; l < k; ++l)
                for (int i = 0; i < n; ++i)
                    A[j][l] += J[i][j] * J[i][l];
            A[j][j] *= (1.0 + lam);
            for (int i = 0; i < n; ++i)
                b[j] += J[i][j] * r[i];
        }

        // Solve A*dp = b via Gaussian elimination (k <= 5, tiny matrix)
        std::vector<double> dp(k, 0);
        for (int pivot = 0; pivot < k; ++pivot) {
            int best = pivot;
            for (int row = pivot + 1; row < k; ++row)
                if (std::abs(A[row][pivot]) > std::abs(A[best][pivot])) best = row;
            std::swap(A[pivot], A[best]); std::swap(b[pivot], b[best]);
            if (std::abs(A[pivot][pivot]) < 1e-30) return false;
            for (int row = pivot + 1; row < k; ++row) {
                double f = A[row][pivot] / A[pivot][pivot];
                for (int col = pivot; col < k; ++col) A[row][col] -= f * A[pivot][col];
                b[row] -= f * b[pivot];
            }
        }
        for (int i = k - 1; i >= 0; --i) {
            dp[i] = b[i];
            for (int j = i + 1; j < k; ++j) dp[i] -= A[i][j] * dp[j];
            dp[i] /= A[i][i];
        }

        std::vector<double> p_new = p;
        for (int j = 0; j < k; ++j) p_new[j] += dp[j];

        double chi2_new = chi2_fn(p_new);
        if (chi2_new < chi2_cur) {
            p = p_new;
            chi2_cur = chi2_new;
            lam *= 0.5;
            // Check convergence
            double step = 0;
            for (double d : dp) step += d * d;
            if (std::sqrt(step) < TOL_LM) break;
        } else {
            lam *= 2.0;
        }
    }

    chi2ndf = chi2_cur / (n - k);
    return std::isfinite(chi2_cur) && chi2_cur > 0;
}

// ---------------------------------------------------------------------------
// Model functions
// ---------------------------------------------------------------------------
static double m1_model(double d, const std::vector<double>& p) {
    // p = {N0, lambda}
    return p[0] * std::exp(-d / p[1]);
}
static double m2_model(double d, const std::vector<double>& p) {
    // p = {As, ls, Al, log_delta};  ll = ls + exp(log_delta)
    double ll = p[1] + std::exp(p[3]);
    return p[0] * std::exp(-d / p[1]) + p[2] * std::exp(-d / ll);
}
static double m3_model(double d, const std::vector<double>& p) {
    // p = {A, lambda, C}
    return p[0] * std::exp(-d / p[1]) + p[2];
}
// M4 uses same M1 model but on a subset

// ---------------------------------------------------------------------------
// Fit a single dataset with all 4 models
// ---------------------------------------------------------------------------
static FitResult fit_dataset(const std::vector<DataPoint>& data) {
    FitResult res{};

    // M1 — single exponential
    {
        std::vector<double> p = {data[0].npe, 30.0};
        res.m1_ok = lm_fit(data, p, res.m1_chi2ndf, m1_model);
        res.m1_lam = p[1];
        if (res.m1_lam <= 0 || res.m1_chi2ndf < 0) res.m1_ok = false;
    }

    // M2 — constrained double exponential
    {
        // Initial: short ~10 cm, long ~40 cm
        // log_delta = log(ll - ls) = log(30)
        std::vector<double> p = {2000.0, 10.0, 400.0, std::log(30.0)};
        res.m2_ok = lm_fit(data, p, res.m2_chi2ndf, m2_model);
        res.m2_As = p[0];  res.m2_ls = p[1];
        res.m2_Al = p[2];  res.m2_ll = p[1] + std::exp(p[3]);
        if (res.m2_ls <= 0 || res.m2_ll <= res.m2_ls ||
            res.m2_As <= 0 || res.m2_Al <= 0) res.m2_ok = false;
    }

    // M3 — exp + constant floor
    {
        std::vector<double> p = {data[0].npe, 20.0, 5.0};
        res.m3_ok = lm_fit(data, p, res.m3_chi2ndf, m3_model);
        res.m3_lam = p[1];  res.m3_C = p[2];
        if (res.m3_lam <= 0) res.m3_ok = false;
    }

    // M4 — tail slope only (d > D_TAIL_CM)
    {
        std::vector<DataPoint> tail;
        for (const auto& pt : data)
            if (pt.d_cm > D_TAIL_CM) tail.push_back(pt);
        if ((int)tail.size() > 2) {
            std::vector<double> p = {100.0, 40.0};
            res.m4_ok = lm_fit(data, p, res.m4_chi2ndf, m1_model, &tail);
            res.m4_lam = p[1];
            if (res.m4_lam <= 0) res.m4_ok = false;
        }
    }

    return res;
}

// ---------------------------------------------------------------------------
// CSV I/O
// ---------------------------------------------------------------------------
static std::vector<DataPoint> load_csv(const std::string& path,
                                       const std::string& col_d,
                                       const std::string& col_npe,
                                       const std::string& col_err) {
    std::ifstream f(path);
    if (!f) { std::cerr << "Cannot open " << path << "\n"; std::exit(1); }

    std::string line;
    std::getline(f, line);   // header
    std::vector<std::string> headers;
    {
        std::istringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) headers.push_back(tok);
    }
    auto col_idx = [&](const std::string& name) -> int {
        for (int i = 0; i < (int)headers.size(); ++i)
            if (headers[i] == name) return i;
        std::cerr << "Column '" << name << "' not found\n"; std::exit(1);
    };
    int id  = col_idx(col_d);
    int inp = col_idx(col_npe);
    int ier = col_idx(col_err);

    std::vector<DataPoint> pts;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        std::vector<std::string> toks;
        std::istringstream ss(line);
        std::string tok;
        while (std::getline(ss, tok, ',')) toks.push_back(tok);
        double d_cm = std::stod(toks[id]);
        if (d_cm < 0) d_cm = -d_cm;  // convert signed x_mm to |distance| later
        pts.push_back({d_cm, std::stod(toks[inp]), std::stod(toks[ier])});
    }
    // Sort by ascending distance
    std::sort(pts.begin(), pts.end(),
              [](const DataPoint& a, const DataPoint& b){ return a.d_cm < b.d_cm; });
    return pts;
}

static void write_replica_csv(const std::string& path,
                              const std::vector<FitResult>& results,
                              int n_boot) {
    std::ofstream f(path);
    f << "replica,m1_lam,m1_chi2ndf,m1_ok,"
         "m2_ls,m2_ll,m2_As,m2_Al,m2_chi2ndf,m2_ok,"
         "m3_lam,m3_C,m3_chi2ndf,m3_ok,"
         "m4_lam,m4_chi2ndf,m4_ok\n";
    for (int i = 0; i < n_boot; ++i) {
        const auto& r = results[i];
        f << i+1 << ","
          << r.m1_lam << "," << r.m1_chi2ndf << "," << r.m1_ok << ","
          << r.m2_ls  << "," << r.m2_ll << "," << r.m2_As << "," << r.m2_Al << ","
          << r.m2_chi2ndf << "," << r.m2_ok << ","
          << r.m3_lam << "," << r.m3_C << "," << r.m3_chi2ndf << "," << r.m3_ok << ","
          << r.m4_lam << "," << r.m4_chi2ndf << "," << r.m4_ok << "\n";
    }
    std::cout << "Written: " << path << "\n";
}

struct Summary { double median, p16, p84, mean, std_dev; int n_ok; };

static Summary summarise(const std::vector<FitResult>& results, int n_boot,
                         std::function<std::pair<double,bool>(const FitResult&)> extractor) {
    std::vector<double> vals;
    for (int i = 0; i < n_boot; ++i) {
        auto [v, ok] = extractor(results[i]);
        if (ok && std::isfinite(v)) vals.push_back(v);
    }
    Summary s{};
    s.n_ok = (int)vals.size();
    if (vals.empty()) return s;
    std::sort(vals.begin(), vals.end());
    s.median = vals[vals.size()/2];
    s.p16    = vals[std::max(0, (int)(0.16 * vals.size()))];
    s.p84    = vals[std::min((int)vals.size()-1, (int)(0.84 * vals.size()))];
    double sum = 0, sum2 = 0;
    for (double v : vals) { sum += v; sum2 += v*v; }
    s.mean    = sum / vals.size();
    s.std_dev = std::sqrt(sum2/vals.size() - s.mean*s.mean);
    return s;
}

static void write_summary_csv(const std::string& path,
                              const std::vector<FitResult>& results,
                              int n_boot) {
    std::ofstream f(path);
    f << "parameter,n_ok,median,p16,p84,mean,std\n";
    auto write_row = [&](const std::string& name,
                         std::function<std::pair<double,bool>(const FitResult&)> ex) {
        auto s = summarise(results, n_boot, ex);
        f << name << "," << s.n_ok << "," << s.median << ","
          << s.p16 << "," << s.p84 << "," << s.mean << "," << s.std_dev << "\n";
    };
    write_row("M1_lam_cm",   [](const FitResult& r){ return std::make_pair(r.m1_lam, r.m1_ok); });
    write_row("M2_lam_s_cm", [](const FitResult& r){ return std::make_pair(r.m2_ls,  r.m2_ok); });
    write_row("M2_lam_l_cm", [](const FitResult& r){ return std::make_pair(r.m2_ll,  r.m2_ok); });
    write_row("M2_As",       [](const FitResult& r){ return std::make_pair(r.m2_As,  r.m2_ok); });
    write_row("M2_Al",       [](const FitResult& r){ return std::make_pair(r.m2_Al,  r.m2_ok); });
    write_row("M3_lam_cm",   [](const FitResult& r){ return std::make_pair(r.m3_lam, r.m3_ok); });
    write_row("M3_C_PE",     [](const FitResult& r){ return std::make_pair(r.m3_C,   r.m3_ok); });
    write_row("M4_lam_cm",   [](const FitResult& r){ return std::make_pair(r.m4_lam, r.m4_ok); });
    std::cout << "Written: " << path << "\n";
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    // --- Parse args ---
    int n_boot = 200;
    uint64_t global_seed = 230123;  // EJ-230 reproducible seed
    std::string input_csv  = "attenuation_curve.csv";
    std::string col_d      = "d_cm";
    std::string col_npe    = "npe_left";
    std::string col_err    = "npe_left_err";
    std::string out_replicas = "bootstrap_results_openmp.csv";
    std::string out_summary  = "bootstrap_summary_openmp.csv";

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "--quick")    { n_boot = 50; }
        else if (a == "--standard") { n_boot = 200; }
        else if (a == "--full")     { n_boot = 1000; }
        else if (a == "--seed"  && i+1 < argc) { global_seed = std::stoull(argv[++i]); }
        else if (a == "--input" && i+1 < argc) { input_csv = argv[++i]; }
        else if (a == "--col-d"   && i+1 < argc) { col_d   = argv[++i]; }
        else if (a == "--col-npe" && i+1 < argc) { col_npe = argv[++i]; }
        else if (a == "--col-err" && i+1 < argc) { col_err = argv[++i]; }
        else if (a == "--out-replicas" && i+1 < argc) { out_replicas = argv[++i]; }
        else if (a == "--out-summary"  && i+1 < argc) { out_summary  = argv[++i]; }
        else {
            std::cerr << "Unknown argument: " << a << "\n";
            return 1;
        }
    }

    // --- Load data ---
    auto data = load_csv(input_csv, col_d, col_npe, col_err);
    int n_pts = (int)data.size();
    std::cout << "Loaded " << n_pts << " data points from " << input_csv << "\n";

    // --- Fit nominal (on full dataset, single thread) ---
    auto t0_nom = std::chrono::steady_clock::now();
    FitResult nominal = fit_dataset(data);
    auto t1_nom = std::chrono::steady_clock::now();
    double ms_nom = std::chrono::duration<double,std::milli>(t1_nom - t0_nom).count();

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "\nNominal fit (" << ms_nom << " ms):\n";
    std::cout << "  M1  lambda = " << nominal.m1_lam << " cm   chi2/ndf = " << nominal.m1_chi2ndf << "\n";
    if (nominal.m2_ok)
        std::cout << "  M2  ls=" << nominal.m2_ls << " ll=" << nominal.m2_ll
                  << " chi2/ndf=" << nominal.m2_chi2ndf << "\n";
    if (nominal.m3_ok)
        std::cout << "  M3  lam=" << nominal.m3_lam << " C=" << nominal.m3_C
                  << " chi2/ndf=" << nominal.m3_chi2ndf << "\n";
    if (nominal.m4_ok)
        std::cout << "  M4  lam_tail=" << nominal.m4_lam
                  << " chi2/ndf=" << nominal.m4_chi2ndf << "\n";

    // --- Bootstrap (parallel) ---
    std::cout << "\nRunning " << n_boot << " bootstrap replicas";
#ifdef _OPENMP
    int n_threads = omp_get_max_threads();
    std::cout << " on " << n_threads << " OpenMP threads";
#endif
    std::cout << "...\n";

    // Pre-generate all replica indices (thread-safe seeding)
    // Seed per replica: global_seed + replica index (reproducible regardless of thread count)
    std::vector<FitResult> results(n_boot);
    std::vector<std::vector<DataPoint>> replicas(n_boot);
    for (int b = 0; b < n_boot; ++b) {
        std::mt19937_64 rng(global_seed + (uint64_t)b);
        std::uniform_int_distribution<int> dist(0, n_pts - 1);
        replicas[b].resize(n_pts);
        for (int i = 0; i < n_pts; ++i)
            replicas[b][i] = data[dist(rng)];
    }

    auto t0 = std::chrono::steady_clock::now();

    #pragma omp parallel for schedule(dynamic, 4)
    for (int b = 0; b < n_boot; ++b) {
        results[b] = fit_dataset(replicas[b]);
    }

    auto t1 = std::chrono::steady_clock::now();
    double ms_total = std::chrono::duration<double,std::milli>(t1 - t0).count();

    // Count M2 failures
    int m2_fail = 0;
    for (int b = 0; b < n_boot; ++b)
        if (!results[b].m2_ok) ++m2_fail;
    std::cout << "Bootstrap complete: " << ms_total << " ms  ("
              << ms_total / n_boot << " ms/replica)\n";
    std::cout << "M2 failure fraction: " << m2_fail << "/" << n_boot
              << " = " << std::setprecision(3) << (100.0 * m2_fail / n_boot) << "%\n";

    // --- Write outputs ---
    write_replica_csv(out_replicas, results, n_boot);
    write_summary_csv(out_summary, results, n_boot);

    return 0;
}
