#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr int kOrderMax = 20;

struct InputCell {
  std::string id;
  std::string name;
  int material = -1;
  int x = 0;
  int n = 0;
  std::string path;
  std::string recordedHash;
};

std::vector<std::string> splitTab(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string field;
  while (std::getline(ss, field, '\t')) out.push_back(field);
  return out;
}

std::vector<InputCell> loadInputs(const char* path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error(std::string("Cannot open ") + path);
  std::string line;
  std::getline(input, line);
  std::vector<InputCell> cells;
  while (std::getline(input, line)) {
    auto f = splitTab(line);
    if (f.size() != 6) throw std::runtime_error("Malformed input row: " + line);
    InputCell c;
    c.id = f[0]; c.name = f[1]; c.x = std::stoi(f[2]); c.n = std::stoi(f[3]);
    c.path = f[4]; c.recordedHash = f[5];
    c.material = c.name == "EJ-200" ? 0 : c.name == "EJ-204" ? 1 : 2;
    cells.push_back(c);
  }
  return cells;
}

void keepFirst(std::vector<double>& values) {
  if (values.size() > kOrderMax) {
    std::nth_element(values.begin(), values.begin() + kOrderMax, values.end());
    values.resize(kOrderMax);
  }
  std::sort(values.begin(), values.end());
}
}

void build_timing_dataset() {
  const char* inventory = "sources/input_roots.tsv";
  const auto cells = loadInputs(inventory);
  if (cells.size() != 21) throw std::runtime_error("Expected 21 cells");

  TFile output("sources/timing_events.root", "RECREATE");
  if (output.IsZombie()) throw std::runtime_error("Cannot create output ROOT");
  TTree events("timing_events", "Event-level light and detected-photon order statistics");

  char cellId[32] = {};
  int material = -1, xMm = 0, eventId = -1;
  int npeLeft = 0, npeRight = 0, npeTop = 0, producedScint = 0;
  double edepMeV = 0.0, topFirstNs = 0.0;
  double tLeftNs[kOrderMax] = {}, tRightNs[kOrderMax] = {};
  events.Branch("cell_id", cellId, "cell_id/C");
  events.Branch("material", &material, "material/I");
  events.Branch("x_mm", &xMm, "x_mm/I");
  events.Branch("event_id", &eventId, "event_id/I");
  events.Branch("npe_left", &npeLeft, "npe_left/I");
  events.Branch("npe_right", &npeRight, "npe_right/I");
  events.Branch("npe_top", &npeTop, "npe_top/I");
  events.Branch("produced_scint", &producedScint, "produced_scint/I");
  events.Branch("edep_MeV", &edepMeV, "edep_MeV/D");
  events.Branch("t_left_ns", tLeftNs, "t_left_ns[20]/D");
  events.Branch("t_right_ns", tRightNs, "t_right_ns[20]/D");
  events.Branch("t_top_first_ns", &topFirstNs, "t_top_first_ns/D");

  std::ofstream audit("sources/timing_extraction.csv");
  audit << "cell_id,material,x_mm,N,hit_entries,valid_events,count_mismatches,"
           "mean_npe_left,mean_npe_right,mean_npe_top,mean_edep_MeV,mean_produced_scint\n";
  audit << std::setprecision(12);

  for (const auto& cell : cells) {
    std::cout << "Reading " << cell.id << " from " << cell.path << std::endl;
    TFile input(cell.path.c_str(), "READ");
    if (input.IsZombie()) throw std::runtime_error("Zombie ROOT: " + cell.path);
    auto* obs = dynamic_cast<TTree*>(input.Get("event_observables"));
    auto* hits = dynamic_cast<TTree*>(input.Get("sipm_hits"));
    if (!obs || !hits) throw std::runtime_error("Required tree absent: " + cell.path);
    if (obs->GetEntries() != cell.n) throw std::runtime_error("Wrong event count: " + cell.id);

    std::vector<int> nL(cell.n, -1), nR(cell.n, -1), nT(cell.n, -1), nProd(cell.n, -1);
    std::vector<double> eDep(cell.n, std::numeric_limits<double>::quiet_NaN());
    int oe = -1, ol = 0, oright = 0, ot = 0, op = 0;
    double oedep = 0;
    obs->SetBranchStatus("*", 0);
    for (const char* b : {"event_id", "detected_left", "detected_right", "detected_top",
                          "produced_scint", "edep_total_MeV"}) obs->SetBranchStatus(b, 1);
    obs->SetBranchAddress("event_id", &oe);
    obs->SetBranchAddress("detected_left", &ol);
    obs->SetBranchAddress("detected_right", &oright);
    obs->SetBranchAddress("detected_top", &ot);
    obs->SetBranchAddress("produced_scint", &op);
    obs->SetBranchAddress("edep_total_MeV", &oedep);
    for (Long64_t i = 0; i < obs->GetEntries(); ++i) {
      obs->GetEntry(i);
      if (oe < 0 || oe >= cell.n || nL[oe] >= 0) throw std::runtime_error("Bad event ids: " + cell.id);
      nL[oe] = ol; nR[oe] = oright; nT[oe] = ot; nProd[oe] = op; eDep[oe] = oedep;
    }

    std::vector<std::vector<double>> left(cell.n), right(cell.n);
    std::vector<double> top(cell.n, std::numeric_limits<double>::infinity());
    for (int e = 0; e < cell.n; ++e) {
      if (nL[e] < 0 || nR[e] < 0 || nT[e] < 0) throw std::runtime_error("Missing event: " + cell.id);
      left[e].reserve(nL[e]); right[e].reserve(nR[e]);
    }

    int he = -1, gid = -1;
    double timeNs = 0.0;
    hits->SetBranchStatus("*", 0);
    for (const char* b : {"event_id", "global_id", "time_ns"}) hits->SetBranchStatus(b, 1);
    hits->SetBranchAddress("event_id", &he);
    hits->SetBranchAddress("global_id", &gid);
    hits->SetBranchAddress("time_ns", &timeNs);
    hits->SetCacheSize(256LL * 1024 * 1024);
    hits->AddBranchToCache("event_id", true);
    hits->AddBranchToCache("global_id", true);
    hits->AddBranchToCache("time_ns", true);
    for (Long64_t i = 0; i < hits->GetEntries(); ++i) {
      hits->GetEntry(i);
      if (he < 0 || he >= cell.n || !std::isfinite(timeNs))
        throw std::runtime_error("Invalid hit: " + cell.id);
      if (gid >= 0 && gid < 8) left[he].push_back(timeNs);
      else if (gid >= 8 && gid < 16) right[he].push_back(timeNs);
      else if (gid >= 16 && gid < 86 && timeNs < top[he]) top[he] = timeNs;
      else if (gid < 0 || gid >= 86) throw std::runtime_error("Invalid sensor id: " + cell.id);
    }

    long mismatches = 0, valid = 0;
    double sumL = 0, sumR = 0, sumT = 0, sumE = 0, sumP = 0;
    for (int e = 0; e < cell.n; ++e) {
      if (int(left[e].size()) != nL[e] || int(right[e].size()) != nR[e]) ++mismatches;
      keepFirst(left[e]); keepFirst(right[e]);
      if (left[e].size() < kOrderMax || right[e].size() < kOrderMax || !std::isfinite(top[e]))
        throw std::runtime_error("Insufficient timing hits: " + cell.id);
      std::snprintf(cellId, sizeof(cellId), "%s", cell.id.c_str());
      material = cell.material; xMm = cell.x; eventId = e;
      npeLeft = nL[e]; npeRight = nR[e]; npeTop = nT[e]; producedScint = nProd[e];
      edepMeV = eDep[e]; topFirstNs = top[e];
      for (int k = 0; k < kOrderMax; ++k) {
        tLeftNs[k] = left[e][k]; tRightNs[k] = right[e][k];
      }
      events.Fill(); ++valid;
      sumL += nL[e]; sumR += nR[e]; sumT += nT[e]; sumE += eDep[e]; sumP += nProd[e];
    }
    audit << cell.id << ',' << cell.name << ',' << cell.x << ',' << cell.n << ','
          << hits->GetEntries() << ',' << valid << ',' << mismatches << ','
          << sumL/cell.n << ',' << sumR/cell.n << ',' << sumT/cell.n << ','
          << sumE/cell.n << ',' << sumP/cell.n << '\n';
  }

  output.cd();
  events.Write();
  TNamed source("source", "21 validated current-transport ROOT files listed in sources/input_roots.tsv");
  TNamed definition("definition", "END left global_id 0..7; END right 8..15; TOP 16..85; k=1..20 sorted detected-PE time_ns per END side; SPTR=0; no electronics");
  source.Write(); definition.Write();
  output.Close();
  std::cout << "Wrote " << events.GetEntries() << " event rows" << std::endl;
}
