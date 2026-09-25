#include <TFile.h>
#include <TNamed.h>
#include <TTree.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::vector<std::string> split(const std::string& s) {
  std::stringstream ss(s); std::string x; std::vector<std::string> out;
  while (std::getline(ss, x, ',')) out.push_back(x);
  return out;
}
using Row = std::map<std::string,std::string>;
std::vector<Row> readCsv(const std::string& path) {
  std::ifstream in(path); if (!in) throw std::runtime_error("cannot open "+path);
  std::string line; std::getline(in,line); const auto h=split(line); std::vector<Row> rows;
  while(std::getline(in,line)){if(line.empty())continue;auto v=split(line);Row r;for(size_t i=0;i<h.size();i++)r[h[i]]=v.at(i);rows.push_back(r);} return rows;
}
double d(const Row&r,const std::string&k){return std::stod(r.at(k));}
int i(const Row&r,const std::string&k){return std::stoi(r.at(k));}
}

void summarize_diagnostics(){
  const auto rows=readCsv("sources/timing_widths.csv");
  std::map<std::string,Row> by;
  for(const auto&r:rows)by[r.at("material")+"/"+r.at("estimator")+"/"+r.at("x_mm")]=r;
  TFile out("sources/symmetry_metrics.root","RECREATE"); TTree tree("symmetry","mirrored timing-width differences");
  char material[16]={}, estimator[24]={}; int abs_x=0; double dg=0,deg=0,zg=0,dr=0,der=0,zr=0,dq=0,deq=0,zq=0;
  tree.Branch("material",material,"material/C");tree.Branch("estimator",estimator,"estimator/C");tree.Branch("abs_x_mm",&abs_x);
  tree.Branch("delta_sigma_G_ns",&dg);tree.Branch("delta_sigma_G_err_ns",&deg);tree.Branch("z_G",&zg);
  tree.Branch("delta_RMS_ns",&dr);tree.Branch("delta_RMS_err_ns",&der);tree.Branch("z_RMS",&zr);
  tree.Branch("delta_sigma68_ns",&dq);tree.Branch("delta_sigma68_err_ns",&deq);tree.Branch("z_68",&zq);
  std::ofstream csv("sources/symmetry_metrics.csv");csv<<std::setprecision(12)<<"material,estimator,abs_x_mm,delta_sigma_G_ns,delta_sigma_G_err_ns,z_G,delta_RMS_ns,delta_RMS_err_ns,z_RMS,delta_sigma68_ns,delta_sigma68_err_ns,z_68\n";
  for(auto mat:{std::string("EJ-200"),std::string("EJ-204"),std::string("EJ-230")})for(auto est:{std::string("first_pe"),std::string("mean_first_5")})for(int a:{200,500,650}){
    const auto& n=by.at(mat+"/"+est+"/-"+std::to_string(a)); const auto&p=by.at(mat+"/"+est+"/"+std::to_string(a));
    dg=d(p,"baseline_sigma_G_ns")-d(n,"baseline_sigma_G_ns");deg=std::hypot(d(p,"baseline_sigma_G_err_ns"),d(n,"baseline_sigma_G_err_ns"));zg=dg/deg;
    dr=d(p,"RMS_ns")-d(n,"RMS_ns");der=std::hypot(d(p,"RMS_bootstrap_err_ns"),d(n,"RMS_bootstrap_err_ns"));zr=dr/der;
    dq=d(p,"sigma68_ns")-d(n,"sigma68_ns");deq=std::hypot(d(p,"sigma68_bootstrap_err_ns"),d(n,"sigma68_bootstrap_err_ns"));zq=dq/deq;
    snprintf(material,sizeof(material),"%s",mat.c_str());snprintf(estimator,sizeof(estimator),"%s",est.c_str());abs_x=a;tree.Fill();
    csv<<mat<<','<<est<<','<<a<<','<<dg<<','<<deg<<','<<zg<<','<<dr<<','<<der<<','<<zr<<','<<dq<<','<<deq<<','<<zq<<'\n';
  }
  tree.Write();TNamed method("method","Differences are +|x| minus -|x|; independent-cell errors added in quadrature; RMS and sigma68 errors from 500 event bootstraps; Gaussian errors from ROOT covariance");method.Write();out.Close();
}
