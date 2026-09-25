#include "figure_common.C"
#include <TRandom3.h>
#include <TDirectory.h>
#include <TTree.h>

#include <array>
#include <numeric>

namespace {
struct Boot { double rmsErr=0,qErr=0; };
std::pair<double,double> rmsQ(std::vector<double> v) {
  double mean=std::accumulate(v.begin(),v.end(),0.0)/v.size(),ss=0;
  for(double x:v)ss+=(x-mean)*(x-mean);
  std::sort(v.begin(),v.end());
  auto q=[&](double p){double z=p*(v.size()-1),a=std::floor(z),u=z-a;size_t n=a;return n+1<v.size()?v[n]*(1-u)+v[n+1]*u:v[n];};
  return {std::sqrt(ss/v.size()),0.5*(q(.84)-q(.16))};
}
Boot bootstrap(const std::vector<double>& values,UInt_t seed,int replicas=300) {
  TRandom3 rng(seed);std::vector<double> sample(values.size()),rv,qv;rv.reserve(replicas);qv.reserve(replicas);
  for(int b=0;b<replicas;b++){for(size_t j=0;j<values.size();j++)sample[j]=values[rng.Integer(values.size())];auto w=rmsQ(sample);rv.push_back(w.first);qv.push_back(w.second);}
  auto sd=[](const std::vector<double>&v){double m=std::accumulate(v.begin(),v.end(),0.0)/v.size(),s=0;for(double x:v)s+=(x-m)*(x-m);return std::sqrt(s/(v.size()-1));};
  return {sd(rv),sd(qv)};
}
}

void build_symmetry_diagnostics(){
  using namespace v9;
  auto tab=readCsv("sources/timing_summary.csv");
  TFile events("sources/timing_events.root","READ");auto* tr=(TTree*)events.Get("timing_events");
  char cell[32]={};double l[20],r[20];tr->SetBranchAddress("cell_id",cell);tr->SetBranchAddress("t_left_ns",l);tr->SetBranchAddress("t_right_ns",r);
  std::map<std::string,std::vector<double>> values;
  for(Long64_t j=0;j<tr->GetEntries();j++){tr->GetEntry(j);values[cell].push_back(.5*(l[0]+r[0]));}
  std::map<std::string,Boot> boots;UInt_t seed=26091443;for(auto&x:values)boots[x.first]=bootstrap(x.second,seed++);

  TFile src("sources/timing_summary.root","READ"),out("sources/t0_fit_diagnostics.root","RECREATE");
  char oid[32],omat[16];int ox=0,on=0,ondf=0,ostatus=0;double os=0,ose=0,om=0,ome=0,ochi=0,olo=0,ohi=0,orms=0,ormse=0,oq=0,oqe=0;
  TTree tree("t0_fit_diagnostics","T0 fit and width diagnostics");
  tree.Branch("cell_id",oid,"cell_id/C");tree.Branch("material",omat,"material/C");tree.Branch("x_mm",&ox,"x_mm/I");tree.Branch("N",&on,"N/I");tree.Branch("sigma_G_ns",&os,"sigma_G_ns/D");tree.Branch("sigma_G_err_ns",&ose,"sigma_G_err_ns/D");tree.Branch("fit_mean_ns",&om,"fit_mean_ns/D");tree.Branch("fit_mean_err_ns",&ome,"fit_mean_err_ns/D");tree.Branch("chi2",&ochi,"chi2/D");tree.Branch("ndf",&ondf,"ndf/I");tree.Branch("fit_status",&ostatus,"fit_status/I");tree.Branch("fit_lo_ns",&olo,"fit_lo_ns/D");tree.Branch("fit_hi_ns",&ohi,"fit_hi_ns/D");tree.Branch("global_rms_ns",&orms,"global_rms_ns/D");tree.Branch("global_rms_bootstrap_err_ns",&ormse,"global_rms_bootstrap_err_ns/D");tree.Branch("qwidth_ns",&oq,"qwidth_ns/D");tree.Branch("qwidth_bootstrap_err_ns",&oqe,"qwidth_bootstrap_err_ns/D");
  std::ofstream csv("sources/t0_fit_diagnostics.csv");csv<<std::setprecision(12)<<"cell_id,material,x_mm,N,sigma_G_ns,sigma_G_err_ns,fit_mean_ns,fit_mean_err_ns,chi2,ndf,chi2_ndf,fit_status,fit_lo_ns,fit_hi_ns,global_rms_ns,global_rms_bootstrap_err_ns,qwidth_ns,qwidth_bootstrap_err_ns\n";
  std::map<std::pair<std::string,int>,std::map<std::string,std::string>> rows;
  for(auto&z:tab.rows){std::string id=z.at("cell_id"),mat=z.at("material");int x=i(z,"x_mm");rows[{mat,x}]=z;std::snprintf(oid,sizeof(oid),"%s",id.c_str());std::snprintf(omat,sizeof(omat),"%s",mat.c_str());ox=x;on=i(z,"N");os=d(z,"t0_sigma_ns");ose=d(z,"t0_sigma_err");om=d(z,"t0_fit_mean_ns");ome=d(z,"t0_fit_mean_err");ochi=d(z,"t0_chi2");ondf=i(z,"t0_ndf");ostatus=i(z,"t0_fit_status");olo=d(z,"t0_fit_lo_ns");ohi=d(z,"t0_fit_hi_ns");orms=d(z,"t0_rms_ns");ormse=boots[id].rmsErr;oq=d(z,"t0_qwidth_ns");oqe=boots[id].qErr;tree.Fill();csv<<id<<','<<mat<<','<<x<<','<<on<<','<<os<<','<<ose<<','<<om<<','<<ome<<','<<ochi<<','<<ondf<<','<<ochi/ondf<<','<<ostatus<<','<<olo<<','<<ohi<<','<<orms<<','<<ormse<<','<<oq<<','<<oqe<<'\n';out.cd();auto*dir=out.mkdir(id.c_str());dir->cd();auto*h=(TH1D*)src.Get(("h_"+id+"_t0").c_str());auto*f=(TF1*)src.Get(("fit_"+id+"_t0").c_str());if(!h||!f)throw std::runtime_error("Missing fit objects for "+id);h->Clone("t0_histogram")->Write();f->Clone("gaussian_fit")->Write();}
  out.cd();tree.Write();TNamed method("method","Exact histograms and TF1 fits used by timing_summary; qwidth=(q84-q16)/2; 300 deterministic event-bootstrap replicas for RMS/qwidth errors; bootstrap seed starts at 26091443");method.Write();

  char smat[16];int ax=0;double ds=0,dse=0,zs=0,dr=0,dre=0,zr=0,dq=0,dqe=0,zq=0;TTree sym("symmetry_diagnostics","plus-x minus minus-x width diagnostics");sym.Branch("material",smat,"material/C");sym.Branch("abs_x_mm",&ax,"abs_x_mm/I");sym.Branch("delta_sigma_G_ns",&ds,"delta_sigma_G_ns/D");sym.Branch("delta_sigma_G_err_ns",&dse,"delta_sigma_G_err_ns/D");sym.Branch("z_sigma_G",&zs,"z_sigma_G/D");sym.Branch("delta_rms_ns",&dr,"delta_rms_ns/D");sym.Branch("delta_rms_err_ns",&dre,"delta_rms_err_ns/D");sym.Branch("z_rms",&zr,"z_rms/D");sym.Branch("delta_qwidth_ns",&dq,"delta_qwidth_ns/D");sym.Branch("delta_qwidth_err_ns",&dqe,"delta_qwidth_err_ns/D");sym.Branch("z_qwidth",&zq,"z_qwidth/D");
  std::ofstream scsv("sources/symmetry_diagnostics.csv");scsv<<std::setprecision(12)<<"material,abs_x_mm,delta_sigma_G_ps,delta_sigma_G_err_ps,z_sigma_G,delta_RMS_ps,delta_RMS_err_ps,z_RMS,delta_qwidth_ps,delta_qwidth_err_ps,z_qwidth\n";
  for(auto mat:{"EJ-200","EJ-204","EJ-230"})for(int a:{200,500,650}){auto&p=rows[{mat,a}];auto&n=rows[{mat,-a}];std::string pid=p.at("cell_id"),nid=n.at("cell_id");std::snprintf(smat,sizeof(smat),"%s",mat);ax=a;ds=d(p,"t0_sigma_ns")-d(n,"t0_sigma_ns");dse=std::hypot(d(p,"t0_sigma_err"),d(n,"t0_sigma_err"));zs=ds/dse;dr=d(p,"t0_rms_ns")-d(n,"t0_rms_ns");dre=std::hypot(boots[pid].rmsErr,boots[nid].rmsErr);zr=dr/dre;dq=d(p,"t0_qwidth_ns")-d(n,"t0_qwidth_ns");dqe=std::hypot(boots[pid].qErr,boots[nid].qErr);zq=dq/dqe;sym.Fill();scsv<<mat<<','<<a<<','<<1000*ds<<','<<1000*dse<<','<<zs<<','<<1000*dr<<','<<1000*dre<<','<<zr<<','<<1000*dq<<','<<1000*dqe<<','<<zq<<'\n';}
  out.cd();sym.Write();out.Close();std::cout<<"Wrote 21 fit diagnostics and 9 mirrored-pair comparisons"<<std::endl;
}
