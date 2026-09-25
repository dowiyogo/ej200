#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TMath.h>
#include <TNamed.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
constexpr int K = 20;
constexpr double cMmNs = 299.792458;
constexpr double refractiveIndex = 1.58;

struct Event {
  int npeL, npeR, npeT, produced;
  double edep, top;
  std::array<double,K> left, right;
};

struct FitResult {
  double fitMean=0, fitMeanErr=0, sigma=0, sigmaErr=0, chi2=0;
  int ndf=0, status=-1;
  double rangeLo=0, rangeHi=0, globalMean=0, globalRms=0, globalRmsErr=0, qwidth=0;
};

struct Metric {
  std::string id;
  int material=0, x=0, n=0;
  double npeL=0,npeLse=0,npeR=0,npeRse=0,npeEnd=0,npeEndSe=0,npeT=0,npeTse=0,npeAll=0,npeAllSe=0;
  double edep=0,edepSe=0,produced=0,producedSe=0,eff=0,effSe=0;
  double meanL=0,meanLse=0,meanR=0,meanRse=0,meanDt=0,meanDtSe=0;
  FitResult t0,dt,top;
  double sigmaX=0,sigmaXErr=0;
};

struct MaterialFit {
  double intercept=0,interceptErr=0,slope=0,slopeErr=0,chi2=0;
  int ndf=0,status=-1;
  double veff=0,veffErr=0,cn=0,betaDeg=0,betaErrDeg=0;
};

std::pair<double,double> meanSem(const std::vector<double>& v) {
  if (v.size()<2) throw std::runtime_error("Too few values");
  double mean=std::accumulate(v.begin(),v.end(),0.0)/v.size(), ss=0;
  for(double x:v) ss+=(x-mean)*(x-mean);
  return {mean,std::sqrt(ss/(v.size()-1))/std::sqrt(double(v.size()))};
}

double quantile(std::vector<double> v, double p) {
  if(v.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(v.begin(),v.end());
  double z=p*(v.size()-1), f=std::floor(z), d=z-f;
  size_t i=size_t(f); return i+1<v.size()?v[i]*(1-d)+v[i+1]*d:v[i];
}

FitResult fitCentral(const std::vector<double>& v, const std::string& name, TFile& out) {
  FitResult r; r.globalMean=meanSem(v).first;
  double ss=0; for(double x:v)ss+=(x-r.globalMean)*(x-r.globalMean);
  r.globalRms=std::sqrt(ss/v.size());
  r.globalRmsErr=r.globalRms/std::sqrt(2.0*(v.size()-1));
  const double q005=quantile(v,0.005), q16=quantile(v,0.16), q50=quantile(v,0.50), q84=quantile(v,0.84), q995=quantile(v,0.995);
  const double robust=0.5*(q84-q16);
  r.qwidth=robust;
  double lo=q005, hi=q995;
  if(!(hi>lo)) {lo=*std::min_element(v.begin(),v.end());hi=*std::max_element(v.begin(),v.end());}
  auto* h=new TH1D(("h_"+name).c_str(),"",180,lo,hi);
  h->SetDirectory(&out); for(double x:v)h->Fill(x);
  TF1 first(("first_"+name).c_str(),"gaus",q50-2*robust,q50+2*robust);
  first.SetParameters(h->GetMaximum(),q50,robust);
  h->Fit(&first,"QNR");
  double m=first.GetParameter(1),s=std::abs(first.GetParameter(2));
  r.rangeLo=m-2*s; r.rangeHi=m+2*s;
  auto* fit=new TF1(("fit_"+name).c_str(),"gaus",r.rangeLo,r.rangeHi);
  fit->SetParameters(first.GetParameter(0),m,s);
  auto fr=h->Fit(fit,"QRSN");
  r.status=int(fr);r.fitMean=fit->GetParameter(1);r.fitMeanErr=fit->GetParError(1);
  r.sigma=std::abs(fit->GetParameter(2));r.sigmaErr=fit->GetParError(2);
  r.chi2=fit->GetChisquare();r.ndf=fit->GetNDF();
  out.cd();h->Write();fit->Write();
  return r;
}

const char* matName(int m){return m==0?"EJ-200":m==1?"EJ-204":"EJ-230";}
}

void build_timing_summary() {
  TFile input("sources/timing_events.root","READ");
  auto* tree=dynamic_cast<TTree*>(input.Get("timing_events"));
  if(!tree||tree->GetEntries()!=210000)throw std::runtime_error("Invalid timing_events ROOT");
  char cell[32]={};int mat,x,event,nl,nr,nt,np;double edep,top,tL[K],tR[K];
  tree->SetBranchAddress("cell_id",cell);tree->SetBranchAddress("material",&mat);tree->SetBranchAddress("x_mm",&x);tree->SetBranchAddress("event_id",&event);
  tree->SetBranchAddress("npe_left",&nl);tree->SetBranchAddress("npe_right",&nr);tree->SetBranchAddress("npe_top",&nt);tree->SetBranchAddress("produced_scint",&np);tree->SetBranchAddress("edep_MeV",&edep);tree->SetBranchAddress("t_left_ns",tL);tree->SetBranchAddress("t_right_ns",tR);tree->SetBranchAddress("t_top_first_ns",&top);
  std::map<std::string,std::vector<Event>> data;std::map<std::string,std::pair<int,int>> config;
  for(Long64_t i=0;i<tree->GetEntries();++i){tree->GetEntry(i);Event e{nl,nr,nt,np,edep,top,{},{}};for(int k=0;k<K;k++){e.left[k]=tL[k];e.right[k]=tR[k];}data[cell].push_back(e);config[cell]={mat,x};}

  TFile out("sources/timing_summary.root","RECREATE");
  std::vector<Metric> metrics;
  for(auto& item:data){
    Metric m;m.id=item.first;m.material=config[item.first].first;m.x=config[item.first].second;m.n=item.second.size();
    std::vector<double> vl,vr,vend,vt,vall,vedep,vprod,veff,vtimeL,vtimeR,vdt,vt0,vtop;
    for(const auto&e:item.second){
      vl.push_back(e.npeL);vr.push_back(e.npeR);vend.push_back(0.5*(e.npeL+e.npeR));vt.push_back(e.npeT);vall.push_back(e.npeL+e.npeR+e.npeT);vedep.push_back(e.edep);vprod.push_back(e.produced);veff.push_back(double(e.npeL+e.npeR+e.npeT)/e.produced);vtimeL.push_back(e.left[0]);vtimeR.push_back(e.right[0]);vdt.push_back(e.right[0]-e.left[0]);vt0.push_back(0.5*(e.left[0]+e.right[0]));vtop.push_back(e.top);
    }
    std::tie(m.npeL,m.npeLse)=meanSem(vl);std::tie(m.npeR,m.npeRse)=meanSem(vr);std::tie(m.npeEnd,m.npeEndSe)=meanSem(vend);std::tie(m.npeT,m.npeTse)=meanSem(vt);std::tie(m.npeAll,m.npeAllSe)=meanSem(vall);std::tie(m.edep,m.edepSe)=meanSem(vedep);std::tie(m.produced,m.producedSe)=meanSem(vprod);std::tie(m.eff,m.effSe)=meanSem(veff);std::tie(m.meanL,m.meanLse)=meanSem(vtimeL);std::tie(m.meanR,m.meanRse)=meanSem(vtimeR);std::tie(m.meanDt,m.meanDtSe)=meanSem(vdt);
    m.t0=fitCentral(vt0,m.id+"_t0",out);m.dt=fitCentral(vdt,m.id+"_dt",out);m.top=fitCentral(vtop,m.id+"_top",out);metrics.push_back(m);
  }
  std::sort(metrics.begin(),metrics.end(),[](const Metric&a,const Metric&b){return a.material!=b.material?a.material<b.material:a.x<b.x;});

  std::array<MaterialFit,3> mf;
  for(int material=0;material<3;++material){
    std::vector<double> vx,vy,vey;for(auto&m:metrics)if(m.material==material){vx.push_back(m.x);vy.push_back(m.meanDt);vey.push_back(m.meanDtSe);}
    auto*g=new TGraphErrors(vx.size(),vx.data(),vy.data(),nullptr,vey.data());g->SetName((std::string("g_dt_")+matName(material)).c_str());
    auto*f=new TF1((std::string("fit_dt_")+matName(material)).c_str(),"pol1",-650,650);auto fr=g->Fit(f,"QRSN");
    auto&r=mf[material];r.status=int(fr);r.intercept=f->GetParameter(0);r.interceptErr=f->GetParError(0);r.slope=f->GetParameter(1);r.slopeErr=f->GetParError(1);r.chi2=f->GetChisquare();r.ndf=f->GetNDF();r.veff=-2.0/r.slope;r.veffErr=2.0*r.slopeErr/(r.slope*r.slope);r.cn=cMmNs/refractiveIndex;double ratio=r.veff/r.cn;r.betaDeg=std::acos(ratio)*180/TMath::Pi();r.betaErrDeg=(r.veffErr/r.cn)/std::sqrt(1-ratio*ratio)*180/TMath::Pi();out.cd();g->Write();f->Write();
  }
  for(auto&m:metrics){auto&r=mf[m.material];m.sigmaX=0.5*r.veff*m.dt.sigma;m.sigmaXErr=0.5*std::hypot(r.veff*m.dt.sigmaErr,m.dt.sigma*r.veffErr);}

  std::ofstream csv("sources/timing_summary.csv");csv<<std::setprecision(12)<<"cell_id,material,x_mm,N,npe_left,npe_left_sem,npe_right,npe_right_sem,npe_end,npe_end_sem,npe_top,npe_top_sem,npe_total,npe_total_sem,edep_MeV,edep_sem,produced_scint,produced_sem,detection_fraction,detection_fraction_sem,mean_tL_ns,mean_tL_sem,mean_tR_ns,mean_tR_sem,mean_dt_ns,mean_dt_sem,t0_sigma_ns,t0_sigma_err,t0_rms_ns,t0_rms_err,t0_qwidth_ns,t0_fit_mean_ns,t0_fit_mean_err,t0_fit_lo_ns,t0_fit_hi_ns,t0_chi2,t0_ndf,t0_fit_status,dt_sigma_ns,dt_sigma_err,dt_rms_ns,dt_rms_err,dt_qwidth_ns,dt_chi2,dt_ndf,top_sigma_ns,top_sigma_err,top_rms_ns,top_rms_err,top_qwidth_ns,top_chi2,top_ndf,sigma_x_mm,sigma_x_err_mm\n";
  for(auto&m:metrics)csv<<m.id<<','<<matName(m.material)<<','<<m.x<<','<<m.n<<','<<m.npeL<<','<<m.npeLse<<','<<m.npeR<<','<<m.npeRse<<','<<m.npeEnd<<','<<m.npeEndSe<<','<<m.npeT<<','<<m.npeTse<<','<<m.npeAll<<','<<m.npeAllSe<<','<<m.edep<<','<<m.edepSe<<','<<m.produced<<','<<m.producedSe<<','<<m.eff<<','<<m.effSe<<','<<m.meanL<<','<<m.meanLse<<','<<m.meanR<<','<<m.meanRse<<','<<m.meanDt<<','<<m.meanDtSe<<','<<m.t0.sigma<<','<<m.t0.sigmaErr<<','<<m.t0.globalRms<<','<<m.t0.globalRmsErr<<','<<m.t0.qwidth<<','<<m.t0.fitMean<<','<<m.t0.fitMeanErr<<','<<m.t0.rangeLo<<','<<m.t0.rangeHi<<','<<m.t0.chi2<<','<<m.t0.ndf<<','<<m.t0.status<<','<<m.dt.sigma<<','<<m.dt.sigmaErr<<','<<m.dt.globalRms<<','<<m.dt.globalRmsErr<<','<<m.dt.qwidth<<','<<m.dt.chi2<<','<<m.dt.ndf<<','<<m.top.sigma<<','<<m.top.sigmaErr<<','<<m.top.globalRms<<','<<m.top.globalRmsErr<<','<<m.top.qwidth<<','<<m.top.chi2<<','<<m.top.ndf<<','<<m.sigmaX<<','<<m.sigmaXErr<<'\n';

  std::ofstream mcsv("sources/material_summary.csv");mcsv<<std::setprecision(12)<<"material,decay_ns,attenuation_m,light_yield_ph_per_MeV,dt_intercept_ns,dt_intercept_err,dt_slope_ns_per_mm,dt_slope_err,fit_chi2,fit_ndf,fit_status,v_eff_mm_per_ns,v_eff_err,c_over_n_mm_per_ns,beta_eff_deg,beta_eff_err_deg\n";const double decay[3]={2.1,1.8,1.5},att[3]={3.8,1.6,1.2},yield[3]={10000,10400,9700};for(int m=0;m<3;m++){auto&r=mf[m];mcsv<<matName(m)<<','<<decay[m]<<','<<att[m]<<','<<yield[m]<<','<<r.intercept<<','<<r.interceptErr<<','<<r.slope<<','<<r.slopeErr<<','<<r.chi2<<','<<r.ndf<<','<<r.status<<','<<r.veff<<','<<r.veffErr<<','<<r.cn<<','<<r.betaDeg<<','<<r.betaErrDeg<<'\n';}

  auto& center=data.at("EJ230_xp0");std::ofstream kcsv("sources/order_scan_ej230_x0.csv");kcsv<<std::setprecision(12)<<"k,N,sigma_t0_ns,sigma_t0_err,rms_t0_ns,rms_t0_err,fit_mean_ns,fit_mean_err,fit_lo_ns,fit_hi_ns,chi2,ndf,fit_status,mean_first_m_sigma_ns,mean_first_m_sigma_err,mean_first_m_rms_ns,mean_first_m_rms_err,mean_first_m_chi2,mean_first_m_ndf,mean_first_m_fit_status\n";for(int k=0;k<K;k++){std::vector<double> values,averages;for(auto&e:center){values.push_back(0.5*(e.left[k]+e.right[k]));double sl=0,sr=0;for(int j=0;j<=k;j++){sl+=e.left[j];sr+=e.right[j];}averages.push_back(0.5*(sl+sr)/(k+1));}auto r=fitCentral(values,"EJ230_xp0_t0_k"+std::to_string(k+1),out);auto a=fitCentral(averages,"EJ230_xp0_t0_mean_first_"+std::to_string(k+1),out);kcsv<<k+1<<','<<values.size()<<','<<r.sigma<<','<<r.sigmaErr<<','<<r.globalRms<<','<<r.globalRmsErr<<','<<r.fitMean<<','<<r.fitMeanErr<<','<<r.rangeLo<<','<<r.rangeHi<<','<<r.chi2<<','<<r.ndf<<','<<r.status<<','<<a.sigma<<','<<a.sigmaErr<<','<<a.globalRms<<','<<a.globalRmsErr<<','<<a.chi2<<','<<a.ndf<<','<<a.status<<'\n';}

  out.cd();TNamed source("source","sources/timing_events.root; all timing widths and fits produced with ROOT");TNamed method("method","tL/tR are kth detected PE over 8 SiPMs per END; baseline k=1; scan also includes the arithmetic mean of the first m PE on each END; central Gaussian: initial median +/- 2*(q84-q16)/2 then refit mean +/- 2 sigma; global RMS retained");source.Write();method.Write();out.Close();
  std::cout<<"Wrote "<<metrics.size()<<" cell summaries and 20 order-statistic fits"<<std::endl;
}
