#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TMatrixDSym.h>
#include <TNamed.h>
#include <TRandom3.h>
#include <TTree.h>
#include <Math/MinimizerOptions.h>

#include <RooArgSet.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooGaussian.h>
#include <RooRealVar.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

using namespace RooFit;

namespace {
constexpr int kOrder=20;
constexpr int kBootstrap=500;
constexpr UInt_t kBootstrapSeed=26091443;
constexpr double kHistHalfRangeNs=0.32;

struct Event { std::array<double,kOrder> l,r; };
struct Cell { std::string id,material;int x=0;std::vector<Event> events; };
struct Width { double mean=0,median=0,rms=0,q=0,rmsErr=0,qErr=0; };
struct Fit { double amp=0,ampErr=0,mu0=0,sigma0=0,mu=0,muErr=0,sigma=0,sigmaErr=0,chi2=0,edm=0,minFcn=0,lo=0,hi=0;int ndf=0,status=-1,covStatus=-1;TMatrixDSym cov{3}; };

double quantile(const std::vector<double>& sorted,double p){double z=p*(sorted.size()-1),a=std::floor(z),u=z-a;size_t n=a;return n+1<sorted.size()?sorted[n]*(1-u)+sorted[n+1]*u:sorted[n];}
std::pair<double,double> rmsQ(std::vector<double> v){double m=std::accumulate(v.begin(),v.end(),0.0)/v.size(),ss=0;for(double x:v)ss+=(x-m)*(x-m);std::sort(v.begin(),v.end());return {std::sqrt(ss/v.size()),.5*(quantile(v,.84)-quantile(v,.16))};}
double sampleSd(const std::vector<double>&v){double m=std::accumulate(v.begin(),v.end(),0.0)/v.size(),ss=0;for(double x:v)ss+=(x-m)*(x-m);return std::sqrt(ss/(v.size()-1));}
Width measure(const std::vector<double>& values,UInt_t seed){Width w;w.mean=std::accumulate(values.begin(),values.end(),0.0)/values.size();auto s=values;std::sort(s.begin(),s.end());w.median=quantile(s,.5);auto rq=rmsQ(values);w.rms=rq.first;w.q=rq.second;TRandom3 rng(seed);std::vector<double> br,bq,sample(values.size());br.reserve(kBootstrap);bq.reserve(kBootstrap);for(int b=0;b<kBootstrap;b++){for(size_t j=0;j<values.size();j++)sample[j]=values[rng.Integer(values.size())];auto x=rmsQ(sample);br.push_back(x.first);bq.push_back(x.second);}w.rmsErr=sampleSd(br);w.qErr=sampleSd(bq);return w;}

TH1D* makeHist(const std::string&name,const std::vector<double>&v,double center,double binWidthNs,double offsetFraction){int nb=int(std::lround(2*kHistHalfRangeNs/binWidthNs));double lo=center-kHistHalfRangeNs+offsetFraction*binWidthNs,hi=lo+nb*binWidthNs;auto*h=new TH1D(name.c_str(),"",nb,lo,hi);h->Sumw2();for(double x:v)h->Fill(x);return h;}

Fit fitHist(TH1D& h,double med,double q,double windowScale,char initKind,bool likelihood,double forcedLo=999,double forcedHi=999){Fit r;r.lo=forcedLo<900?forcedLo:med-windowScale*q;r.hi=forcedHi<900?forcedHi:med+windowScale*q;double mode=h.GetBinCenter(h.GetMaximumBin());if(initKind=='A'){r.mu0=med;r.sigma0=q;}else if(initKind=='B'){r.mu0=h.GetMean();r.sigma0=h.GetRMS();}else{r.mu0=mode;r.sigma0=q;}TF1 f("tmp_gaus","gaus",r.lo,r.hi);f.SetParameters(h.GetMaximum(),r.mu0,r.sigma0);std::string opt=likelihood?"QLRS0":"QRS0";TFitResultPtr fr=h.Fit(&f,opt.c_str());r.status=int(fr);r.covStatus=fr.Get()?fr->CovMatrixStatus():-1;r.amp=f.GetParameter(0);r.ampErr=f.GetParError(0);r.mu=f.GetParameter(1);r.muErr=f.GetParError(1);r.sigma=std::abs(f.GetParameter(2));r.sigmaErr=f.GetParError(2);r.chi2=f.GetChisquare();r.ndf=f.GetNDF();r.edm=fr.Get()?fr->Edm():-1;r.minFcn=fr.Get()?fr->MinFcnValue():-1;if(fr.Get())r.cov=fr->GetCovarianceMatrix();return r;}

Fit fitUnbinned(const std::vector<double>&v,double med,double q){Fit r;r.lo=med-2*q;r.hi=med+2*q;r.mu0=med;r.sigma0=q;RooRealVar x("t0","t0",r.lo,r.hi,"ns");RooDataSet data("data","data",RooArgSet(x));for(double z:v)if(z>=r.lo&&z<=r.hi){x.setVal(z);data.add(RooArgSet(x));}RooRealVar mu("mu","mu",r.mu0,r.lo,r.hi);RooRealVar sigma("sigma","sigma",r.sigma0,.01,.25);RooGaussian pdf("gaussian","gaussian",x,mu,sigma);auto fr=std::unique_ptr<RooFitResult>(pdf.fitTo(data,Save(),PrintLevel(-1),Warnings(false),Verbose(false),Strategy(1)));r.status=fr->status();r.covStatus=fr->covQual();r.mu=mu.getVal();r.muErr=mu.getError();r.sigma=sigma.getVal();r.sigmaErr=sigma.getError();r.edm=fr->edm();r.minFcn=fr->minNll();r.chi2=-1;r.ndf=-1;r.cov.ResizeTo(2,2);r.cov=fr->covarianceMatrix();return r;}

std::vector<double> estimator(const Cell&c,bool mean5){std::vector<double>v;v.reserve(c.events.size());for(auto&e:c.events){if(!mean5)v.push_back(.5*(e.l[0]+e.r[0]));else{double sl=0,sr=0;for(int k=0;k<5;k++){sl+=e.l[k];sr+=e.r[k];}v.push_back(.1*(sl+sr));}}return v;}
std::string cellId(const std::string&m,int x){std::string p=m=="EJ-200"?"EJ200_":m=="EJ-204"?"EJ204_":"EJ230_";return p+(x<0?"xm"+std::to_string(-x):"xp"+std::to_string(x));}
}

void timing_fit_stability(){
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Migrad");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.01);
  TFile input("../../presentations/v9/sources/timing_events.root","READ");auto*tree=(TTree*)input.Get("timing_events");if(!tree||tree->GetEntries()!=210000)throw std::runtime_error("Invalid timing event source");
  char cid[32]={};int material=0,x=0;double l[kOrder],r[kOrder];tree->SetBranchAddress("cell_id",cid);tree->SetBranchAddress("material",&material);tree->SetBranchAddress("x_mm",&x);tree->SetBranchAddress("t_left_ns",l);tree->SetBranchAddress("t_right_ns",r);std::map<std::string,Cell> map;
  for(Long64_t j=0;j<tree->GetEntries();j++){tree->GetEntry(j);auto&c=map[cid];c.id=cid;c.material=material==0?"EJ-200":material==1?"EJ-204":"EJ-230";c.x=x;Event e;for(int k=0;k<kOrder;k++){e.l[k]=l[k];e.r[k]=r[k];}c.events.push_back(e);}if(map.size()!=21)throw std::runtime_error("Expected 21 cells");

  TFile out("sources/timing_fit_stability.root","RECREATE");
  std::ofstream csv("sources/timing_fit_stability.csv");csv<<std::setprecision(12)<<"cell_id,material,x_mm,estimator,study,bin_width_ps,bin_offset_fraction,window_scale,initialization,method,common_mirror_window,amplitude,amplitude_err,mu0_ns,sigma0_ns,fit_lo_ns,fit_hi_ns,fit_mu_ns,fit_mu_err_ns,fit_sigma_ns,fit_sigma_err_ns,chi2,ndf,chi2_ndf,fit_status,covariance_status,edm,min_fcn\n";
  std::ofstream widths("sources/timing_widths.csv");widths<<std::setprecision(12)<<"cell_id,material,x_mm,estimator,N,mean_ns,median_ns,RMS_ns,RMS_bootstrap_err_ns,sigma68_ns,sigma68_bootstrap_err_ns,baseline_sigma_G_ns,baseline_sigma_G_err_ns,baseline_chi2,baseline_ndf,baseline_fit_status,baseline_covariance_status\n";
  std::map<std::string,std::map<std::string,Width>> wm;std::map<std::string,std::map<std::string,Fit>> fm;UInt_t bseed=kBootstrapSeed;
  auto emit=[&](const Cell&c,const std::string&est,const std::string&study,double bw,double off,double ws,char init,const std::string&method,bool common,const Fit&f){csv<<c.id<<','<<c.material<<','<<c.x<<','<<est<<','<<study<<','<<1000*bw<<','<<off<<','<<ws<<','<<init<<','<<method<<','<<(common?1:0)<<','<<f.amp<<','<<f.ampErr<<','<<f.mu0<<','<<f.sigma0<<','<<f.lo<<','<<f.hi<<','<<f.mu<<','<<f.muErr<<','<<f.sigma<<','<<f.sigmaErr<<','<<f.chi2<<','<<f.ndf<<','<<(f.ndf>0?f.chi2/f.ndf:-1)<<','<<f.status<<','<<f.covStatus<<','<<f.edm<<','<<f.minFcn<<'\n';};

  for(auto&it:map){auto&c=it.second;for(bool mean5:{false,true}){std::string est=mean5?"mean_first_5":"first_pe";auto v=estimator(c,mean5);auto w=measure(v,bseed++);wm[c.id][est]=w;auto*h=makeHist("baseline_"+c.id+"_"+est,v,w.median,.004,0);auto f=fitHist(*h,w.median,w.q,2.0,'A',false);fm[c.id][est]=f;widths<<c.id<<','<<c.material<<','<<c.x<<','<<est<<','<<v.size()<<','<<w.mean<<','<<w.median<<','<<w.rms<<','<<w.rmsErr<<','<<w.q<<','<<w.qErr<<','<<f.sigma<<','<<f.sigmaErr<<','<<f.chi2<<','<<f.ndf<<','<<f.status<<','<<f.covStatus<<'\n';if(!mean5){out.cd();auto*dir=out.mkdir(c.id.c_str());dir->cd();h->SetDirectory(dir);h->Write("t0_histogram");double fullLo=.004*std::floor(*std::min_element(v.begin(),v.end())/.004),fullHi=.004*std::ceil(*std::max_element(v.begin(),v.end())/.004);int fullBins=std::max(1,int(std::lround((fullHi-fullLo)/.004)));auto*full=new TH1D(("complete_"+c.id).c_str(),"",fullBins,fullLo,fullHi);for(double z:v)full->Fill(z);full->SetDirectory(dir);full->Write("t0_complete_histogram");TF1 gf("gaussian_fit","gaus",f.lo,f.hi);gf.SetParameters(f.amp,f.mu,f.sigma);gf.SetParError(0,f.ampErr);gf.SetParError(1,f.muErr);gf.SetParError(2,f.sigmaErr);gf.SetChisquare(f.chi2);gf.SetNDF(f.ndf);gf.Write();f.cov.Write("fit_covariance");emit(c,est,"baseline",.004,0,2,'A',"chi2",false,f);
        for(double bw:{.002,.004,.008})for(double off:{0.0,.5}){auto*hb=makeHist("hbin",v,w.median,bw,off);auto z=fitHist(*hb,w.median,w.q,2,'A',false);emit(c,est,"binning_offset",bw,off,2,'A',"chi2",false,z);delete hb;}
        auto*hi=makeHist("hinit",v,w.median,.004,0);for(char ik:{'A','B','C'})emit(c,est,"initialization",.004,0,2,ik,"chi2",false,fitHist(*hi,w.median,w.q,2,ik,false));for(double ws:{1.5,2.0,2.5})emit(c,est,"fit_window",.004,0,ws,'A',"chi2",false,fitHist(*hi,w.median,w.q,ws,'A',false));emit(c,est,"fit_method",.004,0,2,'A',"chi2",false,fitHist(*hi,w.median,w.q,2,'A',false));emit(c,est,"fit_method",.004,0,2,'A',"binned_likelihood",false,fitHist(*hi,w.median,w.q,2,'A',true));auto uf=fitUnbinned(v,w.median,w.q);emit(c,est,"fit_method",0,0,2,'A',"RooFit_unbinned",false,uf);delete hi;
      }delete h;}}

  // Mirrored common-window fits use centered samples and one pooled robust scale.
  for(auto mat:{"EJ-200","EJ-204","EJ-230"})for(int a:{200,500,650}){auto&cn=map.at(cellId(mat,-a));auto&cp=map.at(cellId(mat,a));auto vn=estimator(cn,false),vp=estimator(cp,false);double mn=wm[cn.id]["first_pe"].median,mp=wm[cp.id]["first_pe"].median;std::vector<double>pool;pool.reserve(vn.size()+vp.size());for(double&z:vn){z-=mn;pool.push_back(z);}for(double&z:vp){z-=mp;pool.push_back(z);}std::sort(pool.begin(),pool.end());double qc=.5*(quantile(pool,.84)-quantile(pool,.16));for(auto item:{std::make_pair(&cn,&vn),std::make_pair(&cp,&vp)}){auto*h=makeHist("hcommon",*item.second,0,.004,0);auto f=fitHist(*h,0,qc,2,'A',false,-2*qc,2*qc);emit(*item.first,"first_pe","common_mirror_window",.004,0,2,'A',"chi2",true,f);delete h;}}

  out.cd();TNamed method("method","Minuit2/Migrad tolerance 0.01; physical bin widths 2/4/8 ps; offsets 0/half-bin; robust-window scales 1.5/2/2.5; initializations A median/q68, B histogram mean/RMS, C mode/q68; chi2, binned likelihood and unbinned RooFit; 500 bootstraps");method.Write();TNamed source("source","presentations/v9/sources/timing_events.root; current corrected transport; 210000 events");source.Write();out.Close();std::cout<<"Wrote stability diagnostics for "<<map.size()<<" cells; RooFit available and executed"<<std::endl;
}
