#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TMatrixDSym.h>
#include <TNamed.h>
#include <TProfile.h>
#include <TTree.h>
#include <Math/MinimizerOptions.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {
constexpr int kOrder = 20;
constexpr double kHalfLengthM = 0.7;
struct E { int id=0,nl=0,nr=0; double tl=0,tr=0,ts=0,t0=0,dt=0,m5l=0,m5r=0,t0m5=0; };
struct Cell { std::string id,mat; int x=0; std::vector<E> e; };
struct Stat { double mean=0,sem=0,median=0,rms=0; };
struct Fit { std::string model; int status=-1,cov=-1,npar=0,ndf=0; double chi2=0,prob=0,aic=0; std::vector<double>p,pe; TMatrixDSym cm; Fit(int n=1):npar(n),p(n),pe(n),cm(n){} };

Stat stat(std::vector<double> v){Stat s;s.mean=std::accumulate(v.begin(),v.end(),0.)/v.size();double ss=0;for(double x:v)ss+=(x-s.mean)*(x-s.mean);s.rms=std::sqrt(ss/v.size());s.sem=s.rms/std::sqrt(v.size());std::sort(v.begin(),v.end());s.median=.5*(v[(v.size()-1)/2]+v[v.size()/2]);return s;}
double pearson(const std::vector<double>&a,const std::vector<double>&b){double ma=std::accumulate(a.begin(),a.end(),0.)/a.size(),mb=std::accumulate(b.begin(),b.end(),0.)/b.size(),sa=0,sb=0,sab=0;for(size_t i=0;i<a.size();i++){double x=a[i]-ma,y=b[i]-mb;sa+=x*x;sb+=y*y;sab+=x*y;}return sab/std::sqrt(sa*sb);}
std::string safe(std::string s){for(char&c:s)if(c=='-')c='_';return s;}
std::string cid(const std::string&m,int x){std::string p=m=="EJ-200"?"EJ200_":m=="EJ-204"?"EJ204_":"EJ230_";return p+(x<0?"xm"+std::to_string(-x):"xp"+std::to_string(x));}
std::vector<double> take(const Cell&c,const std::string&v){std::vector<double>z;z.reserve(c.e.size());for(auto&e:c.e){if(v=="tl")z.push_back(e.tl);else if(v=="tr")z.push_back(e.tr);else if(v=="tsum")z.push_back(e.ts);else if(v=="t0")z.push_back(e.t0);else if(v=="dt")z.push_back(e.dt);else if(v=="m5l")z.push_back(e.m5l);else if(v=="m5r")z.push_back(e.m5r);else if(v=="t0m5")z.push_back(e.t0m5);}return z;}
Fit doFit(TGraphErrors&g,const std::string&model,const std::string&name){int np=model=="constant"?1:model=="linear"||model=="even_quadratic"?2:3;Fit o(np);o.model=model;std::string form=model=="constant"?"[0]":model=="linear"?"[0]+[1]*x":model=="even_quadratic"?"[0]+[1]*x*x":"[0]+[1]*x+[2]*x*x";TF1 f(name.c_str(),form.c_str(),-.66,1.36);double my=std::accumulate(g.GetY(),g.GetY()+g.GetN(),0.)/g.GetN();f.SetParameter(0,my);if(np>1)f.SetParameter(1,model=="linear"?5.5:0);if(np>2){f.SetParameter(1,5.5);f.SetParameter(2,0);}TFitResultPtr r=g.Fit(&f,"QRS0");o.status=int(r);o.cov=r.Get()?r->CovMatrixStatus():-1;o.chi2=f.GetChisquare();o.ndf=f.GetNDF();o.prob=f.GetProb();o.aic=o.chi2+2*np;for(int j=0;j<np;j++){o.p[j]=f.GetParameter(j);o.pe[j]=f.GetParError(j);}if(r.Get())o.cm=r->GetCovarianceMatrix();f.Write(name.c_str());o.cm.Write((name+"_covariance").c_str());return o;}
}

void analyze_tsum_veff(){
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Migrad"); ROOT::Math::MinimizerOptions::SetDefaultTolerance(.01);
  TFile in("../../presentations/v9/sources/timing_events.root","READ");auto*t=(TTree*)in.Get("timing_events");if(!t||t->GetEntries()!=210000)throw std::runtime_error("bad timing event input");
  char id[32]={};int im=0,x=0,event=0,nl=0,nr=0;double l[kOrder],r[kOrder];t->SetBranchAddress("cell_id",id);t->SetBranchAddress("material",&im);t->SetBranchAddress("x_mm",&x);t->SetBranchAddress("event_id",&event);t->SetBranchAddress("npe_left",&nl);t->SetBranchAddress("npe_right",&nr);t->SetBranchAddress("t_left_ns",l);t->SetBranchAddress("t_right_ns",r);
  std::map<std::string,Cell> cells;
  TFile out("sources/tsum_veff.root","RECREATE");TTree ev("derived_events","Event-level timing and light observables");char oid[32]={},omat[16]={};int ox=0,oe=0,onl=0,onr=0,onend=0;double otl=0,otr=0,ots=0,ot0=0,odt=0,om5l=0,om5r=0,ot0m5=0;
  ev.Branch("cell_id",oid,"cell_id/C");ev.Branch("material",omat,"material/C");ev.Branch("x_mm",&ox);ev.Branch("event_id",&oe);ev.Branch("npe_left",&onl);ev.Branch("npe_right",&onr);ev.Branch("npe_end",&onend);ev.Branch("t_left_ns",&otl);ev.Branch("t_right_ns",&otr);ev.Branch("tsum_ns",&ots);ev.Branch("t0_ns",&ot0);ev.Branch("delta_t_ns",&odt);ev.Branch("mean5_left_ns",&om5l);ev.Branch("mean5_right_ns",&om5r);ev.Branch("t0_mean5_ns",&ot0m5);
  for(Long64_t j=0;j<t->GetEntries();j++){t->GetEntry(j);auto&c=cells[id];c.id=id;c.mat=im==0?"EJ-200":im==1?"EJ-204":"EJ-230";c.x=x;E e;e.id=event;e.nl=nl;e.nr=nr;e.tl=l[0];e.tr=r[0];e.ts=e.tl+e.tr;e.t0=.5*e.ts;e.dt=e.tr-e.tl;for(int k=0;k<5;k++){e.m5l+=l[k]/5.;e.m5r+=r[k]/5.;}e.t0m5=.5*(e.m5l+e.m5r);c.e.push_back(e);snprintf(oid,sizeof(oid),"%s",c.id.c_str());snprintf(omat,sizeof(omat),"%s",c.mat.c_str());ox=x;oe=event;onl=nl;onr=nr;onend=nl+nr;otl=e.tl;otr=e.tr;ots=e.ts;ot0=e.t0;odt=e.dt;om5l=e.m5l;om5r=e.m5r;ot0m5=e.t0m5;ev.Fill();}
  if(cells.size()!=21)throw std::runtime_error("expected 21 cells");out.cd();ev.Write();

  std::ofstream means("sources/means.csv");means<<std::setprecision(12)<<"cell_id,material,x_mm,N,mean_tL_ns,sem_tL_ns,mean_tR_ns,sem_tR_ns,mean_Tsum_ns,sem_Tsum_ns,mean_T0_ns,sem_T0_ns,median_T0_ns,RMS_T0_ns,mean_DeltaT_ns,sem_DeltaT_ns,mean_tL_m5_ns,sem_tL_m5_ns,mean_tR_m5_ns,sem_tR_m5_ns,mean_T0_m5_ns,sem_T0_m5_ns,median_T0_m5_ns,RMS_T0_m5_ns,mean_Npe_L,mean_Npe_R,mean_Npe_END\n";
  std::ofstream corr("sources/light_correlations.csv");corr<<std::setprecision(12)<<"cell_id,material,x_mm,N,pearson_T0_NpeEND,pearson_tL_NpeL,pearson_tR_NpeR,pearson_T0m5_NpeEND\n";
  std::map<std::string,std::map<int,std::map<std::string,Stat>>> sm;
  auto*hdir=out.mkdir("distributions");auto*pdir=out.mkdir("light_profiles");
  for(auto&kv:cells){auto&c=kv.second;for(auto v:{"tl","tr","tsum","t0","dt","m5l","m5r","t0m5"})sm[c.mat][c.x][v]=stat(take(c,v));auto&z=sm[c.mat][c.x];double nml=0,nmr=0;std::vector<double>vt0,vm5,vnl,vnr,vne,vtl,vtr;for(auto&e:c.e){nml+=e.nl;nmr+=e.nr;vt0.push_back(e.t0);vm5.push_back(e.t0m5);vnl.push_back(e.nl);vnr.push_back(e.nr);vne.push_back(e.nl+e.nr);vtl.push_back(e.tl);vtr.push_back(e.tr);}nml/=c.e.size();nmr/=c.e.size();means<<c.id<<','<<c.mat<<','<<c.x<<','<<c.e.size()<<','<<z["tl"].mean<<','<<z["tl"].sem<<','<<z["tr"].mean<<','<<z["tr"].sem<<','<<z["tsum"].mean<<','<<z["tsum"].sem<<','<<z["t0"].mean<<','<<z["t0"].sem<<','<<z["t0"].median<<','<<z["t0"].rms<<','<<z["dt"].mean<<','<<z["dt"].sem<<','<<z["m5l"].mean<<','<<z["m5l"].sem<<','<<z["m5r"].mean<<','<<z["m5r"].sem<<','<<z["t0m5"].mean<<','<<z["t0m5"].sem<<','<<z["t0m5"].median<<','<<z["t0m5"].rms<<','<<nml<<','<<nmr<<','<<nml+nmr<<'\n';corr<<c.id<<','<<c.mat<<','<<c.x<<','<<c.e.size()<<','<<pearson(vt0,vne)<<','<<pearson(vtl,vnl)<<','<<pearson(vtr,vnr)<<','<<pearson(vm5,vne)<<'\n';
    hdir->cd();double lo=*std::min_element(vt0.begin(),vt0.end()),hi=*std::max_element(vt0.begin(),vt0.end());auto*h=new TH1D(("T0_"+c.id).c_str(),"",std::max(50,int(std::ceil((hi-lo)/.004))),lo-.002,hi+.002);for(double q:vt0)h->Fill(q);h->Write();
    pdir->cd();auto makeP=[&](const std::string&name,const std::vector<double>&xx,const std::vector<double>&yy){double a=*std::min_element(xx.begin(),xx.end()),b=*std::max_element(xx.begin(),xx.end());auto*p=new TProfile(name.c_str(),"",40,a-.5,b+.5);for(size_t q=0;q<xx.size();q++)p->Fill(xx[q],yy[q]);p->Write();};makeP("T0_vs_NpeEND_"+c.id,vne,vt0);makeP("tL_vs_NpeL_"+c.id,vnl,vtl);makeP("tR_vs_NpeR_"+c.id,vnr,vtr);makeP("T0m5_vs_NpeEND_"+c.id,vne,vm5);
  }

  std::ofstream fits("sources/model_fits.csv");fits<<std::setprecision(12)<<"observable,material,model,npoints,status,covariance_status,chi2,ndf,chi2_ndf,probability,AIC,p0,p0_err,p1,p1_err,p2,p2_err\n";
  auto*gdir=out.mkdir("position_graphs");auto*fdir=out.mkdir("position_fits");
  for(auto mat:{std::string("EJ-200"),std::string("EJ-204"),std::string("EJ-230")})for(auto obs:{std::string("t0"),std::string("tsum")}){std::vector<double>xx,yy,xe,ye;for(int xi:{-650,-500,-200,0,200,500,650}){xx.push_back(xi/1000.);xe.push_back(0);yy.push_back(sm[mat][xi][obs].mean);ye.push_back(sm[mat][xi][obs].sem);}gdir->cd();TGraphErrors g(xx.size(),xx.data(),yy.data(),xe.data(),ye.data());std::string gn=obs+"_vs_x_"+safe(mat);g.SetName(gn.c_str());g.Write();fdir->cd();for(auto mod:{std::string("constant"),std::string("linear"),std::string("even_quadratic"),std::string("full_quadratic")}){auto ft=doFit(g,mod,gn+"_"+mod);fits<<obs<<','<<mat<<','<<mod<<','<<g.GetN()<<','<<ft.status<<','<<ft.cov<<','<<ft.chi2<<','<<ft.ndf<<','<<(ft.ndf?ft.chi2/ft.ndf:-1)<<','<<ft.prob<<','<<ft.aic;for(int k=0;k<3;k++)fits<<','<<(k<ft.npar?ft.p[k]:0)<<','<<(k<ft.npar?ft.pe[k]:0);fits<<'\n';}}
  // EJ-230 mean-first-five position fits.
  {std::string mat="EJ-230",obs="t0m5";std::vector<double>xx,yy,xe,ye;for(int xi:{-650,-500,-200,0,200,500,650}){xx.push_back(xi/1000.);xe.push_back(0);yy.push_back(sm[mat][xi][obs].mean);ye.push_back(sm[mat][xi][obs].sem);}gdir->cd();TGraphErrors g(xx.size(),xx.data(),yy.data(),xe.data(),ye.data());std::string gn="t0m5_vs_x_EJ_230";g.SetName(gn.c_str());g.Write();fdir->cd();for(auto mod:{std::string("constant"),std::string("linear"),std::string("even_quadratic"),std::string("full_quadratic")}){auto ft=doFit(g,mod,gn+"_"+mod);fits<<obs<<','<<mat<<','<<mod<<','<<g.GetN()<<','<<ft.status<<','<<ft.cov<<','<<ft.chi2<<','<<ft.ndf<<','<<(ft.ndf?ft.chi2/ft.ndf:-1)<<','<<ft.prob<<','<<ft.aic;for(int k=0;k<3;k++)fits<<','<<(k<ft.npar?ft.p[k]:0)<<','<<(k<ft.npar?ft.pe[k]:0);fits<<'\n';}}

  std::ofstream mirror("sources/mirror_tests.csv");mirror<<std::setprecision(12)<<"material,estimator,abs_x_mm,delta_T0_ns,error_ns,z,delta_Tsum_ns,error_Tsum_ns,z_Tsum\n";for(auto mat:{std::string("EJ-200"),std::string("EJ-204"),std::string("EJ-230")})for(auto est:{std::string("first_pe"),std::string("mean_first_5")}){if(est=="mean_first_5"&&mat!="EJ-230")continue;std::string o=est=="first_pe"?"t0":"t0m5";for(int a:{200,500,650}){auto&mn=sm[mat][-a][o];auto&mp=sm[mat][a][o];double d=mp.mean-mn.mean,de=std::hypot(mp.sem,mn.sem);auto&sn=sm[mat][-a]["tsum"];auto&sp=sm[mat][a]["tsum"];double ds=sp.mean-sn.mean,dse=std::hypot(sp.sem,sn.sem);mirror<<mat<<','<<est<<','<<a<<','<<d<<','<<de<<','<<d/de<<','<<(est=="first_pe"?ds:2*d)<<','<<(est=="first_pe"?dse:2*de)<<','<<d/de<<'\n';}}

  // At equal distance d, LEFT(x) must agree with RIGHT(-x) under detector symmetry.
  std::ofstream collapse("sources/propagation_collapse.csv");collapse<<std::setprecision(12)<<"material,estimator,distance_m,left_x_mm,right_x_mm,mean_left_ns,sem_left_ns,mean_right_ns,sem_right_ns,left_minus_right_ns,error_ns,z\n";
  for(auto mat:{std::string("EJ-200"),std::string("EJ-204"),std::string("EJ-230")})for(auto est:{std::string("first_pe"),std::string("mean_first_5")}){if(est=="mean_first_5"&&mat!="EJ-230")continue;std::string vl=est=="first_pe"?"tl":"m5l",vr=est=="first_pe"?"tr":"m5r";for(int xi:{-650,-500,-200,0,200,500,650}){auto&ls=sm[mat][xi][vl];auto&rs=sm[mat][-xi][vr];double delta=ls.mean-rs.mean,err=std::hypot(ls.sem,rs.sem);collapse<<mat<<','<<est<<','<<kHalfLengthM+xi/1000.<<','<<xi<<','<<-xi<<','<<ls.mean<<','<<ls.sem<<','<<rs.mean<<','<<rs.sem<<','<<delta<<','<<err<<','<<delta/err<<'\n';}}

  std::ofstream prop("sources/propagation_fits.csv");prop<<std::setprecision(12)<<"material,estimator,model,npoints,status,covariance_status,chi2,ndf,chi2_ndf,probability,AIC,b0_ns,b0_err_ns,b1_ns_per_m,b1_err_ns_per_m,b2_ns_per_m2,b2_err_ns_per_m2,delta_chi2_vs_linear,veff_at_0p05_mm_per_ns,veff_at_0p7_mm_per_ns,veff_at_1p35_mm_per_ns\n";
  auto*pgdir=out.mkdir("propagation_graphs");auto*pfdir=out.mkdir("propagation_fits");
  for(auto mat:{std::string("EJ-200"),std::string("EJ-204"),std::string("EJ-230")})for(auto est:{std::string("first_pe"),std::string("mean_first_5")}){if(est=="mean_first_5"&&mat!="EJ-230")continue;std::string vl=est=="first_pe"?"tl":"m5l",vr=est=="first_pe"?"tr":"m5r";std::vector<double>d,y,de,ye,side;for(int xi:{-650,-500,-200,0,200,500,650}){d.push_back(kHalfLengthM+xi/1000.);y.push_back(sm[mat][xi][vl].mean);de.push_back(0);ye.push_back(sm[mat][xi][vl].sem);side.push_back(0);d.push_back(kHalfLengthM-xi/1000.);y.push_back(sm[mat][xi][vr].mean);de.push_back(0);ye.push_back(sm[mat][xi][vr].sem);side.push_back(1);}pgdir->cd();TGraphErrors g(d.size(),d.data(),y.data(),de.data(),ye.data());std::string gn="g_"+est+"_"+safe(mat);g.SetName(gn.c_str());g.Write();TTree pts((gn+"_points").c_str(),"LEFT and RIGHT propagation means");double pd=0,py=0,pe=0;int ps=0,px=0;pts.Branch("distance_m",&pd);pts.Branch("mean_time_ns",&py);pts.Branch("sem_ns",&pe);pts.Branch("side",&ps);pts.Branch("x_mm",&px);for(size_t q=0;q<d.size();q++){pd=d[q];py=y[q];pe=ye[q];ps=side[q];px=(q/2==0?-650:q/2==1?-500:q/2==2?-200:q/2==3?0:q/2==4?200:q/2==5?500:650);pts.Fill();}pts.Write();pfdir->cd();auto fl=doFit(g,"linear",gn+"_linear"),fq=doFit(g,"full_quadratic",gn+"_quadratic");double dc=fl.chi2-fq.chi2;for(auto&ft:{fl,fq}){double b1=ft.p[1],b2=ft.npar>2?ft.p[2]:0;auto vv=[&](double dm){return 1000./(b1+2*b2*dm);};prop<<mat<<','<<est<<','<<ft.model<<','<<g.GetN()<<','<<ft.status<<','<<ft.cov<<','<<ft.chi2<<','<<ft.ndf<<','<<(ft.ndf?ft.chi2/ft.ndf:-1)<<','<<ft.prob<<','<<ft.aic<<','<<ft.p[0]<<','<<ft.pe[0]<<','<<ft.p[1]<<','<<ft.pe[1]<<','<<(ft.npar>2?ft.p[2]:0)<<','<<(ft.npar>2?ft.pe[2]:0)<<','<<(ft.model=="full_quadratic"?dc:0)<<','<<vv(.05)<<','<<vv(.7)<<','<<vv(1.35)<<'\n';}}

  out.cd();TNamed source("source","presentations/v9/sources/timing_events.root; 21 current corrected-transport cells; 210000 events; no new simulation");TNamed defs("definitions","first PE over IDs 0..7 LEFT and 8..15 RIGHT; Tsum=tL+tR; T0=Tsum/2; DeltaT=tR-tL; distances dL=0.7+x and dR=0.7-x metres; mean-first-five arithmetic per END");TNamed fitmethod("fit_method","ROOT TGraphErrors + TF1; Minuit2/Migrad tolerance 0.01; x and d in metres; covariance matrices stored beside every fit");source.Write();defs.Write();fitmethod.Write();out.Close();std::cout<<"Wrote event analysis for "<<cells.size()<<" cells"<<std::endl;
}
