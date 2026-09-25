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

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr const char* kInput =
    "/home/rrios/ej200/analysis/tsum_veff_20260914/sources/tsum_veff.root";
constexpr const char* kTreeName = "derived_events";
constexpr int kExpectedEntries = 210000;
constexpr int kExpectedCells = 21;
constexpr int kExpectedEventsPerCell = 10000;
constexpr double kHistogramBinWidthNs = 0.004;
constexpr double kHistogramHalfRangeNs = 0.320;
constexpr double kIqrGaussianScale = 1.349;
constexpr double kIqrLowQuantile = 0.25;
constexpr double kMedianQuantile = 0.50;
constexpr double kIqrHighQuantile = 0.75;
constexpr double kRobustLowQuantile = 0.16;
constexpr double kRobustHighQuantile = 0.84;
constexpr double kFitHalfWidthRobustSigma = 2.0;
constexpr double kUnreliableChi2Ndf = 5.0;
constexpr double kGuardSigmaMultiplier = 3.0;
constexpr double kCenterToleranceNs = 1e-12;
constexpr int kProfileBins = 40;
constexpr double kProfilePaddingPe = 0.5;
constexpr double kRatioThreshold = 0.5;
constexpr double kMmPerM = 1000.0;
constexpr std::array<double, 3> kReportedEvenA2NsPerM2 = {
    -0.0718094078599, -0.0299287544695, 0.00265027733139};
constexpr std::array<double, 3> kReportedEvenA2ErrorNsPerM2 = {
    0.00179749428654, 0.00179473390186, 0.00174633304245};

struct Event { double left=0,right=0,t0=0,npe=0; };
struct Cell { std::string id,material; int x_mm=0; std::vector<Event> events; };
struct BasicWidth { double mean=0,rms=0,iqr=0,q68=0,median=0; };
struct GaussianResult {
  double sigma_ns=0,sigma_error_ns=0,mean_ns=0,mean_error_ns=0;
  double fit_low_ns=0,fit_high_ns=0,chi2=0; int ndf=0,status=-1,covariance_status=-1;
  bool unreliable=true; TMatrixDSym covariance{3};
};
struct ProfileResult {
  double intercept_ns=0,intercept_error_ns=0,slope_ns_per_pe=0;
  double slope_error_ns_per_pe=0,chi2=0; int ndf=0,status=-1,covariance_status=-1;
  bool unreliable=true; TMatrixDSym covariance{2};
};
struct CellResult {
  BasicWidth left,right,t0; GaussianResult fit_left,fit_right,fit_t0;
  double guard_gap_ns=0,guard_uncertainty_ns=0; bool guard_pass=false;
  double rho_t0_npe=0,mean_npe=0,sem_npe=0,mean_t0=0,sem_t0=0;
  ProfileResult profile;
};

double quantile(const std::vector<double>& sorted,double probability){
  double position=probability*(sorted.size()-1); size_t low=std::floor(position);
  size_t high=std::min(low+1,sorted.size()-1); double fraction=position-low;
  return sorted[low]*(1-fraction)+sorted[high]*fraction;
}
BasicWidth basicWidth(std::vector<double> values){
  BasicWidth r; r.mean=std::accumulate(values.begin(),values.end(),0.)/values.size();
  double ss=0; for(double v:values)ss+=(v-r.mean)*(v-r.mean); r.rms=std::sqrt(ss/values.size());
  std::sort(values.begin(),values.end()); r.median=quantile(values,kMedianQuantile);
  r.iqr=(quantile(values,kIqrHighQuantile)-quantile(values,kIqrLowQuantile))/kIqrGaussianScale;
  r.q68=.5*(quantile(values,kRobustHighQuantile)-quantile(values,kRobustLowQuantile)); return r;
}
double pearson(const std::vector<double>& a,const std::vector<double>& b){
  double ma=std::accumulate(a.begin(),a.end(),0.)/a.size(),mb=std::accumulate(b.begin(),b.end(),0.)/b.size();
  double sa=0,sb=0,sab=0; for(size_t i=0;i<a.size();++i){double x=a[i]-ma,y=b[i]-mb;sa+=x*x;sb+=y*y;sab+=x*y;} return sab/std::sqrt(sa*sb);
}
double sem(const std::vector<double>& v,double mean){double ss=0;for(double x:v)ss+=(x-mean)*(x-mean);return std::sqrt(ss/(v.size()*(v.size()-1.)));}
GaussianResult fitGaussian(TH1D& h,const BasicWidth& w,const std::string& name){
  GaussianResult r; r.fit_low_ns=w.median-kFitHalfWidthRobustSigma*w.q68;r.fit_high_ns=w.median+kFitHalfWidthRobustSigma*w.q68;
  double peak=h.GetBinCenter(h.GetMaximumBin());TF1 f(name.c_str(),"gaus",r.fit_low_ns,r.fit_high_ns);f.SetParameters(h.GetMaximum(),peak,w.q68);
  TFitResultPtr fr=h.Fit(&f,"QNRS0");r.status=int(fr);r.covariance_status=fr.Get()?fr->CovMatrixStatus():-1;r.mean_ns=f.GetParameter(1);r.mean_error_ns=f.GetParError(1);r.sigma_ns=std::abs(f.GetParameter(2));r.sigma_error_ns=f.GetParError(2);r.chi2=f.GetChisquare();r.ndf=f.GetNDF();
  r.unreliable=r.status!=0||r.covariance_status<2||r.ndf<=0||r.chi2/r.ndf>kUnreliableChi2Ndf;if(fr.Get())r.covariance=fr->GetCovarianceMatrix();f.Write();r.covariance.Write((name+"_covariance").c_str());return r;
}
double pairedGuardUncertainty(const Cell& c,const BasicWidth& l,const BasicWidth& r,const BasicWidth& t){
  double ss=0;for(const auto&e:c.events){double il=((e.left-l.mean)*(e.left-l.mean)-l.rms*l.rms)/(2*l.rms);double ir=((e.right-r.mean)*(e.right-r.mean)-r.rms*r.rms)/(2*r.rms);double it=((e.t0-t.mean)*(e.t0-t.mean)-t.rms*t.rms)/(2*t.rms);double ig=it-.5*(il+ir);ss+=ig*ig;}return std::sqrt(ss)/c.events.size();
}
ProfileResult fitProfile(TProfile& p,double low,double high,const std::string& name){
  ProfileResult r;TF1 f(name.c_str(),"pol1",low,high);TFitResultPtr fr=p.Fit(&f,"QNRS0");r.status=int(fr);r.covariance_status=fr.Get()?fr->CovMatrixStatus():-1;r.intercept_ns=f.GetParameter(0);r.intercept_error_ns=f.GetParError(0);r.slope_ns_per_pe=f.GetParameter(1);r.slope_error_ns_per_pe=f.GetParError(1);r.chi2=f.GetChisquare();r.ndf=f.GetNDF();r.unreliable=r.status!=0||r.covariance_status<2||r.ndf<=0||r.chi2/r.ndf>kUnreliableChi2Ndf;if(fr.Get())r.covariance=fr->GetCovarianceMatrix();f.Write();r.covariance.Write((name+"_covariance").c_str());return r;
}
std::string cellId(const std::string&m,int x){std::string p=m=="EJ-200"?"EJ200_":m=="EJ-204"?"EJ204_":"EJ230_";return p+(x<0?"xm":"xp")+std::to_string(std::abs(x));}
std::string safe(std::string s){for(char&c:s)if(c=='-')c='_';return s;}
struct PositionFit{double a0=0,a2=0,a4=0,a2_error=0,a4_error=0,chi2=0;int ndf=0,status=-1,covariance_status=-1;};
PositionFit fitEven(const std::vector<double>&x,const std::vector<double>&y,const std::vector<double>&ey,bool quartic,const std::string&name){
  std::vector<double>ex(x.size(),0.);TGraphErrors g(x.size(),x.data(),y.data(),ex.data(),ey.data());g.SetName((name+"_graph").c_str());TF1 f(name.c_str(),quartic?"[0]+[1]*x*x+[2]*x*x*x*x":"[0]+[1]*x*x",-.66,.66);f.SetParameter(0,y[x.size()/2]);f.SetParameter(1,0);if(quartic)f.SetParameter(2,0);TFitResultPtr fr=g.Fit(&f,"QNRS0");PositionFit r;r.a0=f.GetParameter(0);r.a2=f.GetParameter(1);r.a2_error=f.GetParError(1);if(quartic){r.a4=f.GetParameter(2);r.a4_error=f.GetParError(2);}r.chi2=f.GetChisquare();r.ndf=f.GetNDF();r.status=int(fr);r.covariance_status=fr.Get()?fr->CovMatrixStatus():-1;g.Write();f.Write();if(fr.Get())fr->GetCovarianceMatrix().Write((name+"_covariance").c_str());return r;
}
}

void order_stat_weight(){
  TFile input(kInput,"READ");if(input.IsZombie())throw std::runtime_error("Cannot open input ROOT read-only");auto*tree=dynamic_cast<TTree*>(input.Get(kTreeName));if(!tree||tree->GetEntries()!=kExpectedEntries)throw std::runtime_error("Unexpected derived_events tree or count");
  char cid[64]={0},mat[32]={0};int x=0,npe=0;double l=0,r=0,t0=0;tree->SetBranchAddress("cell_id",cid);tree->SetBranchAddress("material",mat);tree->SetBranchAddress("x_mm",&x);tree->SetBranchAddress("t_left_ns",&l);tree->SetBranchAddress("t_right_ns",&r);tree->SetBranchAddress("t0_ns",&t0);tree->SetBranchAddress("npe_end",&npe);
  std::map<std::string,Cell>cells;for(Long64_t i=0;i<tree->GetEntries();++i){tree->GetEntry(i);auto&c=cells[cid];c.id=cid;c.material=mat;c.x_mm=x;c.events.push_back({l,r,t0,double(npe)});}if(cells.size()!=kExpectedCells)throw std::runtime_error("Expected 21 cells");for(const auto&z:cells)if(z.second.events.size()!=kExpectedEventsPerCell)throw std::runtime_error("Expected 10000 events in "+z.first);
  TFile output("sources/order_stat_weight.root","RECREATE");std::ofstream widths("sources/part_a_widths.csv"),guards("sources/part_a_guardrails.csv"),center("sources/part_a_x0_sanity.csv"),profiles("sources/part_b_profiles.csv");
  widths<<std::setprecision(12)<<"cell_id,material,x_mm,N,observable,RMS_ns,sigma_IQR_ns,sigma_G_ns,sigma_G_err_ns,fit_mean_ns,fit_mean_err_ns,fit_lo_ns,fit_hi_ns,chi2,ndf,chi2_ndf,fit_status,covariance_status,fit_reliable,ratio_to_T0_IQR,ratio_to_T0_G\n";
  guards<<std::setprecision(12)<<"cell_id,material,x_mm,N,RMS_tL_ns,RMS_tR_ns,RMS_T0_ns,triangle_bound_ns,gap_ns,paired_combined_uncertainty_ns,three_sigma_margin_ns,status\n";
  center<<std::setprecision(12)<<"cell_id,material,method,sigma_tL_ns,sigma_tR_ns,sigma_T0_ns,T0_le_tL,T0_le_tR,status\n";
  profiles<<std::setprecision(12)<<"cell_id,material,x_mm,N,mean_Npe_END,sem_Npe_END,mean_T0_ns,sem_T0_ns,rho_T0_Npe_END,profile_bins,fit_low_Npe,fit_high_Npe,intercept_ns,intercept_err_ns,slope_ns_per_pe,slope_err_ns_per_pe,chi2,ndf,chi2_ndf,fit_status,covariance_status,fit_reliable\n";
  const std::vector<std::string>materials={"EJ-200","EJ-204","EJ-230"};const std::vector<int>positions={-650,-500,-200,0,200,500,650};std::map<std::string,CellResult>results;auto*cellsdir=output.mkdir("cells");
  for(const auto&m:materials)for(int xpos:positions){std::string id=cellId(m,xpos);const Cell&c=cells.at(id);std::vector<double>vl,vr,vt,vn;for(const auto&e:c.events){vl.push_back(e.left);vr.push_back(e.right);vt.push_back(e.t0);vn.push_back(e.npe);}auto&z=results[id];z.left=basicWidth(vl);z.right=basicWidth(vr);z.t0=basicWidth(vt);z.mean_npe=std::accumulate(vn.begin(),vn.end(),0.)/vn.size();z.sem_npe=sem(vn,z.mean_npe);z.mean_t0=z.t0.mean;z.sem_t0=sem(vt,z.mean_t0);z.rho_t0_npe=pearson(vt,vn);
    cellsdir->cd();auto*d=cellsdir->mkdir(id.c_str());d->cd();int nb=std::lround(2*kHistogramHalfRangeNs/kHistogramBinWidthNs);TH1D hl("hTL",(id+";t_{L} [ns];events").c_str(),nb,z.left.median-kHistogramHalfRangeNs,z.left.median+kHistogramHalfRangeNs),hr("hTR",(id+";t_{R} [ns];events").c_str(),nb,z.right.median-kHistogramHalfRangeNs,z.right.median+kHistogramHalfRangeNs),ht("hT0",(id+";T_{0} [ns];events").c_str(),nb,z.t0.median-kHistogramHalfRangeNs,z.t0.median+kHistogramHalfRangeNs);std::string sel="cell_id==\""+id+"\"";if(tree->Draw("t_left_ns>>hTL",sel.c_str(),"goff")!=kExpectedEventsPerCell||tree->Draw("t_right_ns>>hTR",sel.c_str(),"goff")!=kExpectedEventsPerCell||tree->Draw("t0_ns>>hT0",sel.c_str(),"goff")!=kExpectedEventsPerCell)throw std::runtime_error("TTree::Draw count failure in "+id);z.fit_left=fitGaussian(hl,z.left,"fit_tL");z.fit_right=fitGaussian(hr,z.right,"fit_tR");z.fit_t0=fitGaussian(ht,z.t0,"fit_T0");hl.Write();hr.Write();ht.Write();
    z.guard_gap_ns=z.t0.rms-.5*(z.left.rms+z.right.rms);z.guard_uncertainty_ns=pairedGuardUncertainty(c,z.left,z.right,z.t0);z.guard_pass=z.guard_gap_ns<=kGuardSigmaMultiplier*z.guard_uncertainty_ns;guards<<id<<','<<m<<','<<xpos<<','<<c.events.size()<<','<<z.left.rms<<','<<z.right.rms<<','<<z.t0.rms<<','<<.5*(z.left.rms+z.right.rms)<<','<<z.guard_gap_ns<<','<<z.guard_uncertainty_ns<<','<<kGuardSigmaMultiplier*z.guard_uncertainty_ns<<','<<(z.guard_pass?"PASS":"FAIL")<<'\n';if(!z.guard_pass){output.Write();output.Close();throw std::runtime_error("Corrected L2 guardrail failed in "+id);}
    auto emit=[&](const char*o,const BasicWidth&w,const GaussianResult&f){double ri=std::string(o)=="T0"?1:w.iqr/z.t0.iqr,rg=std::string(o)=="T0"?1:f.sigma_ns/z.fit_t0.sigma_ns;widths<<id<<','<<m<<','<<xpos<<','<<c.events.size()<<','<<o<<','<<w.rms<<','<<w.iqr<<','<<f.sigma_ns<<','<<f.sigma_error_ns<<','<<f.mean_ns<<','<<f.mean_error_ns<<','<<f.fit_low_ns<<','<<f.fit_high_ns<<','<<f.chi2<<','<<f.ndf<<','<<(f.ndf>0?f.chi2/f.ndf:-1)<<','<<f.status<<','<<f.covariance_status<<','<<(f.unreliable?"NO":"YES")<<','<<ri<<','<<rg<<'\n';};emit("tL",z.left,z.fit_left);emit("tR",z.right,z.fit_right);emit("T0",z.t0,z.fit_t0);
    if(xpos==0){auto ec=[&](const char*method,double sl,double sr,double st){bool pl=st<=sl+kCenterToleranceNs,pr=st<=sr+kCenterToleranceNs;center<<id<<','<<m<<','<<method<<','<<sl<<','<<sr<<','<<st<<','<<(pl?"YES":"NO")<<','<<(pr?"YES":"NO")<<','<<(pl&&pr?"PASS":"FLAG")<<'\n';};ec("RMS",z.left.rms,z.right.rms,z.t0.rms);ec("IQR",z.left.iqr,z.right.iqr,z.t0.iqr);ec("Gaussian",z.fit_left.sigma_ns,z.fit_right.sigma_ns,z.fit_t0.sigma_ns);}
    auto mm=std::minmax_element(vn.begin(),vn.end());double plo=*mm.first-kProfilePaddingPe,phi=*mm.second+kProfilePaddingPe;TProfile p("T0_vs_Npe_END",(id+";N_{pe}^{END};<T_{0}> [ns]").c_str(),kProfileBins,plo,phi);for(size_t i=0;i<vn.size();++i)p.Fill(vn[i],vt[i]);z.profile=fitProfile(p,plo,phi,"fit_T0_vs_Npe_END");p.Write();profiles<<id<<','<<m<<','<<xpos<<','<<c.events.size()<<','<<z.mean_npe<<','<<z.sem_npe<<','<<z.mean_t0<<','<<z.sem_t0<<','<<z.rho_t0_npe<<','<<kProfileBins<<','<<plo<<','<<phi<<','<<z.profile.intercept_ns<<','<<z.profile.intercept_error_ns<<','<<z.profile.slope_ns_per_pe<<','<<z.profile.slope_error_ns_per_pe<<','<<z.profile.chi2<<','<<z.profile.ndf<<','<<(z.profile.ndf>0?z.profile.chi2/z.profile.ndf:-1)<<','<<z.profile.status<<','<<z.profile.covariance_status<<','<<(z.profile.unreliable?"NO":"YES")<<'\n';}
  std::ofstream threshold("sources/part_a_thresholds.csv");threshold<<std::setprecision(12)<<"material,method,threshold,below_x_mm,neighbor_x_mm,ratio_below,ratio_neighbor,linear_interpolated_crossing_x_mm,interpretation\n";for(const auto&m:materials)for(const std::string method:{"IQR","Gaussian"}){std::vector<std::pair<int,double>>pts;for(int xpos:positions){const auto&z=results.at(cellId(m,xpos));pts.push_back({xpos,method=="IQR"?z.left.iqr/z.t0.iqr:z.fit_left.sigma_ns/z.fit_t0.sigma_ns});}int bi=-1;for(int i=int(pts.size())-1;i>=0;--i)if(pts[i].second<kRatioThreshold)bi=i;if(bi<0||bi+1>=int(pts.size())){threshold<<m<<','<<method<<','<<kRatioThreshold<<",NA,NA,NA,NA,NA,no sampled crossing\n";continue;}auto below=pts[bi],next=pts[bi+1];double cross=below.first+(kRatioThreshold-below.second)*(next.first-below.first)/(next.second-below.second);threshold<<m<<','<<method<<','<<kRatioThreshold<<','<<below.first<<','<<next.first<<','<<below.second<<','<<next.second<<','<<cross<<",linear interpolation between sparse sampled positions\n";}
  std::ofstream cp("sources/part_b_chain_rule_points.csv"),cs("sources/part_b_chain_rule_summary.csv");cp<<std::setprecision(12)<<"material,x_mm,mean_Npe_END,even_mean_Npe_END,slope_ns_per_pe,even_slope_ns_per_pe,observed_T0_ns,observed_even_shift_from_x0_ns,predicted_even_shift_from_x0_ns\n";cs<<std::setprecision(12)<<"material,reported_even_quadratic_a2_ns_per_m2,reported_a2_err_ns_per_m2,recomputed_a2_ns_per_m2,recomputed_a2_err_ns_per_m2,predicted_a2_ns_per_m2,fraction_explained,predicted_quartic_a2_ns_per_m2,predicted_quartic_a4_ns_per_m4,observed_quartic_a2_ns_per_m2,observed_quartic_a4_ns_per_m4,observed_quadratic_chi2,observed_quadratic_ndf,predicted_quadratic_chi2,predicted_quadratic_ndf,observed_quartic_chi2,observed_quartic_ndf,predicted_quartic_chi2,predicted_quartic_ndf\n";auto*chaindir=output.mkdir("chain_rule");
  for(size_t mi=0;mi<materials.size();++mi){const auto&m=materials[mi];std::map<int,double>en,eb,pred;const auto&zc=results.at(cellId(m,0));en[0]=zc.mean_npe;eb[0]=zc.profile.slope_ns_per_pe;pred[0]=0;for(int a:{200,500,650}){const auto&zm=results.at(cellId(m,-a));const auto&zp=results.at(cellId(m,a));en[a]=.5*(zm.mean_npe+zp.mean_npe);eb[a]=.5*(zm.profile.slope_ns_per_pe+zp.profile.slope_ns_per_pe);}int prev=0;for(int a:{200,500,650}){pred[a]=pred[prev]+.5*(eb[prev]+eb[a])*(en[a]-en[prev]);prev=a;}std::vector<double>xx,obs,prd,err;for(int xpos:positions){const auto&z=results.at(cellId(m,xpos));int a=std::abs(xpos);double oe=xpos==0?zc.mean_t0:.5*(results.at(cellId(m,-a)).mean_t0+results.at(cellId(m,a)).mean_t0);xx.push_back(xpos/kMmPerM);obs.push_back(z.mean_t0);prd.push_back(pred[a]);err.push_back(z.sem_t0);cp<<m<<','<<xpos<<','<<z.mean_npe<<','<<en[a]<<','<<z.profile.slope_ns_per_pe<<','<<eb[a]<<','<<z.mean_t0<<','<<oe-zc.mean_t0<<','<<pred[a]<<'\n';}chaindir->cd();auto oq=fitEven(xx,obs,err,false,"observed_even_quadratic_"+safe(m)),pq=fitEven(xx,prd,err,false,"predicted_even_quadratic_"+safe(m)),o4=fitEven(xx,obs,err,true,"observed_even_quartic_"+safe(m)),p4=fitEven(xx,prd,err,true,"predicted_even_quartic_"+safe(m));cs<<m<<','<<kReportedEvenA2NsPerM2[mi]<<','<<kReportedEvenA2ErrorNsPerM2[mi]<<','<<oq.a2<<','<<oq.a2_error<<','<<pq.a2<<','<<pq.a2/kReportedEvenA2NsPerM2[mi]<<','<<p4.a2<<','<<p4.a4<<','<<o4.a2<<','<<o4.a4<<','<<oq.chi2<<','<<oq.ndf<<','<<pq.chi2<<','<<pq.ndf<<','<<o4.chi2<<','<<o4.ndf<<','<<p4.chi2<<','<<p4.ndf<<'\n';}
  output.cd();TNamed status("analysis_status","COMPLETE"),source("source","tsum_veff.root opened READ; 210000 events; 21 cells; no simulation"),wm("width_method","RMS and IQR/1.349; Gaussian peak seed; 4 ps bins over median +/-0.320 ns; fit median +/-2*q68; unreliable if chi2/ndf>5"),gm("guard_method","L2 guard uses RMS only: RMS(T0) <= [RMS(tL)+RMS(tR)]/2 + 3 paired delta-method SE; IQR and Gaussian core widths do not obey L2 triangle inequality"),cm("chain_method","40-bin TProfile linear fit per cell; measured Npe_END even-symmetrized without exponential fit; trapezoidal integration dT=beta*dN from x=0; same position weights as measured T0");status.Write();source.Write();wm.Write();gm.Write();cm.Write();output.Write();output.Close();widths.close();guards.close();center.close();profiles.close();threshold.close();cp.close();cs.close();std::ofstream meta("sources/order_stat_weight.meta.json");meta<<"{\n  \"status\": \"COMPLETE\",\n  \"input\": \""<<kInput<<"\",\n  \"input_mode\": \"READ\",\n  \"entries\": "<<kExpectedEntries<<",\n  \"cells\": "<<kExpectedCells<<",\n  \"ratio_threshold\": "<<kRatioThreshold<<",\n  \"gaussian_unreliable_chi2_ndf\": "<<kUnreliableChi2Ndf<<",\n  \"part_b_executed\": true,\n  \"new_simulation\": false\n}\n";meta.close();std::cout<<"COMPLETE: corrected guardrail passed 21/21; Part B executed"<<std::endl;
}
