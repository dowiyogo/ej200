#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
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
constexpr int K=20;
double quantile(const std::vector<double>&s,double p){double z=p*(s.size()-1),a=std::floor(z),u=z-a;size_t n=a;return n+1<s.size()?s[n]*(1-u)+s[n+1]*u:s[n];}
std::string id(const std::string&m,int x){std::string p=m=="EJ-200"?"EJ200_":m=="EJ-204"?"EJ204_":"EJ230_";return p+(x<0?"xm"+std::to_string(-x):"xp"+std::to_string(x));}
}

void mirror_diagnostics(){
  TFile input("../../presentations/v9/sources/timing_events.root","READ");auto*t=(TTree*)input.Get("timing_events");if(!t)throw std::runtime_error("timing_events missing");char cell[32]={};double l[K],r[K];t->SetBranchAddress("cell_id",cell);t->SetBranchAddress("t_left_ns",l);t->SetBranchAddress("t_right_ns",r);std::map<std::string,std::vector<double>> data;for(Long64_t j=0;j<t->GetEntries();j++){t->GetEntry(j);data[cell].push_back(.5*(l[0]+r[0]));}
  TFile out("sources/mirror_diagnostics.root","RECREATE");std::ofstream csv("sources/mirror_diagnostics.csv");csv<<std::setprecision(12)<<"material,abs_x_mm,centering,N_minus,N_plus,center_minus_ns,center_plus_ns,hist_low_ns,hist_high_ns,bin_width_ps,nbins,KS_p,AD_p,integral_minus,integration_plus,max_abs_density_difference\n";
  for(auto mat:{"EJ-200","EJ-204","EJ-230"})for(int a:{200,500,650})for(auto centering:{"median","mean"}){auto vn=data.at(id(mat,-a)),vp=data.at(id(mat,a));auto center=[&](std::vector<double>v){if(std::string(centering)=="mean")return std::accumulate(v.begin(),v.end(),0.0)/v.size();std::sort(v.begin(),v.end());return quantile(v,.5);};double cn=center(vn),cp=center(vp);for(double&z:vn)z-=cn;for(double&z:vp)z-=cp;std::string key=std::string(mat)+"_"+std::to_string(a)+"_"+centering;std::replace(key.begin(),key.end(),'-','_');auto*dir=out.mkdir(key.c_str());dir->cd();TH1D hn("minus_counts","",160,-.32,.32),hp("plus_counts","",160,-.32,.32);hn.Sumw2();hp.Sumw2();for(double z:vn)hn.Fill(z);for(double z:vp)hp.Fill(z);double ks=hn.KolmogorovTest(&hp),ad=hn.AndersonDarlingTest(&hp);auto*nn=(TH1D*)hn.Clone("minus_density");auto*np=(TH1D*)hp.Clone("plus_density");nn->Scale(1.0/nn->Integral("width"));np->Scale(1.0/np->Integral("width"));auto*diff=(TH1D*)np->Clone("plus_minus_difference");diff->Add(nn,-1);double md=0;for(int b=1;b<=diff->GetNbinsX();b++)md=std::max(md,std::abs(diff->GetBinContent(b)));hn.Write();hp.Write();nn->Write();np->Write();diff->Write();csv<<mat<<','<<a<<','<<centering<<','<<vn.size()<<','<<vp.size()<<','<<cn<<','<<cp<<','<<-.32<<','<<.32<<','<<4<<','<<160<<','<<ks<<','<<ad<<','<<hn.Integral()<<','<<hp.Integral()<<','<<md<<'\n';}
  out.cd();TNamed method("method","T0 centered separately by median or mean; identical [-0.32,+0.32] ns range, 160 bins (4 ps), raw-count ROOT Kolmogorov and Anderson-Darling tests; density histograms and plus-minus difference stored");method.Write();out.Close();std::cout<<"Wrote 18 mirror shape diagnostics"<<std::endl;
}
