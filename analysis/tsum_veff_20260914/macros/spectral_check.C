#include <TFile.h>
#include <TH1D.h>
#include <TNamed.h>
#include <TTree.h>

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
struct Input { std::string id,mat,path; int x=0,n=0; };
std::vector<std::string> split(const std::string&s,char sep){std::stringstream ss(s);std::string q;std::vector<std::string>v;while(std::getline(ss,q,sep))v.push_back(q);return v;}
std::vector<Input> inputs(){std::ifstream f("../../presentations/v9/sources/input_roots.tsv");if(!f)throw std::runtime_error("missing inventory");std::string s;std::getline(f,s);std::vector<Input>v;while(std::getline(f,s)){auto q=split(s,'\t');Input z;z.id=q.at(0);z.mat=q.at(1);z.x=std::stoi(q.at(2));z.n=std::stoi(q.at(3));z.path=q.at(4);v.push_back(z);}return v;}
std::string safe(std::string s){for(char&c:s)if(c=='-')c='_';return s;}
std::string key(const std::string&m,int dmm,const std::string&type){return safe(m)+"_d"+std::to_string(dmm)+"_"+type;}
double median(TH1D&h){double p=.5,q=0;h.GetQuantiles(1,&q,&p);return q;}
}

void spectral_check(){
  const auto files=inputs();if(files.size()!=21)throw std::runtime_error("expected 21 inputs");TFile out("sources/spectral_check.root","RECREATE");
  std::map<std::string,TH1D*> hall,hfirst;
  for(auto m:{std::string("EJ-200"),std::string("EJ-204"),std::string("EJ-230")})for(int d:{50,200,500,700,900,1200,1350}){out.cd();hall[m+std::to_string(d)]=new TH1D(key(m,d,"all_detected").c_str(),"",400,250,650);hfirst[m+std::to_string(d)]=new TH1D(key(m,d,"first_arrival").c_str(),"",400,250,650);}
  for(size_t fi=0;fi<files.size();fi++){
    const auto&z=files[fi];std::cout<<"["<<fi+1<<"/"<<files.size()<<"] "<<z.id<<std::endl;TFile f(z.path.c_str(),"READ");auto*t=(TTree*)f.Get("sipm_hits");if(!t)throw std::runtime_error("missing sipm_hits "+z.id);int ev=0,gid=0;double time=0,wl=0;t->SetBranchStatus("*",0);for(auto b:{"event_id","global_id","time_ns","wl_nm"})t->SetBranchStatus(b,1);t->SetBranchAddress("event_id",&ev);t->SetBranchAddress("global_id",&gid);t->SetBranchAddress("time_ns",&time);t->SetBranchAddress("wl_nm",&wl);t->SetCacheSize(512LL*1024*1024);for(auto b:{"event_id","global_id","time_ns","wl_nm"})t->AddBranchToCache(b,true);
    std::vector<double>tl(z.n,std::numeric_limits<double>::infinity()),tr(z.n,std::numeric_limits<double>::infinity()),wlL(z.n,0),wlR(z.n,0);int dl=700+z.x,dr=700-z.x;auto*haL=hall[z.mat+std::to_string(dl)],*haR=hall[z.mat+std::to_string(dr)];
    for(Long64_t j=0;j<t->GetEntries();j++){t->GetEntry(j);if(ev<0||ev>=z.n)throw std::runtime_error("bad event");if(gid>=0&&gid<8){haL->Fill(wl);if(time<tl[ev]){tl[ev]=time;wlL[ev]=wl;}}else if(gid>=8&&gid<16){haR->Fill(wl);if(time<tr[ev]){tr[ev]=time;wlR[ev]=wl;}}}
    auto*hfL=hfirst[z.mat+std::to_string(dl)],*hfR=hfirst[z.mat+std::to_string(dr)];for(int e=0;e<z.n;e++){if(!std::isfinite(tl[e])||!std::isfinite(tr[e]))throw std::runtime_error("missing END first hit");hfL->Fill(wlL[e]);hfR->Fill(wlR[e]);}
  }
  std::ofstream csv("sources/spectral_summary.csv");csv<<std::setprecision(12)<<"material,distance_mm,N_all,mean_wl_all_nm,sem_wl_all_nm,median_wl_all_nm,N_first,mean_wl_first_nm,sem_wl_first_nm,median_wl_first_nm\n";for(auto m:{std::string("EJ-200"),std::string("EJ-204"),std::string("EJ-230")})for(int d:{50,200,500,700,900,1200,1350}){auto*a=hall[m+std::to_string(d)],*q=hfirst[m+std::to_string(d)];csv<<m<<','<<d<<','<<a->GetEntries()<<','<<a->GetMean()<<','<<a->GetRMS()/std::sqrt(a->GetEntries())<<','<<median(*a)<<','<<q->GetEntries()<<','<<q->GetMean()<<','<<q->GetRMS()/std::sqrt(q->GetEntries())<<','<<median(*q)<<'\n';}
  std::ofstream tests("sources/spectral_short_long_tests.csv");tests<<std::setprecision(12)<<"material,population,short_distance_mm,long_distance_mm,KS_p,mean_shift_long_minus_short_nm,median_shift_long_minus_short_nm\n";for(auto m:{std::string("EJ-200"),std::string("EJ-204"),std::string("EJ-230")})for(auto type:{std::string("all_detected"),std::string("first_arrival")}){auto*a=type=="all_detected"?hall[m+"50"]:hfirst[m+"50"],*b=type=="all_detected"?hall[m+"1350"]:hfirst[m+"1350"];tests<<m<<','<<type<<",50,1350,"<<a->KolmogorovTest(b)<<','<<b->GetMean()-a->GetMean()<<','<<median(*b)-median(*a)<<'\n';}
  out.cd();for(auto&x:hall)x.second->Write();for(auto&x:hfirst)x.second->Write();TNamed src("source","sipm_hits from the 21 current corrected-transport production ROOTs; detected wl_nm only");TNamed limits("limits","No accumulated optical path, boundary-count, track id, or creation wavelength is stored with sipm_hits; detected spectral evolution can be tested but microscopic decomposition cannot");src.Write();limits.Write();out.Close();std::cout<<"Wrote detected-wavelength check"<<std::endl;
}
