#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TNamed.h>
#include <TTree.h>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

void summarize_decision(){
  TFile in("sources/tsum_veff.root","READ"),out("sources/decision_summary.root","RECREATE");std::ofstream csv("sources/decision_summary.csv");csv<<std::setprecision(12)<<"material,estimator,center_ns,endpoint_average_ns,endpoint_minus_center_ns,peak_to_peak_ns,constant_chi2,constant_ndf,constant_probability,odd_a1_ns_per_m,odd_a1_err_ns_per_m,odd_z,even_a2_ns_per_m2,even_a2_err_ns_per_m2,even_z\n";TTree tree("decision","Position-independence summary");char material[16]={},estimator[24]={};double center=0,endavg=0,effect=0,peak=0,cchi=0,cprob=0,a1=0,ea1=0,z1=0,a2=0,ea2=0,z2=0;int cndf=0;tree.Branch("material",material,"material/C");tree.Branch("estimator",estimator,"estimator/C");tree.Branch("center_ns",&center);tree.Branch("endpoint_average_ns",&endavg);tree.Branch("endpoint_minus_center_ns",&effect);tree.Branch("peak_to_peak_ns",&peak);tree.Branch("constant_chi2",&cchi);tree.Branch("constant_ndf",&cndf);tree.Branch("constant_probability",&cprob);tree.Branch("odd_a1_ns_per_m",&a1);tree.Branch("odd_a1_err_ns_per_m",&ea1);tree.Branch("odd_z",&z1);tree.Branch("even_a2_ns_per_m2",&a2);tree.Branch("even_a2_err_ns_per_m2",&ea2);tree.Branch("even_z",&z2);
  for(auto slug:{std::string("EJ_200"),std::string("EJ_204"),std::string("EJ_230")})for(auto est:{std::string("first_pe"),std::string("mean_first_5")}){if(est=="mean_first_5"&&slug!="EJ_230")continue;std::string mat=slug=="EJ_200"?"EJ-200":slug=="EJ_204"?"EJ-204":"EJ-230",obs=est=="first_pe"?"t0":"t0m5",gn=obs+"_vs_x_"+slug;auto*g=(TGraphErrors*)in.Get(("position_graphs/"+gn).c_str());auto*fc=(TF1*)in.Get(("position_fits/"+gn+"_constant").c_str());auto*ff=(TF1*)in.Get(("position_fits/"+gn+"_full_quadratic").c_str());auto*fe=(TF1*)in.Get(("position_fits/"+gn+"_even_quadratic").c_str());if(!g||!fc||!ff||!fe)throw std::runtime_error("missing objects");double mn=1e9,mx=-1e9,ym=0,yp=0;for(int j=0;j<g->GetN();j++){double x=g->GetX()[j],y=g->GetY()[j];mn=std::min(mn,y);mx=std::max(mx,y);if(std::abs(x)<1e-9)center=y;if(x<-.64)ym=y;if(x>.64)yp=y;}endavg=.5*(ym+yp);effect=endavg-center;peak=mx-mn;cchi=fc->GetChisquare();cndf=fc->GetNDF();cprob=fc->GetProb();a1=ff->GetParameter(1);ea1=ff->GetParError(1);z1=a1/ea1;a2=fe->GetParameter(1);ea2=fe->GetParError(1);z2=a2/ea2;snprintf(material,sizeof(material),"%s",mat.c_str());snprintf(estimator,sizeof(estimator),"%s",est.c_str());tree.Fill();csv<<mat<<','<<est<<','<<center<<','<<endavg<<','<<effect<<','<<peak<<','<<cchi<<','<<cndf<<','<<cprob<<','<<a1<<','<<ea1<<','<<z1<<','<<a2<<','<<ea2<<','<<z2<<'\n';}
  out.cd();tree.Write();TNamed note("definitions","endpoint effect = average of x=+/-650 mm means minus x=0 mean; odd coefficient from full quadratic; even coefficient from even quadratic; x in metres");note.Write();out.Close();
}
