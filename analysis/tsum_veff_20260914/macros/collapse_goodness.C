#include <TFile.h>
#include <TMath.h>
#include <TNamed.h>
#include <TTree.h>

#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {std::vector<std::string> sp(const std::string&s){std::stringstream ss(s);std::string q;std::vector<std::string>v;while(std::getline(ss,q,','))v.push_back(q);return v;}}
void collapse_goodness(){std::ifstream in("sources/propagation_collapse.csv");std::string line;std::getline(in,line);std::map<std::string,std::pair<double,int>> sums;while(std::getline(in,line)){auto v=sp(line);if(v.size()<12||v[4]=="0")continue;std::string k=v[0]+","+v[1];double z=std::stod(v[11]);sums[k].first+=z*z;sums[k].second++;}TFile out("sources/collapse_goodness.root","RECREATE");TTree t("collapse","Equal-distance LEFT/RIGHT collapse goodness excluding correlated x=0 pair");char material[16]={},estimator[24]={};double chi2=0,p=0;int ndf=0;t.Branch("material",material,"material/C");t.Branch("estimator",estimator,"estimator/C");t.Branch("chi2",&chi2);t.Branch("ndf",&ndf);t.Branch("probability",&p);std::ofstream csv("sources/collapse_goodness.csv");csv<<std::setprecision(12)<<"material,estimator,chi2,ndf,probability\n";for(auto&kv:sums){auto q=sp(kv.first);snprintf(material,sizeof(material),"%s",q[0].c_str());snprintf(estimator,sizeof(estimator),"%s",q[1].c_str());chi2=kv.second.first;ndf=kv.second.second;p=TMath::Prob(chi2,ndf);t.Fill();csv<<q[0]<<','<<q[1]<<','<<chi2<<','<<ndf<<','<<p<<'\n';}t.Write();TNamed note("method","Sum of squared equal-distance LEFT-minus-RIGHT z values; x=0 excluded because its two sides share events and covariance is not stored");note.Write();out.Close();}
