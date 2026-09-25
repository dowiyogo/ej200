#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TGraphErrors.h>
#include <TMatrixDSym.h>
#include <TNamed.h>
#include <Math/MinimizerOptions.h>

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

void extended_models(){
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2","Migrad");ROOT::Math::MinimizerOptions::SetDefaultTolerance(.01);TFile in("sources/tsum_veff.root","READ"),out("sources/extended_models.root","RECREATE");std::ofstream csv("sources/extended_models.csv");csv<<std::setprecision(12)<<"observable,material,model,status,covariance_status,chi2,ndf,chi2_ndf,probability,AIC,p0,p0_err,p2,p2_err,p4,p4_err,delta_chi2_vs_constant\n";
  for(auto slug:{std::string("EJ_200"),std::string("EJ_204"),std::string("EJ_230")})for(auto obs:{std::string("t0"),std::string("t0m5")}){if(obs=="t0m5"&&slug!="EJ_230")continue;std::string material=slug=="EJ_200"?"EJ-200":slug=="EJ_204"?"EJ-204":"EJ-230",gn=obs+"_vs_x_"+slug;auto*g=(TGraphErrors*)in.Get(("position_graphs/"+gn).c_str());auto*fc=(TF1*)in.Get(("position_fits/"+gn+"_constant").c_str());if(!g||!fc)throw std::runtime_error("missing graph");TF1 f((gn+"_even_quartic").c_str(),"[0]+[1]*x*x+[2]*x*x*x*x",-.66,.66);f.SetParameters(fc->GetParameter(0),0,0);TFitResultPtr r=g->Fit(&f,"QRS0");int status=int(r),cov=r.Get()?r->CovMatrixStatus():-1;double chi=f.GetChisquare();out.cd();f.Write();if(r.Get())r->GetCovarianceMatrix().Write((gn+"_even_quartic_covariance").c_str());csv<<obs<<','<<material<<",even_quartic,"<<status<<','<<cov<<','<<chi<<','<<f.GetNDF()<<','<<chi/f.GetNDF()<<','<<f.GetProb()<<','<<chi+6<<','<<f.GetParameter(0)<<','<<f.GetParError(0)<<','<<f.GetParameter(1)<<','<<f.GetParError(1)<<','<<f.GetParameter(2)<<','<<f.GetParError(2)<<','<<fc->GetChisquare()-chi<<'\n';}
  out.cd();TNamed note("scope","Even quartic is a diagnostic added because the preregistered even quadratic has poor goodness of fit; decision does not rely on AIC alone");note.Write();out.Close();
}
