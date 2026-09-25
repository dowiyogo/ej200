#ifndef V9_FIGURE_COMMON_C
#define V9_FIGURE_COMMON_C

#include <TArrow.h>
#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TNamed.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace v9 {
struct Table {
  std::vector<std::string> header;
  std::vector<std::map<std::string,std::string>> rows;
};

std::vector<std::string> split(const std::string& line) {
  std::vector<std::string> out; std::stringstream ss(line); std::string s;
  while(std::getline(ss,s,',')) out.push_back(s);
  return out;
}

Table readCsv(const std::string& path) {
  std::ifstream in(path); if(!in) throw std::runtime_error("Cannot open "+path);
  Table t; std::string line; std::getline(in,line); t.header=split(line);
  while(std::getline(in,line)) {
    if(line.empty()) continue; auto fields=split(line); std::map<std::string,std::string> row;
    for(size_t i=0;i<t.header.size() && i<fields.size();++i) row[t.header[i]]=fields[i];
    t.rows.push_back(row);
  }
  return t;
}

double d(const std::map<std::string,std::string>& r,const std::string& key){return std::stod(r.at(key));}
int i(const std::map<std::string,std::string>& r,const std::string& key){return std::stoi(r.at(key));}

void style() {
  gStyle->SetOptStat(0); gStyle->SetOptFit(0); gStyle->SetTitleFont(42,"XYZ");
  gStyle->SetLabelFont(42,"XYZ"); gStyle->SetLegendFont(42); gStyle->SetTextFont(42);
  gStyle->SetTitleSize(0.050,"XYZ"); gStyle->SetLabelSize(0.043,"XYZ");
  gStyle->SetPadLeftMargin(0.13); gStyle->SetPadBottomMargin(0.13); gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin(0.04); gStyle->SetLineWidth(2); gStyle->SetEndErrorSize(4);
}

int color(const std::string& material) {
  if(material=="EJ-200") return kBlue+1;
  if(material=="EJ-204") return kGreen+2;
  return kRed+1;
}
int marker(const std::string& material) {
  if(material=="EJ-200") return 20;
  if(material=="EJ-204") return 21;
  return 22;
}

void writeMeta(const std::string& name,const std::string& question,const std::string& method) {
  std::ofstream o("figures/"+name+".meta.json");
  o << "{\n  \"figure\": \""<<name<<"\",\n"
    << "  \"created_utc\": \"2026-09-14\",\n"
    << "  \"question\": \""<<question<<"\",\n"
    << "  \"inputs\": [\"sources/timing_summary.csv\", \"sources/material_summary.csv\", \"sources/timing_events.root\"],\n"
    << "  \"configuration\": \"EndTop; 16 END and 70 TOP SiPMs; vertical 1 GeV mu-; SPTR=0; no electronics; N=10000 per cell; corrected current transport\",\n"
    << "  \"method\": \""<<method<<"\",\n"
    << "  \"generator\": \"CERN ROOT macro macros/"<<name<<".C\"\n}\n";
}

void provenance(TFile& f,const std::string& method) {
  f.cd();
  TNamed source("source","sources/timing_events.root derived from the 21 current corrected-transport ROOT cells");
  TNamed config("configuration","EndTop; 16 END + 70 TOP SiPMs; vertical 1 GeV mu-; SPTR=0; no electronics; N=10000/cell");
  TNamed m("method",method.c_str()); source.Write(); config.Write(); m.Write();
}

void makeGeometry() {
  style(); TCanvas c("geometry_current","",1100,620); c.Range(0,0,1,1);
  TBox bar(0.12,0.36,0.88,0.62); bar.SetFillColor(kAzure-9);bar.SetLineColor(kBlue+2);bar.SetLineWidth(3);bar.Draw();
  for(int j=0;j<8;j++){double y=.37+j*.03;TBox*l=new TBox(.085,y,.115,y+.022);l->SetFillColor(kOrange+1);l->Draw();TBox*r=new TBox(.885,y,.915,y+.022);r->SetFillColor(kOrange+1);r->Draw();}
  for(int j=0;j<14;j++){double x=.145+j*.052;for(int q=0;q<5;q++){TBox*t=new TBox(x,.625+q*.023,x+.031,.642+q*.023);t->SetFillColor(kGreen+1);t->SetLineColor(kGreen+3);t->Draw();}}
  TArrow mu(.50,.93,.50,.67,.025,"|>");mu.SetLineWidth(4);mu.SetLineColor(kMagenta+2);mu.SetFillColor(kMagenta+2);mu.Draw();
  TLatex tx;tx.SetTextAlign(22);tx.SetTextFont(42);tx.SetTextSize(.046);tx.DrawLatex(.50,.49,"scintillator: 1400 #times 60 #times 10 mm^{3}");
  tx.SetTextSize(.040);tx.DrawLatex(.10,.30,"8 END");tx.DrawLatex(.90,.30,"8 END");tx.DrawLatex(.50,.79,"70 TOP SiPMs");tx.DrawLatex(.56,.91,"1 GeV #mu^{-}");
  tx.SetTextSize(.034);tx.DrawLatex(.50,.19,"x = 0, #pm200, #pm500, #pm650 mm   |   optical transport only");
  c.SaveAs("figures/geometry_current.pdf");TFile f("figures/geometry_current.root","RECREATE");c.Write();provenance(f,"ROOT primitives; geometry values from verified per-cell configuration");f.Close();writeMeta("geometry_current","What detector geometry is evaluated?","Schematic drawn with ROOT primitives; dimensions are physical configuration values");
}

void makeNpeVsX() {
  style();auto t=readCsv("sources/timing_summary.csv");TCanvas c("npe_vs_x","",920,640);TMultiGraph mg;TLegend leg(.16,.69,.38,.89);leg.SetBorderSize(0);leg.SetFillStyle(0);
  TFile f("figures/npe_vs_x.root","RECREATE");
  for(auto mat:{"EJ-200","EJ-204","EJ-230"}){std::vector<double>x,y,ex,ey;for(auto&r:t.rows)if(r.at("material")==mat){x.push_back(d(r,"x_mm"));y.push_back(d(r,"npe_end"));ex.push_back(0);ey.push_back(d(r,"npe_end_sem"));}auto*g=new TGraphErrors(x.size(),x.data(),y.data(),ex.data(),ey.data());g->SetName((std::string("npe_end_")+mat).c_str());g->SetMarkerStyle(marker(mat));g->SetMarkerSize(1.35);g->SetMarkerColor(color(mat));g->SetLineColor(color(mat));g->SetLineWidth(3);mg.Add(g,"LP");leg.AddEntry(g,mat,"lp");f.cd();g->Write();}
  mg.Draw("A");mg.GetXaxis()->SetTitle("Muon position x [mm]");mg.GetYaxis()->SetTitle("Detected photoelectrons per END");mg.GetYaxis()->SetRangeUser(250,1410);leg.Draw();
  TLatex tx;tx.SetNDC();tx.SetTextSize(.035);tx.DrawLatex(.54,.85,"mean of left and right END yields");tx.DrawLatex(.54,.80,"error bars: event-level SEM");
  c.SaveAs("figures/npe_vs_x.pdf");f.cd();mg.Write("npe_vs_x");c.Write();provenance(f,"Event mean Npe/end=(Npe_left+Npe_right)/2; error is event-level SEM");f.Close();writeMeta("npe_vs_x","How much END light is detected versus position?","Mean Npe/end with event-level SEM, built and drawn in ROOT");
}

void makeDeltaTVsX() {
  style();auto t=readCsv("sources/timing_summary.csv");auto mt=readCsv("sources/material_summary.csv");TCanvas c("delta_t_vs_x","",920,640);TMultiGraph mg;TLegend leg(.16,.66,.43,.89);leg.SetBorderSize(0);leg.SetFillStyle(0);TFile f("figures/delta_t_vs_x.root","RECREATE");
  std::vector<TF1*> fits;
  for(auto mat:{"EJ-200","EJ-204","EJ-230"}){std::vector<double>x,y,ex,ey;for(auto&r:t.rows)if(r.at("material")==mat){x.push_back(d(r,"x_mm"));y.push_back(d(r,"mean_dt_ns"));ex.push_back(0);ey.push_back(d(r,"mean_dt_sem"));}auto*g=new TGraphErrors(x.size(),x.data(),y.data(),ex.data(),ey.data());g->SetName((std::string("mean_dt_")+mat).c_str());g->SetMarkerStyle(marker(mat));g->SetMarkerSize(1.25);g->SetMarkerColor(color(mat));g->SetLineColor(color(mat));mg.Add(g,"P");double a=0,b=0,v=0,ve=0;for(auto&r:mt.rows)if(r.at("material")==mat){a=d(r,"dt_intercept_ns");b=d(r,"dt_slope_ns_per_mm");v=d(r,"v_eff_mm_per_ns");ve=d(r,"v_eff_err");}auto*fit=new TF1((std::string("linear_")+mat).c_str(),"[0]+[1]*x",-670,670);fit->SetParameters(a,b);fit->SetLineColor(color(mat));fit->SetLineWidth(3);fits.push_back(fit);leg.AddEntry(g,(std::string(mat)+Form("  v_{eff}=%.2f#pm%.02f mm/ns",v,ve)).c_str(),"p");f.cd();g->Write();fit->Write();}
  mg.Draw("A");mg.GetXaxis()->SetTitle("Muon position x [mm]");mg.GetYaxis()->SetTitle("<#Deltat = t_{R}-t_{L}> [ns]");mg.GetYaxis()->SetRangeUser(-8.3,8.3);for(auto*fit:fits)fit->Draw("SAME");leg.Draw();TLatex tx;tx.SetNDC();tx.SetTextSize(.037);tx.DrawLatex(.58,.18,"#Deltat #approx -2x/v_{eff}");
  c.SaveAs("figures/delta_t_vs_x.pdf");f.cd();mg.Write("delta_t_vs_x");c.Write();provenance(f,"ROOT weighted linear fit to mean Delta t; point errors are SEM; v_eff=-2/slope");f.Close();writeMeta("delta_t_vs_x","How fast does timing information propagate longitudinally?","ROOT TGraphErrors and linear TF1 fit; v_eff=-2/slope");
}

void makePropagationTimes() {
  style();auto t=readCsv("sources/timing_summary.csv");TCanvas c("propagation_times","",920,640);TMultiGraph mg;TLegend leg(.16,.70,.38,.89);leg.SetBorderSize(0);leg.SetFillStyle(0);TFile f("figures/propagation_times.root","RECREATE");std::vector<double>x,l,r,ex,el,er;for(auto&z:t.rows)if(z.at("material")=="EJ-230"){x.push_back(d(z,"x_mm"));l.push_back(d(z,"mean_tL_ns"));r.push_back(d(z,"mean_tR_ns"));ex.push_back(0);el.push_back(d(z,"mean_tL_sem"));er.push_back(d(z,"mean_tR_sem"));}auto*gl=new TGraphErrors(x.size(),x.data(),l.data(),ex.data(),el.data());auto*gr=new TGraphErrors(x.size(),x.data(),r.data(),ex.data(),er.data());gl->SetName("mean_t_left");gr->SetName("mean_t_right");gl->SetMarkerStyle(20);gr->SetMarkerStyle(22);gl->SetMarkerColor(kBlue+1);gr->SetMarkerColor(kRed+1);gl->SetLineColor(kBlue+1);gr->SetLineColor(kRed+1);gl->SetLineWidth(3);gr->SetLineWidth(3);mg.Add(gl,"LP");mg.Add(gr,"LP");mg.Draw("A");mg.GetXaxis()->SetTitle("Muon position x [mm]");mg.GetYaxis()->SetTitle("Mean first-PE time [ns]");leg.AddEntry(gl,"left END","lp");leg.AddEntry(gr,"right END","lp");leg.Draw();TLatex tx;tx.SetNDC();tx.SetTextSize(.036);tx.DrawLatex(.56,.18,"EJ-230; eight SiPMs per END");c.SaveAs("figures/propagation_times.pdf");f.cd();gl->Write();gr->Write();mg.Write("propagation_times");c.Write();provenance(f,"Mean first detected PE time on each END; errors are event-level SEM");f.Close();writeMeta("propagation_times","How do left and right propagation times vary with x?","Mean first-PE time over the eight SiPMs at each END; event-level SEM");
}

void makeT0Fit(const std::string& cellName,const std::string& outName,const std::string& label) {
  style();TFile in("sources/timing_events.root","READ");auto*tr=(TTree*)in.Get("timing_events");char cell[32]={};double L[20],R[20];tr->SetBranchAddress("cell_id",cell);tr->SetBranchAddress("t_left_ns",L);tr->SetBranchAddress("t_right_ns",R);std::vector<double>v;for(Long64_t j=0;j<tr->GetEntries();j++){tr->GetEntry(j);if(cellName==cell)v.push_back(.5*(L[0]+R[0]));}if(v.size()!=10000)throw std::runtime_error("Expected 10000 events for "+cellName);std::sort(v.begin(),v.end());auto q=[&](double p){double z=p*(v.size()-1),a=std::floor(z),u=z-a;size_t n=a;return n+1<v.size()?v[n]*(1-u)+v[n+1]*u:v[n];};double mean=0;for(double z:v)mean+=z;mean/=v.size();double ss=0;for(double z:v)ss+=(z-mean)*(z-mean);double rms=std::sqrt(ss/v.size());double rob=.5*(q(.84)-q(.16));TH1D hall("t0_full","",190,v.front(),v.back());TH1D hz("t0_central","",180,q(.005),q(.995));for(double z:v){hall.Fill(z);hz.Fill(z);}TF1 first("initial_gaussian","gaus",q(.5)-2*rob,q(.5)+2*rob);first.SetParameters(hz.GetMaximum(),q(.5),rob);hz.Fit(&first,"QNR");double fm=first.GetParameter(1),fs=std::abs(first.GetParameter(2));TF1 fit("central_gaussian","gaus",fm-2*fs,fm+2*fs);fit.SetParameters(first.GetParameter(0),fm,fs);hz.Fit(&fit,"QRS");fit.SetLineColor(kRed+1);fit.SetLineWidth(3);
  TCanvas c(outName.c_str(),"",1120,570);c.Divide(2,1);c.cd(1);gPad->SetLogy();hall.SetLineColor(kBlue+1);hall.SetLineWidth(2);hall.GetXaxis()->SetTitle("T_{0} [ns]");hall.GetYaxis()->SetTitle("Events / bin");hall.Draw("HIST");TPaveText left(.15,.73,.61,.91,"NDC");left.SetFillColor(kWhite);left.SetFillStyle(1001);left.SetBorderSize(1);left.SetTextAlign(12);left.SetTextSize(.040);left.AddText((label+"; full range").c_str());left.AddText(Form("RMS = %.2f ps",1000*rms));left.Draw();c.cd(2);hz.SetMarkerStyle(20);hz.SetMarkerSize(.55);hz.GetXaxis()->SetTitle("T_{0} [ns]");hz.GetYaxis()->SetTitle("Events / bin");hz.Draw("E");fit.Draw("SAME");TPaveText right(.15,.50,.64,.91,"NDC");right.SetFillColor(kWhite);right.SetFillStyle(1001);right.SetBorderSize(1);right.SetTextAlign(12);right.SetTextSize(.034);right.AddText((label+"; central 99%").c_str());right.AddText(Form("#sigma_{G} = %.2f #pm %.2f ps",1000*fit.GetParameter(2),1000*fit.GetParError(2)));right.AddText(Form("#mu = %.4f #pm %.4f ns",fit.GetParameter(1),fit.GetParError(1)));right.AddText(Form("#chi^{2}/ndf = %.1f/%d = %.2f",fit.GetChisquare(),fit.GetNDF(),fit.GetChisquare()/fit.GetNDF()));right.AddText(Form("entries = %zu",v.size()));right.AddText("fit window: #mu #pm 2#sigma");right.Draw();c.SaveAs(("figures/"+outName+".pdf").c_str());TFile f(("figures/"+outName+".root").c_str(),"RECREATE");hall.Write();hz.Write();fit.Write();c.Write();provenance(f,"T0=(first END-left PE + first END-right PE)/2; ROOT Gaussian refit over mean +/- 2 sigma; full RMS retained");f.Close();writeMeta(outName,"What does the event-level T0 distribution look like?","Full distribution plus ROOT central Gaussian fit; RMS and fit width shown separately");
}
void makeT0FitCenter(){makeT0Fit("EJ230_xp0","t0_fit_ej230_x0","EJ-230, x = 0 mm");}
void makeT0FitEdge(){makeT0Fit("EJ230_xp650","t0_fit_ej230_x650","EJ-230, x = +650 mm");}

void makeSigmaT0VsX() {
  style();auto t=readCsv("sources/timing_summary.csv");TCanvas c("sigma_t0_vs_x","",920,640);TMultiGraph mg;TLegend leg(.16,.68,.39,.89);leg.SetBorderSize(0);leg.SetFillStyle(0);TFile f("figures/sigma_t0_vs_x.root","RECREATE");for(auto mat:{"EJ-200","EJ-204","EJ-230"}){std::vector<double>x,y,ex,ey;for(auto&r:t.rows)if(r.at("material")==mat){x.push_back(d(r,"x_mm"));y.push_back(1000*d(r,"t0_sigma_ns"));ex.push_back(0);ey.push_back(1000*d(r,"t0_sigma_err"));}auto*g=new TGraphErrors(x.size(),x.data(),y.data(),ex.data(),ey.data());g->SetName((std::string("sigma_t0_")+mat).c_str());g->SetMarkerStyle(marker(mat));g->SetMarkerSize(1.35);g->SetMarkerColor(color(mat));g->SetLineColor(color(mat));g->SetLineWidth(3);mg.Add(g,"LP");leg.AddEntry(g,mat,"lp");f.cd();g->Write();}mg.Draw("A");mg.GetXaxis()->SetTitle("Muon position x [mm]");mg.GetYaxis()->SetTitle("Central Gaussian #sigma(T_{0}) [ps]");mg.GetYaxis()->SetRangeUser(62,89);leg.Draw();TLatex tx;tx.SetNDC();tx.SetTextSize(.034);tx.DrawLatex(.49,.84,"T_{0}=(t_{L}^{(1)}+t_{R}^{(1)})/2");tx.DrawLatex(.49,.79,"error bars: ROOT fit error on #sigma");c.SaveAs("figures/sigma_t0_vs_x.pdf");f.cd();mg.Write("sigma_t0_vs_x");c.Write();provenance(f,"Central ROOT Gaussian sigma using first detected PE on each END; fit range mean +/- 2 sigma");f.Close();writeMeta("sigma_t0_vs_x","What optical timing resolution is reached versus position?","Central Gaussian sigma from ROOT fit; error is fit uncertainty on sigma");
}

void makeSigmaXVsX() {
  style();auto t=readCsv("sources/timing_summary.csv");TCanvas c("sigma_x_vs_x","",920,640);TMultiGraph mg;TLegend leg(.16,.68,.39,.89);leg.SetBorderSize(0);leg.SetFillStyle(0);TFile f("figures/sigma_x_vs_x.root","RECREATE");for(auto mat:{"EJ-200","EJ-204","EJ-230"}){std::vector<double>x,y,ex,ey;for(auto&r:t.rows)if(r.at("material")==mat){x.push_back(d(r,"x_mm"));y.push_back(d(r,"sigma_x_mm"));ex.push_back(0);ey.push_back(d(r,"sigma_x_err_mm"));}auto*g=new TGraphErrors(x.size(),x.data(),y.data(),ex.data(),ey.data());g->SetName((std::string("sigma_x_")+mat).c_str());g->SetMarkerStyle(marker(mat));g->SetMarkerSize(1.35);g->SetMarkerColor(color(mat));g->SetLineColor(color(mat));g->SetLineWidth(3);mg.Add(g,"LP");leg.AddEntry(g,mat,"lp");f.cd();g->Write();}mg.Draw("A");mg.GetXaxis()->SetTitle("Muon position x [mm]");mg.GetYaxis()->SetTitle("Longitudinal resolution #sigma_{x} [mm]");mg.GetYaxis()->SetRangeUser(10.5,16.0);leg.Draw();TLatex tx;tx.SetNDC();tx.SetTextSize(.037);tx.DrawLatex(.54,.19,"#sigma_{x} = (v_{eff}/2) #sigma(#Deltat)");c.SaveAs("figures/sigma_x_vs_x.pdf");f.cd();mg.Write("sigma_x_vs_x");c.Write();provenance(f,"sigma_x=v_eff*sigma_DeltaT/2; central ROOT Gaussian sigma_DeltaT; propagated fit errors");f.Close();writeMeta("sigma_x_vs_x","What longitudinal resolution follows from differential timing?","sigma_x=v_eff sigma_DeltaT/2 with propagated fit uncertainties");
}

void makeMaterialComparison() {
  style();auto t=readCsv("sources/timing_summary.csv");auto mt=readCsv("sources/material_summary.csv");double dec[3],att[3],npe[3],sig[3],ne[3],se[3];std::string mats[3]={"EJ-200","EJ-204","EJ-230"};for(int j=0;j<3;j++){for(auto&r:mt.rows)if(r.at("material")==mats[j]){dec[j]=d(r,"decay_ns");att[j]=d(r,"attenuation_m");}for(auto&r:t.rows)if(r.at("material")==mats[j]&&i(r,"x_mm")==0){npe[j]=d(r,"npe_end");ne[j]=d(r,"npe_end_sem");sig[j]=1000*d(r,"t0_sigma_ns");se[j]=1000*d(r,"t0_sigma_err");}}
  TCanvas c("material_comparison","",1050,730);c.Divide(2,2,0.01,0.01);TFile f("figures/material_comparison.root","RECREATE");const char*titles[4]={"Scintillation decay time [ns]","Bulk attenuation length [m]","Detected PE per END at x=0","Central #sigma(T_{0}) at x=0 [ps]"};double*arr[4]={dec,att,npe,sig};double*err[4]={nullptr,nullptr,ne,se};double ymax[4]={2.45,4.3,590,84};for(int p=0;p<4;p++){c.cd(p+1);gPad->SetLeftMargin(.17);TH1D*h=new TH1D(Form("panel_%d",p),titles[p],3,0,3);h->SetMinimum(0);h->SetMaximum(ymax[p]);h->GetYaxis()->SetTitle(titles[p]);for(int j=0;j<3;j++){h->GetXaxis()->SetBinLabel(j+1,mats[j].c_str());h->SetBinContent(j+1,arr[p][j]);if(err[p])h->SetBinError(j+1,err[p][j]);h->SetFillColor(j+1);}h->SetFillColor(kAzure-4);h->SetBarWidth(.62);h->SetBarOffset(.19);h->Draw(err[p]?"BAR E1":"BAR");for(int j=0;j<3;j++){TBox*b=new TBox(j+.19,0,j+.81,arr[p][j]);b->SetFillColor(color(mats[j]));b->SetLineColor(color(mats[j]));b->Draw();if(err[p]){TGraphErrors*g=new TGraphErrors(1);g->SetPoint(0,j+.5,arr[p][j]);g->SetPointError(0,0,err[p][j]);g->SetMarkerStyle(20);g->Draw("P");}}h->Draw("AXIS SAME");f.cd();h->Write();}
  c.SaveAs("figures/material_comparison.pdf");f.cd();c.Write();provenance(f,"Material constants from the configured material model; x=0 event mean Npe/end and ROOT central Gaussian timing width");f.Close();writeMeta("material_comparison","Why do the three scintillators give different timing?","Configured decay/attenuation compared with measured central light and timing");
}

void makeOrderScan() {
  style();auto t=readCsv("sources/order_scan_ej230_x0.csv");std::vector<double>x,y,ex,ey,ya,ea;for(auto&r:t.rows){x.push_back(d(r,"k"));y.push_back(1000*d(r,"sigma_t0_ns"));ey.push_back(1000*d(r,"sigma_t0_err"));ya.push_back(1000*d(r,"mean_first_m_sigma_ns"));ea.push_back(1000*d(r,"mean_first_m_sigma_err"));ex.push_back(0);}TCanvas c("order_scan","",920,640);TGraphErrors g(x.size(),x.data(),y.data(),ex.data(),ey.data());TGraphErrors ga(x.size(),x.data(),ya.data(),ex.data(),ea.data());g.SetName("kth_photon_central_gaussian_sigma");ga.SetName("mean_first_m_central_gaussian_sigma");g.SetMarkerStyle(20);g.SetMarkerColor(kRed+1);g.SetLineColor(kRed+1);g.SetLineWidth(3);ga.SetMarkerStyle(24);ga.SetMarkerColor(kBlue+1);ga.SetLineColor(kBlue+1);ga.SetLineStyle(2);ga.SetLineWidth(3);g.SetTitle("");g.GetXaxis()->SetTitle("Order k or number m of PE used at each END");g.GetYaxis()->SetTitle("Central Gaussian width [ps]");g.GetYaxis()->SetRangeUser(54,102);g.Draw("ALP");ga.Draw("LP SAME");TLegend leg(.17,.69,.53,.89);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.AddEntry(&g,"k-th PE: T_{0}(k)","lp");leg.AddEntry(&ga,"mean of first m PE per END","lp");leg.Draw();TLine mark(5,54,5,57.59);mark.SetLineColor(kGreen+2);mark.SetLineWidth(4);mark.Draw();TLatex tx;tx.SetNDC();tx.SetTextSize(.036);tx.SetTextColor(kGreen+2);tx.DrawLatex(.58,.85,"best tested: mean first 5");tx.DrawLatex(.58,.79,"57.59 #pm 0.66 ps");tx.SetTextColor(kBlack);tx.DrawLatex(.58,.71,"first PE: 67.14 #pm 0.77 ps");tx.DrawLatex(.58,.65,"gain: 9.55 ps (14.2%)");c.SaveAs("figures/order_scan.pdf");TFile f("figures/order_scan.root","RECREATE");g.Write();ga.Write();c.Write();provenance(f,"EJ-230 x=0; comparison of kth PE and arithmetic mean of first m PE on each END; central ROOT Gaussian fits");f.Close();writeMeta("order_scan","Is the first detected photon the best END estimator?","k=1..20 order statistic and mean-first-m scans at EJ-230 x=0; central ROOT Gaussian widths");
}

void makeNpeVsSigma() {
  style();auto t=readCsv("sources/timing_summary.csv");TCanvas c("npe_vs_sigma","",920,640);TLegend leg(.61,.68,.85,.89);leg.SetBorderSize(0);leg.SetFillStyle(0);TFile f("figures/npe_vs_sigma.root","RECREATE");TH2D frame("frame","",10,250,1400,10,62,89);frame.GetXaxis()->SetTitle("Detected photoelectrons per END");frame.GetYaxis()->SetTitle("Central #sigma(T_{0}) [ps]");frame.Draw();for(auto mat:{"EJ-200","EJ-204","EJ-230"}){std::vector<double>x,y,ex,ey;for(auto&r:t.rows)if(r.at("material")==mat){x.push_back(d(r,"npe_end"));y.push_back(1000*d(r,"t0_sigma_ns"));ex.push_back(d(r,"npe_end_sem"));ey.push_back(1000*d(r,"t0_sigma_err"));}auto*g=new TGraphErrors(x.size(),x.data(),y.data(),ex.data(),ey.data());g->SetName((std::string("npe_sigma_")+mat).c_str());g->SetMarkerStyle(marker(mat));g->SetMarkerSize(1.35);g->SetMarkerColor(color(mat));g->SetLineColor(color(mat));g->Draw("P SAME");leg.AddEntry(g,mat,"p");f.cd();g->Write();}leg.Draw();TLatex tx;tx.SetNDC();tx.SetTextSize(.037);tx.DrawLatex(.17,.87,"Near an end: more PE, but the far arm controls T_{0}");tx.DrawLatex(.17,.81,"Decay time and early paths break a universal 1/#sqrt{N_{pe}} law");c.SaveAs("figures/npe_vs_sigma.pdf");f.cd();frame.Write();c.Write();provenance(f,"All 21 cells; x coordinate event mean Npe/end, y coordinate ROOT central Gaussian timing width");f.Close();writeMeta("npe_vs_sigma","Does more detected light automatically improve timing?","All 21 current cells; measured Npe/end versus central Gaussian timing width");
}

void makeEndVsTop() {
  style();auto t=readCsv("sources/timing_summary.csv");double nend[3],ntop[3],send[3],stop[3],ne[3],nte[3],se[3],ste[3];std::string mats[3]={"EJ-200","EJ-204","EJ-230"};for(int j=0;j<3;j++)for(auto&r:t.rows)if(r.at("material")==mats[j]&&i(r,"x_mm")==0){nend[j]=d(r,"npe_end");ntop[j]=d(r,"npe_top");send[j]=1000*d(r,"t0_sigma_ns");stop[j]=1000*d(r,"top_sigma_ns");ne[j]=d(r,"npe_end_sem");nte[j]=d(r,"npe_top_sem");se[j]=1000*d(r,"t0_sigma_err");ste[j]=1000*d(r,"top_sigma_err");}TCanvas c("end_vs_top","",1050,600);c.Divide(2,1);TFile f("figures/end_vs_top.root","RECREATE");c.cd(1);gPad->SetLeftMargin(.17);TH1D hne("npe_end_center","",3,0,3),hnt("npe_top_center","",3,0,3);for(int j=0;j<3;j++){hne.GetXaxis()->SetBinLabel(j+1,mats[j].c_str());hne.SetBinContent(j+1,nend[j]);hne.SetBinError(j+1,ne[j]);hnt.SetBinContent(j+1,ntop[j]);hnt.SetBinError(j+1,nte[j]);}hne.SetTitle("");hne.SetMinimum(0);hne.SetMaximum(5600);hne.GetYaxis()->SetTitle("Detected PE at x=0");hne.SetMarkerStyle(20);hne.SetMarkerColor(kBlue+1);hne.SetLineColor(kBlue+1);hnt.SetMarkerStyle(22);hnt.SetMarkerColor(kGreen+2);hnt.SetLineColor(kGreen+2);hne.Draw("E1 P");hnt.Draw("E1 P SAME");TLegend l1(.18,.72,.52,.90);l1.SetBorderSize(0);l1.AddEntry(&hne,"per END (8 SiPMs)","p");l1.AddEntry(&hnt,"TOP total (70 SiPMs)","p");l1.Draw();c.cd(2);gPad->SetLeftMargin(.17);TH1D hse("sigma_end_center","",3,0,3),hst("sigma_top_center","",3,0,3);for(int j=0;j<3;j++){hse.GetXaxis()->SetBinLabel(j+1,mats[j].c_str());hse.SetBinContent(j+1,send[j]);hse.SetBinError(j+1,se[j]);hst.SetBinContent(j+1,stop[j]);hst.SetBinError(j+1,ste[j]);}hse.SetTitle("");hse.SetMinimum(0);hse.SetMaximum(90);hse.GetYaxis()->SetTitle("Central optical width [ps]");hse.SetMarkerStyle(20);hse.SetMarkerColor(kBlue+1);hse.SetLineColor(kBlue+1);hst.SetMarkerStyle(22);hst.SetMarkerColor(kGreen+2);hst.SetLineColor(kGreen+2);hse.Draw("E1 P");hst.Draw("E1 P SAME");TLegend l2(.18,.72,.76,.90);l2.SetBorderSize(0);l2.AddEntry(&hse,"END T_{0}: first PE on each END","p");l2.AddEntry(&hst,"TOP: first PE among 70 channels","p");l2.Draw();TLatex tx;tx.SetNDC();tx.SetTextSize(.032);tx.DrawLatex(.18,.64,"TOP is an ideal optical lower bound,");tx.DrawLatex(.18,.59,"not a validated readout estimator.");c.SaveAs("figures/end_vs_top.pdf");f.cd();hne.Write();hnt.Write();hse.Write();hst.Write();c.Write();provenance(f,"x=0 comparison; END T0 uses two eight-SiPM arms; TOP time is global first PE among 70 sensors");f.Close();writeMeta("end_vs_top","What information is available from END and TOP readout?","x=0: END T0 versus ideal global first TOP PE; estimators are intentionally distinguished");
}

void makeFitGridEJ230() {
  style();TFile in("sources/t0_fit_diagnostics.root","READ");std::vector<int> xs={-650,-500,-200,0,200,500,650};
  TCanvas c("fit_grid_ej230_v9","",1450,760);c.Divide(4,2,.002,.002);TFile out("figures/fit_grid_ej230_v9.root","RECREATE");
  for(size_t j=0;j<xs.size();j++){c.cd(j+1);gPad->SetLeftMargin(.15);gPad->SetBottomMargin(.15);std::string id=xs[j]<0?"EJ230_xm"+std::to_string(-xs[j]):"EJ230_xp"+std::to_string(xs[j]);auto*h=(TH1D*)in.Get((id+"/t0_histogram").c_str());auto*f=(TF1*)in.Get((id+"/gaussian_fit").c_str());if(!h||!f)throw std::runtime_error("Missing "+id);auto*hc=(TH1D*)h->Clone(("hist_"+id).c_str());auto*fc=(TF1*)f->Clone(("fit_"+id).c_str());hc->SetDirectory(&out);hc->SetMarkerStyle(20);hc->SetMarkerSize(.35);hc->GetXaxis()->SetTitle("T_{0} [ns]");hc->GetYaxis()->SetTitle("events / bin");fc->SetLineColor(kRed+1);fc->SetLineWidth(3);hc->Draw("E");fc->Draw("SAME");auto*box=new TPaveText(.14,.69,.72,.91,"NDC");box->SetFillColor(kWhite);box->SetBorderSize(1);box->SetTextAlign(12);box->SetTextSize(.043);box->AddText(Form("x = %+d mm",xs[j]));box->AddText(Form("#sigma_{G}=%.2f#pm%.2f ps",1000*fc->GetParameter(2),1000*fc->GetParError(2)));box->AddText(Form("#chi^{2}/ndf=%.2f",fc->GetChisquare()/fc->GetNDF()));box->Draw();out.cd();hc->Write();fc->Write();}
  c.cd(8);gPad->Range(0,0,1,1);TLatex tx;tx.SetTextFont(42);tx.SetTextSize(.065);tx.DrawLatex(.10,.78,"EJ-230: all seven positions");tx.SetTextSize(.052);tx.DrawLatex(.10,.62,"Exact central histograms used");tx.DrawLatex(.10,.53,"for the quoted resolution.");tx.DrawLatex(.10,.37,"Common two-pass fit policy:");tx.DrawLatex(.10,.28,"final window = #mu #pm 2#sigma.");
  c.SaveAs("figures/fit_grid_ej230_v9.pdf");out.cd();c.Write();provenance(out,"Seven exact T0 histograms and stored ROOT TF1 fits used in timing_summary for EJ-230");out.Close();writeMeta("fit_grid_ej230_v9","Are the seven EJ-230 Gaussian fits stable?","Exact stored T0 histograms and ROOT TF1 Gaussian fits, one panel per x");
}

void makeMirrorOverlaysEJ230() {
  style();TFile in("sources/timing_events.root","READ");auto*tr=(TTree*)in.Get("timing_events");char cell[32]={};double l[20],r[20];tr->SetBranchAddress("cell_id",cell);tr->SetBranchAddress("t_left_ns",l);tr->SetBranchAddress("t_right_ns",r);std::map<int,std::vector<double>> vals;for(Long64_t j=0;j<tr->GetEntries();j++){tr->GetEntry(j);std::string id=cell;if(id.rfind("EJ230_",0)!=0)continue;int x=id.find("xm")!=std::string::npos?-std::stoi(id.substr(id.find("xm")+2)):std::stoi(id.substr(id.find("xp")+2));if(x!=0)vals[x].push_back(.5*(l[0]+r[0]));}
  TCanvas c("mirror_overlays_ej230","",1400,460);c.Divide(3,1,.005,.005);TFile out("figures/mirror_overlays_ej230.root","RECREATE");int pad=0;for(int a:{200,500,650}){c.cd(++pad);gPad->SetLeftMargin(.16);gPad->SetBottomMargin(.15);double mn=std::accumulate(vals[-a].begin(),vals[-a].end(),0.0)/vals[-a].size(),mp=std::accumulate(vals[a].begin(),vals[a].end(),0.0)/vals[a].size();auto*hn=new TH1D(Form("centered_xm%d",a),"",120,-.30,.30);auto*hp=new TH1D(Form("centered_xp%d",a),"",120,-.30,.30);hn->SetDirectory(&out);hp->SetDirectory(&out);for(double z:vals[-a])hn->Fill(z-mn);for(double z:vals[a])hp->Fill(z-mp);hn->Scale(1.0/hn->Integral("width"));hp->Scale(1.0/hp->Integral("width"));hn->SetLineColor(kBlue+1);hp->SetLineColor(kRed+1);hn->SetLineWidth(3);hp->SetLineWidth(3);hn->GetXaxis()->SetTitle("T_{0}-<T_{0}> [ns]");hn->GetYaxis()->SetTitle("normalized density");hn->SetMaximum(1.13*std::max(hn->GetMaximum(),hp->GetMaximum()));hn->Draw("HIST");hp->Draw("HIST SAME");auto*leg=new TLegend(.61,.72,.90,.90);leg->SetBorderSize(0);leg->SetFillStyle(0);leg->AddEntry(hn,Form("x=-%d mm",a),"l");leg->AddEntry(hp,Form("x=+%d mm",a),"l");leg->Draw();auto*tx=new TLatex();tx->SetNDC();tx->SetTextSize(.043);tx->DrawLatex(.18,.87,Form("|x| = %d mm",a));out.cd();hn->Write();hp->Write();}
  c.SaveAs("figures/mirror_overlays_ej230.pdf");out.cd();c.Write();provenance(out,"EJ-230 mirrored pairs; each T0 distribution centered at its own mean; common range [-0.30,0.30] ns and 120 bins; unit area");out.Close();writeMeta("mirror_overlays_ej230","Do mirrored EJ-230 T0 shapes agree after removing their means?","Normalized centered overlays; identical range and binning for all mirrored pairs");
}

void makeWidthEstimatorsEJ230() {
  style();auto t=readCsv("sources/t0_fit_diagnostics.csv");std::vector<double>x,ex,sg,sge,rms,rmse,q,qe;for(auto&r:t.rows)if(r.at("material")=="EJ-230"){x.push_back(d(r,"x_mm"));ex.push_back(0);sg.push_back(1000*d(r,"sigma_G_ns"));sge.push_back(1000*d(r,"sigma_G_err_ns"));rms.push_back(1000*d(r,"global_rms_ns"));rmse.push_back(1000*d(r,"global_rms_bootstrap_err_ns"));q.push_back(1000*d(r,"qwidth_ns"));qe.push_back(1000*d(r,"qwidth_bootstrap_err_ns"));}TCanvas c("width_estimators_ej230","",930,650);TGraphErrors gs(x.size(),x.data(),sg.data(),ex.data(),sge.data()),gr(x.size(),x.data(),rms.data(),ex.data(),rmse.data()),gq(x.size(),x.data(),q.data(),ex.data(),qe.data());gs.SetName("gaussian_sigma");gr.SetName("global_rms");gq.SetName("robust_qwidth");gs.SetMarkerStyle(20);gs.SetMarkerColor(kRed+1);gs.SetLineColor(kRed+1);gr.SetMarkerStyle(21);gr.SetMarkerColor(kBlue+1);gr.SetLineColor(kBlue+1);gq.SetMarkerStyle(22);gq.SetMarkerColor(kGreen+2);gq.SetLineColor(kGreen+2);for(auto*g:{&gs,&gr,&gq}){g->SetLineWidth(3);g->SetMarkerSize(1.2);}gs.SetTitle("");gs.GetXaxis()->SetTitle("Muon position x [mm]");gs.GetYaxis()->SetTitle("Width of T_{0} [ps]");gs.GetYaxis()->SetRangeUser(63,90);gs.Draw("ALP");gr.Draw("LP SAME");gq.Draw("LP SAME");TLegend leg(.16,.68,.43,.89);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.AddEntry(&gs,"central Gaussian #sigma","lp");leg.AddEntry(&gr,"global RMS","lp");leg.AddEntry(&gq,"(q_{84}-q_{16})/2","lp");leg.Draw();TLatex tx;tx.SetNDC();tx.SetTextSize(.031);tx.DrawLatex(.47,.18,"RMS/q-width errors: 300 event bootstraps");c.SaveAs("figures/width_estimators_ej230.pdf");TFile out("figures/width_estimators_ej230.root","RECREATE");gs.Write();gr.Write();gq.Write();c.Write();provenance(out,"EJ-230 comparison of central Gaussian sigma, global RMS, and robust half-width; bootstrap seed 26091443 series");out.Close();writeMeta("width_estimators_ej230","Does the apparent x asymmetry depend on the width definition?","EJ-230 Gaussian sigma, global RMS and robust q-width; 300-replica ROOT bootstrap errors for RMS/q-width");
}
}
#endif
