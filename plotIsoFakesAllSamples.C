#include "tdrstyle.C"
#include "PlotStyles.C"
#include "CMS_lumi.C"
#include "plotTimeReturnTGraphAsymmErrors.C"
#include "plotTimeReturnTGraphAsymmErrorsRange2.C"
#include "plotTimeReturnTGraphAsymmErrors1Line.C"
#include "plotDistributions.C"

void plotIsoFakesAllSamples(){

  TString isoCut = "2";

  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "       Phase-2 Simulation Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = "14 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 0; //0 is out of the frame

  gStyle->SetOptStat(0);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(4);

  TCanvas *c1 = new TCanvas("c","c",800,800);
  setCanvasStyle(c1);
  c1->cd();

  // pu 200
  TChain pu200_gaus("PFChargedBased/jetNtuple");  
  pu200_gaus.Add("/data/uhussain/TauTiming_Feb20_hadd/Ztt_pre4_pu200.root");
  
  //pu 140
  TChain pu140_gaus("PFChargedBased/jetNtuple");  
  pu140_gaus.Add("/data/uhussain/TauTiming_Feb20_hadd/Ztt_pre4_pu140.root");
  
  TChain pu0_gaus("PFChargedBased/jetNtuple");  
  pu0_gaus.Add("/data/uhussain/TauTiming_Feb20_hadd/Ztt_pre4_pu0.root");
  
  TString plotName = "compilation-jetFakeProbability-ZTT-denom-genJetMatch-reco-V2";
  TString legLabel = "";

  TH1F *basehist = new TH1F("basehist","",100,0,5);
  basehist->SetStats(false);
  TString iso200("0.050"), iso140("0.050");
  TString date = "2_2_17";

  //TString numerator = "jetPt > 22 && jetPt < 400 && genJetMatch > 0 &&  dmf==10 && good3ProngT3 > 0  && abs(jetEta) <2.1 && abs(tauEta)<2.1 && tauPt> 30 && vtxIndex==0";
  //TString denominator = "jetPt > 22 && jetPt < 400 && genJetMatch > 0 && abs(jetEta) <2.1 && vtxIndex==0";
  TString numerator = "jetPt > 22 && jetPt < 400 && genJetMatch>0 && (dmf!=5&&dmf!=6 && dmf > 0) && abs(jetEta) <2.3 && abs(tauEta)<2.3 && tauPt> 27";
  TString denominator = "jetPt > 22 && jetPt < 400 && genJetMatch>0 && abs(jetEta) <2.3 && (dmf!=5&&dmf!=6 && dmf > 0)";

  TString numerator1 = "jetPt > 22 && jetPt < 400 && genJetMatch>0 && (dmf!=5&&dmf!=6 && dmf > 0) && abs(jetEta)>=2.3 && abs(jetEta) <4.0 && abs(tauEta)>=2.3 && abs(tauEta)<4.0 && tauPt> 27";
  TString denominator1 = "jetPt > 22 && jetPt < 400 && genJetMatch>0 && abs(jetEta)>=2.3 && abs(jetEta) <4.0  && (dmf!=5&&dmf!=6 && dmf > 0)";
  
  TString logand = " && ";

  double zs_200[4]   = {0.,0.,0.,0.};
  double vals_200[4] = {0.,0.,0.,0.};
  double erhi_200[4] = {0.,0.,0.,0.};
  double erlo_200[4] = {0.,0.,0.,0.};
  double zs_140[4]   = {0.,0.,0.,0.};
  double vals_140[4] = {0.,0.,0.,0.};
  double erhi_140[4] = {0.,0.,0.,0.};
  double erlo_140[4] = {0.,0.,0.,0.};
  double dens_200[4] = {0.,0.,0.,0.};

  double zs_time_200[4]   = {0.,0.,0.,0.};
  double vals_time_200[4] = {0.,0.,0.,0.};
  double erhi_time_200[4] = {0.,0.,0.,0.};
  double erlo_time_200[4] = {0.,0.,0.,0.};
  double zs_time_140[4]   = {0.,0.,0.,0.};
  double vals_time_140[4] = {0.,0.,0.,0.};
  double erhi_time_140[4] = {0.,0.,0.,0.};
  double erlo_time_140[4] = {0.,0.,0.,0.};

  double exl[4] = {0.,0.,0.,0.};
  double exh[4] = {0.,0.,0.,0.};

  //plotDistributions(pu0_gaus,denominator,"jetPt-denominator-0PU");

  TString numeratorNominal = numerator + "&& PFCharged <" + isoCut;
  //plotDistributions(pu200_gaus,numeratorNominal,"jetPt-nominal");
  TGraphAsymmErrors *pu200_eff = plotTimeReturnTGraphAsymmErrorsRange2(pu200_gaus, numerator,denominator);
  TGraphAsymmErrors *pu140_eff = plotTimeReturnTGraphAsymmErrors(pu140_gaus, numerator,denominator);
  TGraphAsymmErrors *pu0_eff   = plotTimeReturnTGraphAsymmErrors1Line(pu0_gaus, numerator,denominator);

  /////////////////////////////////////////

  TString numeratorNominal1 = numerator1 + "&& PFCharged <"+isoCut;
  //plotDistributions(pu200_gaus,numeratorT4,"jetPt-T4");
  TGraphAsymmErrors *pu200_eff_etaExtension = plotTimeReturnTGraphAsymmErrorsRange2(pu200_gaus, numerator1,denominator1);
  TGraphAsymmErrors *pu140_eff_etaExtension = plotTimeReturnTGraphAsymmErrors(pu140_gaus, numerator1,denominator1);
  TGraphAsymmErrors *pu0_eff_etaExtension   = plotTimeReturnTGraphAsymmErrors(pu0_gaus,   numerator1,denominator1);

  /////////////////////////////////////////
  basehist->GetXaxis()->SetTitle("Density (events / mm)");
  basehist->GetYaxis()->SetTitle("Tau Misidentification Probability ");  
  basehist->GetYaxis()->SetRangeUser(0.0,0.35);
  basehist->GetXaxis()->SetRangeUser(0.0,2);
  basehist->GetYaxis()->SetLabelSize(0.030);

  basehist->Draw("");
  

  setPlotStyleAsymm(  pu200_eff_etaExtension,       color4,            3005,                 33);
  pu200_eff_etaExtension->Draw("P Same");

  setPlotStyleAsymm(  pu140_eff_etaExtension,       color1,            3005,                 27);
  pu140_eff_etaExtension->Draw("P Same");

  //setPlotStyleAsymm(  pu0_eff_etaExtension,       color2,            3005,                 33);
  //pu0_eff_etaExtension->Draw("P Same");

  //setPlotStyleAsymm(                  plot , Int_t  color, Int_t fillStyle,  Int_t MarkerStyle){
  setPlotStyleAsymm(             pu140_eff,       color3,            3005,                 32);
  pu140_eff->Draw("P Same");

//setPlotStyleAsymm(                  plot , Int_t  color, Int_t fillStyle,  Int_t MarkerStyle){
  setPlotStyleAsymm(             pu200_eff,       color5,            3005,                 23);
  pu200_eff->Draw("P Same");
  ////////////////// PU 0
//setPlotStyleAsymm(                  plot , Int_t  color, Int_t fillStyle,  Int_t MarkerStyle){
  //TGraphErrors *errorBand = new TGraphErrors(n,x,y,ex,ey);
  //pu0_eff_etaExtension->SetFillColor(kBlue+4);
  //pu0_eff_etaExtension->SetFillStyle(3005);   
  //setPlotStyleAsymm(  pu0_eff_etaExtension,       color3,            3005,                 23);
  pu0_eff->SetFillColor(kBlue+4); 
  //pu0_eff->SetMarkerStyle(24);
  pu0_eff->SetFillStyle(3005);
  pu0_eff->Draw("4, Same");


  //setLegendStyles options
  //Legend option 0 == manual set
  //              1 == upper left
  //              2 == lower left
  //              3 == lower right
  //              4 == upper right
  //              5 == center upper
  //              6 == center lower
  //
  TLegend *leg = new TLegend(.15, .606, .35, .92,legLabel,"nbNDC");
  setLegendStyles(leg,legLabel, 2);

  leg->SetHeader("Z #rightarrow #tau #tau");
  leg->AddEntry(            pu0_eff,   " Zero PU","F");  
  leg->AddEntry(pu140_eff,"140PU","PL");
  leg->AddEntry(pu140_eff_etaExtension,"140PU,2.3<|#eta|<4.0","PL");

  leg->AddEntry(pu200_eff,"200PU","PL");
  leg->AddEntry(pu200_eff_etaExtension,"200PU, 2.3<|#eta|<4.0","PL");



  leg->Draw();

   c1->cd();
  //Writing the lumi info
  CMS_lumi(c1, iPeriod, iPos );

  c1->Update();
  c1->RedrawAxis();
  c1->GetFrame()->Draw();
  //TLine *line = new TLine(0,0.15,3,0.15);
  //line->Draw();  
  c1->SaveAs("May18Plots/"+plotName+".pdf");
  //c1->SaveAs("May18/+"plotName+".png");
}
