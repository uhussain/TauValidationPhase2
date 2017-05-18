#include "tdrstyle.C"
#include "PlotStyles.C"
#include "CMS_lumi.C"
#include "plotTimeReturnTGraphAsymmErrors.C"
#include "plotTimeReturnTGraphAsymmErrorsRange3.C"
#include "plotTimeReturnTGraphAsymmErrorsRange4.C"
#include "plotTimeReturnTGraphAsymmErrors1Line.C"

void plotIsoEffAllSamples_DM(){
  TString isoCut = "2";
  Int_t color1 = TColor::GetColor("#283593"); //dark blue
  Int_t color2 = TColor::GetColor("#0288D1"); //medium blue

  Int_t color3 = TColor::GetColor("#00695C"); //green blue
  Int_t color4 = TColor::GetColor("#F44336"); //red

  Int_t color5 = TColor::GetColor("#FF9900"); //Mustard Yellow

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
  
  TCanvas *c1 = new TCanvas("c","c",800,800);
  c1->SetGrid();
  setCanvasStyle(c1);
  c1->cd();

  TString plotName = "compilation-tauEfficiency-ZTT-DM0_May18";
  TString legLabel = "";
  TH1F *basehist = new TH1F("basehist","",100,0,2.5);

  TChain pu140_gaus("PFChargedBased/Ntuple");  
  pu140_gaus.Add("/data/uhussain/TauTiming_May18_hadd/Ztt_pre4_pu140.root");

  TChain pu200_gaus("PFChargedBased/Ntuple");  
  pu200_gaus.Add("/data/uhussain/TauTiming_May18_hadd/Ztt_pre4_pu200.root");

  TChain pu0_gaus("PFChargedBased/Ntuple");  
  pu0_gaus.Add("/data/uhussain/TauTiming_May18_hadd/Ztt_pre4_pu0.root");

  basehist->SetStats(false);
  TString iso200("0.050"), iso140("0.050");
  TString date = "2_2_17";


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

  TString numerator = "genTauPt > 22 && abs(genTauEta) < 2.3 && (dmf!=5&&dmf!=6 && dmf > -1) && (dmf == 0)  &&  leadPFChHadrCandPt> 22 && abs(leadPFChHadrCandEta)<2.3";
  TString denominator = "genTauPt > 22 && abs(genTauEta) <2.3 && (dmf!=5&&dmf!=6 && dmf > -1) && (dmf == 0)";
  
  TString numerator1 = "genTauPt > 22 && abs(genTauEta)>=2.3 && abs(genTauEta) < 4.0 && (dmf!=5&&dmf!=6 && dmf > -1) && (dmf == 0)  &&  leadPFChHadrCandPt> 22 && abs(leadPFChHadrCandEta)>=2.3 && abs(leadPFChHadrCandEta < 4.0)";
  TString denominator1 = "genTauPt > 22 && abs(genTauEta)>=2.3 && abs(genTauEta) <4.0 && (dmf!=5&&dmf!=6 && dmf > -1) && (dmf == 0)" ;
  TString logand = " && ";

  TString numerator_dm0 = "genTauPt > 22 && abs(genTauEta) < 4.0 && (dmf!=5&&dmf!=6 && dmf > -1) && (dmf == 0 || dmf == 1 || dmf == 10)  &&  leadPFChHadrCandPt> 22 && abs(leadPFChHadrCandEta)<4.0";
  TString denominator_dm0 = "genTauPt > 22 && abs(genTauEta) <4.0 && (dmf!=5&&dmf!=6 && dmf > -1) && (dmf == 0 || dmf == 1 || dmf == 10)";
  
  TString numeratorNominal = numerator + "&& PFCharged <"+isoCut;
  //TString numeratorNominal_pu0 = numerator_pu0 + "&& PFCharged <"+isoCut;

  TGraphAsymmErrors *pu200_eff = plotTimeReturnTGraphAsymmErrorsRange3(pu200_gaus, numerator, denominator);
  TGraphAsymmErrors *pu140_eff = plotTimeReturnTGraphAsymmErrorsRange4(pu140_gaus, numerator, denominator);
  TGraphAsymmErrors *pu0_eff   = plotTimeReturnTGraphAsymmErrors1Line(pu0_gaus, numerator,denominator);

  /////////////////////////////////////////
  //TString numeratorT4 = numerator + "&& PFChargedT4 <"+isoCut;
  //TString numeratorT4v0 = numerator + "&& PFChargedT4 <5";

  TString numeratorNominal1 = numerator1 + "&& PFCharged <"+isoCut;
  TGraphAsymmErrors *pu200_eff_etaExtension = plotTimeReturnTGraphAsymmErrorsRange3(pu200_gaus, numerator1, denominator1);  
  TGraphAsymmErrors *pu140_eff_etaExtension = plotTimeReturnTGraphAsymmErrorsRange4(pu140_gaus, numerator1, denominator1);  
  TGraphAsymmErrors *pu0_eff_etaExtension   = plotTimeReturnTGraphAsymmErrors(pu0_gaus, numerator1, denominator1);

  /////////////////////////////////////////
  basehist->GetXaxis()->SetTitle("Density (events / mm)");
  basehist->GetYaxis()->SetTitle("Tau Efficiency ");  
  basehist->GetYaxis()->SetRangeUser(0.0,1.0);
  basehist->GetXaxis()->SetRangeUser(0.0,2);
  basehist->Draw("");

  //setPlotStyleAsymm(                  plot , Int_t  color, Int_t fillStyle,  Int_t MarkerStyle){

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

  pu0_eff->SetFillColor(kBlue+4);
  //pu0_eff->SetMarkerStyle(24);
  pu0_eff->SetFillStyle(3005);
  pu0_eff->Draw("4, Same");
  //pu0_eff->SetFillColor(kBlue+4);
  //pu0_eff->SetFillStyle(3005);
  //pu0_eff->Draw("4, Same");


  TLegend *leg = new TLegend(.55, .706, .75, .92,legLabel,"nbNDC");  
  setLegendStyles(leg,legLabel, 4);

 // TLegend *leg = new TLegend(.45, .15, .65, .32,legLabel,"nbNDC");  
 // leg->SetBorderSize(0);
 // leg->SetTextFont(42);
 // leg->SetShadowColor(kWhite);
 // //leg->SetFillColor(kWhite);
 // leg->SetTextSize(0.025);
  leg->SetHeader("Z #rightarrow #tau #tau DecayMode = 0");
  leg->AddEntry(            pu0_eff,   " Zero PU","F");  
  leg->AddEntry(pu140_eff,"140PU","PL");
  leg->AddEntry(pu140_eff_etaExtension,"140PU,2.3<|#eta|<4.0","PL");

  leg->AddEntry(pu200_eff,"200PU","PL");
  leg->AddEntry(pu200_eff_etaExtension,"200PU, 2.3<|#eta|<4.0","PL");

  //leg->AddEntry(pu0_eff_etaExtension,"Zero Pileup + 2.3<|#eta|<4.0","PL");
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
  c1->SaveAs("May18Plots/"+plotName+".png");
}
