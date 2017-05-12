
//#include "TCanvas.h"
//#include "TROOT.h"
//#include "TFile.h"
//#include "TTree.h"
//#include "TChain.h"
//#include "TH1.h"
//#include "TF1.h"
//#include <math.h>
//#include <TEfficiency.h>
//#include <TMath.h>
//#include <TFormula.h>
//#include <iostream>
//#include <string>
//#include <iostream>
//#include <cmath>
//#include "TLegend.h"
//#include "TPad.h"
//#include "TH2.h"
//#include "TF1.h"
//#include "THStack.h"
//#include "TStyle.h"
//#include "TAxis.h"
//#include "TGaxis.h"
//#include "TTree.h"
//#include "TPaveText.h"
#include "tdrstyle.C"
//#include "TStyle.h"
//#include "TAxis.h"
//#include "TGaxis.h"
#include "CMS_lumi.C"
#include "TH1.h"
#include "TH1F.h"
void applyPadStyle(TPad* pad1){
  pad1->SetFillColor(0);
  pad1->Draw();  pad1->cd();  pad1->SetLeftMargin(0.15);  pad1->SetBottomMargin(0.13); pad1->SetRightMargin(0.05);
  //pad1->SetGrid();
  pad1->SetGrid(10,10);
}

void plotEfficiency2(){

  setTDRStyle();
  TString fileName = "test_ntuple_ztt_pu0_2023D13_910_pre3_Phase2eta_v1.root";
  TString fileName1 = "test_ntuple_ztt_pu140_2023D13_910_pre3_Phase2eta_v1.root";
  TString fileName2 = "test_ntuple_ztt_pu200_2023D13_910_pre3_Phase2eta_v1.root";
  TString fileName3 = "test_ntuple_VBF_Phase1_pu30_RuII_13TeV_820_patch1_v1.root";
  TString treePath = "Events";
  int bins = 20;
  double low  = 0;
  double high = 120;
  writeExtraText = true;       // if extra text
  extraText  = "#splitline{9_1_0_pre3 Validation}{RelValZTT 2023D13}";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.1 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_sqrtS = "14 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 11; //0 is out of the frame

  TString isoCut = "2.5";
  TString dmf = "0.5";
  //Plotting Variables
  TString variable = "gen_tauvis_pt";
  TString variable1 = "gen_tauvis_pt_reco";
  TString GenCut= "gen_tauvis_pt>22&&abs(gen_tauvis_eta)<2.4";
  //TString GenCut= "gen_tauvis_pt>22&&abs(gen_tauvis_eta)<4.0";
  //&& (tau_dm==0 || tau_dm==1 || tau_dm==10)";//denominator
  //TString GenCut1 = "gen_tauvis_pt>22&&abs(gen_tauvis_eta)>=2.1&&abs(gen_tauvis_eta)<4.0";
  //TString RecoCut= "gen_tauvis_pt_reco>20&&" + GenCut;//numerator
  //TString RecoCut1= "gen_tauvis_pt_reco>22&&" + GenCut1;
  //TString RecoCut2= "gen_tauvis_pt_reco>22&&tau_chargedIsoPtSum<"+isoCut+"&&" + GenCut;
//  TString etaCut1 =  "((tau_photonIsoPtSum/tau_pt) <0.50 && abs(gen_tauvis_eta)<1.5)";
//  TString etaCut2 =  "((tau_photonIsoPtSum/tau_pt) <0.55 && (abs(gen_tauvis_eta)>1.5 && abs(gen_tauvis_eta)<2.4))";
//  TString etaCut3 =  "((tau_photonIsoPtSum/tau_pt) <0.60 && (abs(gen_tauvis_eta)>2.4 && abs(gen_tauvis_eta)<4.0))";
//  TString etaCutFinal = etaCut1+"||"+ etaCut2 + "||" + etaCut3;
  //TString RecoCut3= "gen_tauvis_pt_reco>22&&abs(gen_tauvis_eta_reco)<4.0&& (tau_dm_match!=5 && tau_dm_match!=6 && tau_dm_match > -1)&&"+ GenCut; 
  TString RecoCut3= "tau_pt_match>22&&abs(tau_eta_match)<2.4  && tau_dm_match > -1";
    //&& (tau_dm_match!=5 && tau_dm_match!=6 && tau_dm_match > -1)";
  //TString RecoCut3= "tau_dm==1&&tau_dmf>"+dmf+"&&" + GenCut;
  //TString RecoCut3 = RecoCut + "&&"+ etaCutFinal + "&&tau_DBIso< 0.8";
  //TString RecoCut3= "tau_dm==1&&tau_dmf>"+dmf+"&&" + GenCut;
  //Style
  TString xaxis = "Gen #tau_{vis} pt";
  int markerstyle = 20;
  int markerstyle1 = 21; 
  int markerstyle2 = 22;
  int markerstyle3 = 24;
  Color_t color = TColor::GetColor("#283593");//dark blue color1
  Color_t color1 = TColor::GetColor("#F44336");//red color4
  //Color_t color3 = TColor::GetColor("#00695C"); //green blue color2
  Color_t color2 = TColor::GetColor("#0288D1"); //medium blue 
  Color_t color3 = TColor::GetColor(1); //black
  TString outFileName = "plotEfficiency-VBF_eta_dmf_new_allpu";

  TString legLabel = "";

//  setTDRStyle();
  //tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");

  TFile *tauFile    = new TFile(fileName);

  if(!tauFile->IsOpen()||tauFile==0){
    std::cout<<"ERROR FILE "<< fileName<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }

  TFile *tauFile1    = new TFile(fileName1);

  if(!tauFile1->IsOpen()||tauFile1==0){
    std::cout<<"ERROR FILE "<< fileName1<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }

  TFile *tauFile2    = new TFile(fileName2);

  if(!tauFile2->IsOpen()||tauFile2==0){
    std::cout<<"ERROR FILE "<< fileName2<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }

  TFile *tauFile3    = new TFile(fileName3);

  if(!tauFile3->IsOpen()||tauFile3==0){
    std::cout<<"ERROR FILE "<< fileName3<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  int W = 800;
  int H = 800;
  
  int H_ref = 800; 
  int W_ref = 800; 
  //references for T,B,L,R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref; 
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;
  TString canvName = "RelValZTT_tdr_PhaseII_allpu_pt_9_1_0_pre3_dm";
  canvName += W;
  canvName += "-";
  canvName += H;
  canvName += "_";  
  canvName += iPeriod;
  if( writeExtraText ) canvName += "-prelim";
  if( iPos%10==0 ) canvName += "-out";
  else if( iPos%10==1 ) canvName += "-left";
  else if( iPos%10==2 )  canvName += "-center";
  else if( iPos%10==3 )  canvName += "-right";

  TCanvas* canv = new TCanvas(canvName,canvName,50,50,W,H);
  canv->SetGrid();
  canv->SetFillColor(0);
  canv->SetBorderMode(0);
  canv->SetFrameFillStyle(0);
  canv->SetFrameBorderMode(0);
  canv->SetLeftMargin( L/W );
  canv->SetRightMargin( R/W );
  canv->SetTopMargin( T/H );
  canv->SetBottomMargin( B/H );
  canv->SetTickx(0);
  canv->SetTicky(0);
//  TCanvas *Tcan= new TCanvas("Tcan","",100,20,600,600);
//  
//  Tcan->cd();  Tcan->SetFillColor(0);
//  TPad* pad1 = new TPad("pad1","The pad",0,0,0.98,0.98);

  //applyPadStyle(pad1);
//  gStyle->SetOptFit(0);
//  gStyle->SetOptStat(0);
//  gROOT->ForceStyle();
  gStyle->SetEndErrorSize(0);

  TTree* tauTree = (TTree*)tauFile->Get(treePath);
  if(tauTree == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }

  TTree* tauTree1 = (TTree*)tauFile1->Get(treePath);
  if(tauTree1 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree1<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }

  TTree* tauTree2 = (TTree*)tauFile2->Get(treePath);
  if(tauTree2 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree2<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }

  TTree* tauTree3 = (TTree*)tauFile3->Get(treePath);
  if(tauTree3 == 0){
    std::cout<<"ERROR Tau Tree is "<< tauTree3<<" NOT FOUND; EXITING"<<std::endl;
    exit(0);
  }
  /// first
  TH1* Denom;
  Denom = new TH1F("Denom","Denom",bins,low,high);
  Denom->Sumw2();
  tauTree->Draw(variable+">>+Denom",GenCut);
  std::cout<<"Denom_pu0: "<<Denom->GetEntries()<<std::endl;
  TH1* Num;
  Num = new TH1F("Num","Num",bins,low,high); 
  tauTree->Draw(variable1+">>+Num",RecoCut3);

  std::cout<<"Num_pu0: "<<Num->GetEntries()<<std::endl;
  Num->Divide(Denom);
  /// second
  TH1* Denom1;
  Denom1 = new TH1F("Denom1","Denom1",bins,low,high);
  Denom1->Sumw2();
  tauTree1->Draw(variable+">>+Denom1",GenCut);

  std::cout<<"Denom_pu140: "<<Denom1->GetEntries()<<std::endl;
  TH1* Num1;
  Num1 = new TH1F("Num1","Num1",bins,low,high);
  tauTree1->Draw(variable1+">>+Num1",RecoCut3);

  std::cout<<"Num_pu140: "<<Num1->GetEntries()<<std::endl;
  Num1->Divide(Denom1);
  Num1->SetMarkerStyle(markerstyle1);
  Num1->SetMarkerColor(color1);
  ////
  /// third
  TH1* Denom2;
  Denom2 = new TH1F("Denom2","Denom2",bins,low,high);
  Denom2->Sumw2();
  tauTree2->Draw(variable+">>+Denom2",GenCut);

  std::cout<<"Denom_pu200: "<<Denom2->GetEntries()<<std::endl;
  TH1* Num2;
  Num2 = new TH1F("Num2","Num2",bins,low,high);
  tauTree2->Draw(variable1+">>+Num2",RecoCut3);

  std::cout<<"Num_pu200: "<<Num2->GetEntries()<<std::endl;
  Num2->Divide(Denom2);

  Num2->SetMarkerStyle(markerstyle2);
  Num2->SetMarkerColor(color2);
  ////
 /// fourth
  TH1* Denom3;
  Denom3 = new TH1F("Denom3","Denom3",bins,low,high);
  Denom3->Sumw2();
  tauTree3->Draw(variable+">>+Denom3",GenCut);

  TH1* Num3;
  Num3 = new TH1F("Num3","Num3",bins,low,high);
  tauTree3->Draw(variable1+">>+Num3",RecoCut3);

  Num3->Divide(Denom3);

  Num3->SetMarkerStyle(markerstyle3);
  //Num3->SetMarkerColor(color3);
  
  gStyle->SetErrorX(0.5);

  Num->GetXaxis()->SetTitle(xaxis);
 // Num->GetYaxis()->SetTitleSize(0.06);
 // Num->GetYaxis()->SetTitleFont(42);
 // Num->GetYaxis()->SetLabelOffset(0.007);
 // Num->GetYaxis()->SetLabelSize(0.05);
  Num->GetYaxis()->SetTitle("Efficiency");
  Num->GetYaxis()->SetTitleOffset(0.9);
  Num->SetMarkerStyle(markerstyle);
  Num->SetMarkerColor(color);
  Num->Draw("P");
  Num1->Draw("P same");
  Num2->Draw("P same");
  //Num3->Draw("P same");

  Num->SetFillStyle(1001);
  Num->SetLineWidth((short)1.5);
  Num->SetMaximum(1.4);
  Num->SetMinimum(0);

 // TLatex latex;
 // latex.SetTextFont(42);
 // latex.SetTextAngle(0);
 // latex.SetTextColor(kBlack);    
 // latex.SetTextSize(0.1);    
 // latex.SetTextAlign(12);
 // latex.DrawLatex(.20,.345,"VBF H");
  TLegend *leg = new TLegend(.35, .15, .55, .32,legLabel,"nbNDC");  
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetShadowColor(kWhite);
  //leg->SetFillColor(kWhite);
  leg->SetTextSize(0.025);
  //leg->AddEntry(pu200_eff,"PF Charged Iso All","PL");
  //leg->AddEntry((TObject*)0,"PFChargedIso < 2 GeV","");
  leg->AddEntry(Num,"PU0","PL");
  leg->AddEntry(Num1,"PU140","PL");
  leg->AddEntry(Num2,"PU200","PL");//MediumIsolationMVArun2v1DBnewDMwLT 
  //leg->AddEntry(Num3,"PU30(VBF Phase1)","PL");
  //leg->SetHeader("PFChargedIso < 2 GeV");
  //leg->AddEntry(Num3,"VBF PhaseII&& ChargedIso<2GeV","PL");
  leg->Draw();
  
   canv->cd();
  //Writing the lumi info
  CMS_lumi(canv, iPeriod, iPos );

  canv->Update();
  canv->RedrawAxis();
  canv->GetFrame()->Draw();
  
  canv->SaveAs("Validation_May9/"+canvName+".pdf",".pdf");
  canv->Print("Validation_May9/"+canvName+".png",".png");
  //canv->SaveAs(outFileName+".pdf");
}
