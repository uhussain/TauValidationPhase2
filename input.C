// -*- C++ -*-
//
// Package:     step3
// File:        input.C
//
// Author:      Martin Flechl
// Created:     2015, June 3rd
//
#define input_cxx
#include "input.h"
#include "util.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TGraphAsymmErrors.h>

#include<iostream>

const int DO62X=0;
//const float VZ_MAX=5;
const float VZ_MAX=1e6;
const float VZ_MIN=0;
const float MAX_ISO=2;
//const float WP=0.05; //find WP with 5% fake id eff
const float WP=0.10; //find WP with 10% fake id eff
const float MAX_DR=0.2;
const float MAX_DR2=MAX_DR*MAX_DR;

const float EFF_MIN_PT=22;
const float EFF_MAX_ETA=4.0;

const double BS_SIZE=33;

const int DEBUG=0;

Int_t input::Loop(double isocut, double vzmax, double vzmin)
{
  if (isocut<0) isocut=MAX_ISO;
  if (vzmax<0) vzmax=VZ_MAX;
  if (vzmin<0) vzmin=VZ_MIN;
  //std::cout << "ISOCUT:  " << isocut << std::endl;
  //std::cout << "VZMAX:   " << vzmax << std::endl;
  //std::cout << "VZMIN:   " << vzmin << std::endl;

  TH1::SetDefaultSumw2(kTRUE);

  float eta1=0,eta2=0,phi1=0,phi2=0;

  TFile *f = new TFile(ofilename,"RECREATE");
  TTree *tree = new TTree("Events","info");

  TH1D *h_gen_tauvis_eta20 = new TH1D("h_gen_tauvis_eta20","",12,-2.4,+2.4);
  TH1D *h_gen_tauvis_eta20_reco = (TH1D*)h_gen_tauvis_eta20->Clone("h_gen_tauvis_eta20_reco");
  TH1D *h_gen_tauvis_eta20_id = (TH1D*)h_gen_tauvis_eta20->Clone("h_gen_tauvis_eta20_id");

  TH1D *h_gen_tauvis_pt23 = new TH1D("h_gen_tauvis_pt23","",12,0,120);
  TH1D *h_gen_tauvis_pt23_reco = (TH1D*)h_gen_tauvis_pt23->Clone("h_gen_tauvis_pt23_reco");
  TH1D *h_gen_tauvis_pt23_id = (TH1D*)h_gen_tauvis_pt23->Clone("h_gen_tauvis_pt23_id");

  TH1D *h_fake_pt23_reco = (TH1D*)h_gen_tauvis_pt23->Clone("h_fake_pt23_reco");
  TH1D *h_fake_pt23_id = (TH1D*)h_gen_tauvis_pt23->Clone("h_fake_pt23_id");
  TH1D *h_fake_eta20_reco = (TH1D*)h_gen_tauvis_eta20->Clone("h_fake_eta20_reco");
  TH1D *h_fake_eta20_id = (TH1D*)h_gen_tauvis_eta20->Clone("h_fake_eta20_id");

  int c_all=0,c_reco=0,c_id=0;
  int c_fake_reco=0, c_fake_id=0;

  tree->Branch("tau_pt",&floats_tau_pt_nt_obj);
  tree->Branch("tau_eta",&floats_tau_eta_nt_obj);
  tree->Branch("tau_phi",&floats_tau_phi_nt_obj);
  if (!DO62X){
    tree->Branch("tau_againstElMVA_l",&floats_tau_againstElectronLooseMVA6_nt_obj);
    tree->Branch("tau_againstElMVA",&floats_tau_againstElectronMVA6Raw_nt_obj);
    tree->Branch("tau_againstElMVA_t",&floats_tau_againstElectronTightMVA6_nt_obj);
  }/* else{
    tree->Branch("tau_againstElMVA_l",&floats_tau_againstElectronLooseMVA5_nt_obj);
    tree->Branch("tau_againstElMVA",&floats_tau_againstElectronMVA5Raw_nt_obj);
    tree->Branch("tau_againstElMVA_t",&floats_tau_againstElectronTightMVA5_nt_obj);
    }*/
  tree->Branch("tau_againstMu_l",&floats_tau_againstMuonLoose3_nt_obj);
  tree->Branch("tau_againstMu_t",&floats_tau_againstMuonTight3_nt_obj);
  tree->Branch("tau_DBIso",&floats_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits_nt_obj);
  tree->Branch("tau_MVAIso_dr03old",&floats_tau_byIsolationMVArun2v1DBdR03oldDMwLTraw_nt_obj);
  tree->Branch("tau_MVAIso",&floats_tau_byIsolationMVArun2v1DBnewDMwLTraw_nt_obj);
  tree->Branch("tau_DBIso_l",&floats_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits_nt_obj);
  tree->Branch("tau_DBIso_m",&floats_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits_nt_obj);
  tree->Branch("tau_MVAIso_m_dr03old",&floats_tau_byMediumIsolationMVArun2v1DBdR03oldDMwLT_nt_obj);
  tree->Branch("tau_MVAIso_m",&floats_tau_byMediumIsolationMVArun2v1DBnewDMwLT_nt_obj);
  tree->Branch("tau_byPhotonIsoPt",&floats_tau_byPhotonPtSumOutsideSignalCone_nt_obj);
  tree->Branch("tau_chargedIsoPtSum",&floats_tau_chargedIsoPtSum_nt_obj);
  tree->Branch("tau_chargedIsoPtSim_dr03",&floats_tau_chargedIsoPtSumdR03_nt_obj);
  tree->Branch("tau_dm",&floats_tau_decayMode_nt_obj);
  tree->Branch("tau_dmf",&floats_tau_decayModeFinding_nt_obj);
  tree->Branch("tau_dmf_new",&floats_tau_decayModeFindingNewDMs_nt_obj);
  tree->Branch("tau_footprintCorr",&floats_tau_footprintCorrection_nt_obj);
  tree->Branch("tau_footprintCorr_dr03",&floats_tau_footprintCorrectiondR03_nt_obj);
  tree->Branch("tau_neutralIsoPtSum",&floats_tau_neutralIsoPtSum_nt_obj);
  tree->Branch("tau_neutralIsoPtSum_dr03",&floats_tau_neutralIsoPtSumdR03_nt_obj);
  tree->Branch("tau_pdgId",&floats_tau_pdgId_nt_obj);
  tree->Branch("tau_genSize",&floats_tau_genSize_nt_obj);
  tree->Branch("tau_photonIsoPtSum",&floats_tau_photonPtSumOutsideSignalCone_nt_obj);
  tree->Branch("tau_photonIsoPtSum_dr03",&floats_tau_photonPtSumOutsideSignalConedR03_nt_obj);
  tree->Branch("tau_puCorrPtSum",&floats_tau_puCorrPtSum_nt_obj);
  tree->Branch("tau_nChHadr",&floats_tau_NsignalTauChargedHadronCandidates_nt_obj);
  tree->Branch("tau_nPi0",&floats_tau_NsignalTauChargedHadronCandidates_nt_obj);
  tree->Branch("tau_vz",&floats_tau_vz_nt_obj);

  Int_t pu_num_true; tree->Branch("pu_num_true",&pu_num_true); //getPU_NumInteractions()
  Int_t pu_num_sample; tree->Branch("pu_num_sample",&pu_num_sample); //getTrueNumInteractions()
  //  tree->Branch("pu_num_goodvtx",&recoVertexs_goodprimvtx__nt_obj_);
  //  tree->Branch("pu_num_goodvtx_bs",&recoVertexs_goodprimvtxbs__nt_obj_);
  tree->Branch("pu_num_allvtx",&recoVertexs_allprimvtx__nt_obj_);
  //  tree->Branch("pu_num_allvtx_bs",&recoVertexs_allprimvtxbs__nt_obj_);

  tree->Branch("vtx1_z",&recoVertexs_allprimvtx__nt_obj_position__fCoordinates_fZ[0]);

  tree->Branch("gen_pt",&floats_gen_pt_nt_obj);
  tree->Branch("gen_eta",&floats_gen_eta_nt_obj);
  tree->Branch("gen_phi",&floats_gen_phi_nt_obj);
  tree->Branch("gen_m",&floats_gen_m_nt_obj);
  tree->Branch("gen_vz",&floats_gen_vz_nt_obj);
  tree->Branch("gen_pdgId",&floats_gen_pdgId_nt_obj);
  tree->Branch("gen_status",&floats_gen_status_nt_obj);
  tree->Branch("gen_isLast",&floats_gen_isLastCopy_nt_obj);

  std::vector<float> gen_tau_eta, gen_tau_pt;
  std::vector<float> gen_tauvis_eta, gen_tauvis_eta20, gen_tauvis_pt, gen_tauvis_pt23, gen_tauvis_phi;
  std::vector<float> gen_tauvis_eta_reco, gen_tauvis_eta20_reco, gen_tauvis_pt_reco, gen_tauvis_pt23_reco;
  std::vector<float> gen_tauvis_eta_id, gen_tauvis_eta20_id, gen_tauvis_pt_id, gen_tauvis_pt23_id;
  std::vector<float> gen_tauvis_eta_allvz, gen_tauvis_phi_allvz;
  std::vector<float> gen_lep_eta, gen_lep_phi;
  std::vector<float> fake_pt23_id,fake_pt23_reco,fake_eta20_id, fake_eta20_reco;
  std::vector<float> tau_eta_match, tau_pt_match;
  std::vector<float> tau_dm_match, tau_dmf_match, tau_dmf_new_match;
  std::vector<float> tau_chargedIsoPtSum_match, tau_neutralIsoPtSum_match, tau_DBIso_match, tau_DBIso_m_match, tau_MVAIso_match, tau_MVAIso_m_match;
  std::vector<float> tau_againstElMVA_match,tau_againstElMVA_l_match;
  std::vector<float> gen_tau_mom_pdgid, gen_tau_mom_pt, gen_tau_mom_eta, gen_tau_mom_phi, gen_tau_mom_m;
  std::vector<std::vector<float> > gen_tau_dau_pdgid, gen_tau_dau_pt, gen_tau_dau_eta, gen_tau_dau_phi, gen_tau_dau_m;

  tree->Branch("gen_tau_dau_pdgid",&gen_tau_dau_pdgid);
  tree->Branch("gen_tau_dau_pt",&gen_tau_dau_pt);
  tree->Branch("gen_tau_dau_eta",&gen_tau_dau_eta);
  tree->Branch("gen_tau_dau_phi",&gen_tau_dau_phi);
  tree->Branch("gen_tau_dau_m",&gen_tau_dau_m);

  tree->Branch("gen_tau_mom_pdgid",&gen_tau_mom_pdgid);
  tree->Branch("gen_tau_mom_pt",&gen_tau_mom_pt);
  tree->Branch("gen_tau_mom_eta",&gen_tau_mom_eta);
  tree->Branch("gen_tau_mom_phi",&gen_tau_mom_phi);
  tree->Branch("gen_tau_mom_m",&gen_tau_mom_m);

  tree->Branch("gen_tau_eta",&gen_tau_eta);
  tree->Branch("gen_tau_pt",&gen_tau_pt);
  tree->Branch("gen_tauvis_eta",&gen_tauvis_eta);
  tree->Branch("gen_tauvis_eta20",&gen_tauvis_eta20);
  tree->Branch("gen_tauvis_pt",&gen_tauvis_pt);
  tree->Branch("gen_tauvis_pt23",&gen_tauvis_pt23);

  tree->Branch("gen_tauvis_eta_reco",&gen_tauvis_eta_reco);
  tree->Branch("gen_tauvis_eta20_reco",&gen_tauvis_eta20_reco);
  tree->Branch("gen_tauvis_pt_reco",&gen_tauvis_pt_reco);
  tree->Branch("gen_tauvis_pt23_reco",&gen_tauvis_pt23_reco);
  tree->Branch("gen_tauvis_eta_id",&gen_tauvis_eta_id);
  tree->Branch("gen_tauvis_eta20_id",&gen_tauvis_eta20_id);
  tree->Branch("gen_tauvis_pt_id",&gen_tauvis_pt_id);
  tree->Branch("gen_tauvis_pt23_id",&gen_tauvis_pt23_id);

  tree->Branch("fake_pt23_id",&fake_pt23_id);
  tree->Branch("fake_pt23_reco",&fake_pt23_reco);
  tree->Branch("fake_eta20_id",&fake_eta20_id);
  tree->Branch("fake_eta20_reco",&fake_eta20_reco);

  tree->Branch("tau_pt_match",&tau_pt_match);
  tree->Branch("tau_eta_match",&tau_eta_match);
  tree->Branch("tau_dm_match",&tau_dm_match);
  tree->Branch("tau_dmf_match",&tau_dmf_match);
  tree->Branch("tau_dmf_new_match",&tau_dmf_new_match);
  tree->Branch("tau_chargedIsoPtSum_match",&tau_chargedIsoPtSum_match);
  tree->Branch("tau_neutralIsoPtSum_match",&tau_neutralIsoPtSum_match);
  tree->Branch("tau_DBIso_match",&tau_DBIso_match);
  tree->Branch("tau_DBIso_m_match",&tau_DBIso_m_match);
  tree->Branch("tau_MVAIso_match",&tau_MVAIso_match);
  tree->Branch("tau_MVAIso_m_match",&tau_MVAIso_m_match);
  tree->Branch("tau_againstElMVA_match",&tau_againstElMVA_match);
  tree->Branch("tau_againstElMVA_l_match",&tau_againstElMVA_l_match);
  //  tree->Branch("",&);

  tree->Branch("jet_pt",&floats_jet_pt_nt_obj);
  tree->Branch("jet_eta",&floats_jet_eta_nt_obj);
  tree->Branch("jet_nChHad",&floats_jet_chargedHadronMultiplicity_nt_obj);
  tree->Branch("jet_hadronFlavour",&floats_jet_hadronFlavour_nt_obj);
  tree->Branch("jet_partonFlavour",&floats_jet_partonFlavour_nt_obj);
  tree->Branch("jet_nNeHad",&floats_jet_neutralHadronMultiplicity_nt_obj);
  tree->Branch("jet_nPhoton",&floats_jet_photonMultiplicity_nt_obj);

  //  std::cout << "Local density at z=" << BS_SIZE <<"mm: " << this->get_local_density(BS_SIZE) << std::endl;

  Long64_t nentries = fChain->GetEntries();

  std::cout << "Events: " << nentries << std::endl;
   
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    if ( (jentry+1)%20000 == 0 ) std::cout << "Event " << jentry+1 << std::endl;
    //    if (jentry>=20) break;
    if (DEBUG) std::cout << "A" << std::endl;
    fChain->GetEntry(jentry);
    if (DEBUG) std::cout << "B" << std::endl;

    if (!DO62X){
      for (Int_t ib=0; ib<kMaxPileupSummaryInfos_truemu__nt_obj; ib++){
	if (PileupSummaryInfos_truemu__nt_obj_bunchCrossing_[ib]==0){ 
	  pu_num_true=PileupSummaryInfos_truemu__nt_obj_num_PU_vertices_[ib]; 
	  pu_num_sample=PileupSummaryInfos_truemu__nt_obj_TrueNumInteractions_[ib];
	}
      }
    }

    tau_pt_match.clear();
    tau_eta_match.clear();
    tau_dm_match.clear();
    tau_dmf_match.clear();
    tau_dmf_new_match.clear();
    tau_chargedIsoPtSum_match.clear();
    tau_neutralIsoPtSum_match.clear();
    tau_DBIso_match.clear();
    tau_DBIso_m_match.clear();
    tau_MVAIso_match.clear();
    tau_MVAIso_m_match.clear();
    tau_againstElMVA_match.clear();
    tau_againstElMVA_l_match.clear();
    gen_tau_eta.clear();
    gen_tau_pt.clear();
    gen_tau_dau_pdgid.clear();
    gen_tau_dau_pt.clear();
    gen_tau_dau_eta.clear();
    gen_tau_dau_phi.clear();
    gen_tau_dau_m.clear();
    gen_tau_mom_pdgid.clear();
    gen_tau_mom_pt.clear();
    gen_tau_mom_eta.clear();
    gen_tau_mom_phi.clear();
    gen_tau_mom_m.clear();
    gen_tauvis_eta.clear();
    gen_tauvis_phi.clear();
    gen_tauvis_pt.clear();
    gen_tauvis_eta20.clear();
    gen_tauvis_pt23.clear();
    gen_tauvis_eta_reco.clear();
    gen_tauvis_pt_reco.clear();
    gen_tauvis_eta20_reco.clear();
    gen_tauvis_pt23_reco.clear();
    gen_tauvis_eta_id.clear();
    gen_tauvis_pt_id.clear();
    gen_tauvis_eta20_id.clear();
    gen_tauvis_pt23_id.clear();
    gen_tauvis_eta_allvz.clear();
    gen_tauvis_phi_allvz.clear();
    gen_lep_eta.clear();
    gen_lep_phi.clear();
    fake_pt23_reco.clear();
    fake_eta20_reco.clear();
    fake_pt23_id.clear();
    fake_eta20_id.clear();


    if (DEBUG) std::cout << "C" << std::endl;

    //get true taus
    //    std:vector<TLorentzVector> v_tauvis;
    for (unsigned it=0; it<floats_gen_pdgId_nt_obj.size(); it++){
      if (DEBUG) std::cout << "C " << it << std::endl;
      //    if ( std::abs(floats_gen_pdgId_nt_obj.at(it))==15 && floats_gen_isLastCopy_nt_obj.at(it) && floats_gen_isFromHardProcessDecayed_nt_obj.at(it)  ){

      if ( ( std::abs(floats_gen_pdgId_nt_obj.at(it))==11 || std::abs(floats_gen_pdgId_nt_obj.at(it))==13 ) ){ 
	//&& ( DO62X || ( floats_gen_isLastCopy_nt_obj.at(it) && floats_gen_isFromHardProcessDecayed_nt_obj.at(it) ) ) ){
	if ( floats_gen_status_nt_obj.at(it)==1 ){
	  gen_lep_eta.push_back( floats_gen_eta_nt_obj.at(it) );
	  gen_lep_phi.push_back( floats_gen_phi_nt_obj.at(it) );
	}
      }

      if ( std::abs(floats_gen_pdgId_nt_obj.at(it))==15 && (floats_gen_status_nt_obj.at(it)<=2) && ( DO62X || ( floats_gen_isLastCopy_nt_obj.at(it) && floats_gen_isFromHardProcessDecayed_nt_obj.at(it) ) ) ){
	if (DEBUG) std::cout << "CA " << it << std::endl;

	std::vector<float> dau_pdgid, dau_pt, dau_eta, dau_phi, dau_m;
	get_tau_vec(it, dau_pdgid, dau_pt, dau_eta, dau_phi, dau_m); 
	gen_tau_dau_pdgid.push_back(dau_pdgid);
	gen_tau_dau_pt.push_back(dau_pt);
	gen_tau_dau_eta.push_back(dau_eta);
	gen_tau_dau_phi.push_back(dau_phi);
	gen_tau_dau_m.push_back(dau_m);
	if (DEBUG) std::cout << "CAA " << it << std::endl;

	if ( ! DO62X ){
	  gen_tau_mom_pdgid.push_back( floats_gen_mother0pdgId_nt_obj.at(it) );
	  if (DEBUG) std::cout << "CAB " << it << std::endl;
	  gen_tau_mom_pt.push_back(    floats_gen_mother0pt_nt_obj.at(it) )  ;
	  gen_tau_mom_eta.push_back(   floats_gen_mother0eta_nt_obj.at(it)   );
	  gen_tau_mom_phi.push_back(   floats_gen_mother0phi_nt_obj.at(it)   );
	  gen_tau_mom_m.push_back(     floats_gen_mother0m_nt_obj.at(it)     );
	}

	if (DEBUG && jentry<5) for (unsigned j=0; j<dau_pdgid.size(); j++){
	  std::cout << j << " " << dau_pdgid.at(j) << " X " << dau_pt.at(j) <<" X " << dau_eta.at(j) <<" X " << dau_phi.at(j) <<" X " << dau_m.at(j) <<std::endl;
	  std::cout << j << "|" << gen_tau_dau_pdgid.back().at(j) << " X " << gen_tau_dau_pt.back().at(j) <<" X " << gen_tau_dau_eta.back().at(j) <<" X " << gen_tau_dau_phi.back().at(j) <<" X " << gen_tau_dau_m.back().at(j) <<std::endl;
	}
	if (DEBUG) std::cout << "CB " << it << std::endl;

	int lepDecay=0;
	for (unsigned it2=0; it2<dau_pdgid.size(); it2++){
	  if ( std::abs(dau_pdgid.at(it2))==12 || ( std::abs(dau_pdgid.at(it2))==14 ) ){ lepDecay=1; break; }
	}
	if (DEBUG) std::cout << "CC " << it << std::endl;

	if (lepDecay) continue;

	for (unsigned it2=0; it2<dau_pdgid.size(); it2++){
	  if (DEBUG) std::cout << "CD " << it2 << std::endl;

	  if ( std::abs(dau_pdgid.at(it2))==16 &&  (dau_pdgid.at(it2)*floats_gen_pdgId_nt_obj.at(it)) >0   ){
	    TLorentzVector tau; tau.SetPtEtaPhiM( floats_gen_pt_nt_obj.at(it)  , floats_gen_eta_nt_obj.at(it)  , floats_gen_phi_nt_obj.at(it)  , floats_gen_m_nt_obj.at(it)  );
            TLorentzVector nu;  nu.SetPtEtaPhiM(  dau_pt.at(it2) , dau_eta.at(it2) , dau_phi.at(it2) , dau_m.at(it2) );
            TLorentzVector tauvis=tau-nu;
	    //	    v_tauvis.push_back(tauvis);

	    //	    if (jentry<500){ std::cout << fabs(floats_gen_vz_nt_obj.at(it)) << " < " << vzmax/10. << "  \t" << fabs(floats_gen_vz_nt_obj.at(it)) << " > " << vzmin/10. << std::endl; }
	    if ( fabs(floats_gen_vz_nt_obj.at(it)) > vzmax/10. || fabs(floats_gen_vz_nt_obj.at(it)) < vzmin/10. ){
	      gen_tauvis_eta_allvz.push_back( tauvis.Eta() );
	      gen_tauvis_phi_allvz.push_back( tauvis.Phi() );
	      continue;
	    }

            gen_tau_eta.push_back(floats_gen_eta_nt_obj.at(it));
            gen_tau_pt.push_back(floats_gen_pt_nt_obj.at(it));
            gen_tauvis_pt.push_back( tauvis.Pt() );
            gen_tauvis_eta.push_back( tauvis.Eta() );
            gen_tauvis_phi.push_back( tauvis.Phi() );
            if ( tauvis.Pt()>20 ){         gen_tauvis_eta20.push_back( tauvis.Eta() ); h_gen_tauvis_eta20->Fill( tauvis.Eta() ); }
            if ( fabs(tauvis.Eta())<4.0 ){ gen_tauvis_pt23.push_back( tauvis.Pt() );   h_gen_tauvis_pt23->Fill(  tauvis.Pt() ); }
            if ( tauvis.Pt()>EFF_MIN_PT && fabs(tauvis.Eta())<EFF_MAX_ETA ){ c_all++; }
	    break;
	  }
	}
      }
    }//end: for over gen particles
    if (DEBUG) std::cout << "D" << std::endl;

    //for efficiency plots: pick out true taus that have been reconstructed   
    for (unsigned it=0; it<gen_tauvis_eta.size(); it++){
      if(gen_tauvis_pt.at(it) > 15 && abs(gen_tauvis_eta.at(it)) < 4.0){
      int found_match =0;
      for (unsigned it2=0; it2<floats_tau_eta_nt_obj.size(); it2++){
	if ( deltaR2( gen_tauvis_eta.at(it) , gen_tauvis_phi.at(it) , floats_tau_eta_nt_obj.at(it2) , floats_tau_phi_nt_obj.at(it2) ) < MAX_DR2 ){
    if(floats_tau_pt_nt_obj.at(it2) > 18 && abs(floats_tau_eta_nt_obj.at(it2)) < 4.0) {
	  gen_tauvis_eta_reco.push_back(  gen_tauvis_eta.at(it) );
	  gen_tauvis_pt_reco.push_back(   gen_tauvis_pt.at(it)  );
	  tau_pt_match.push_back( floats_tau_pt_nt_obj.at(it2) );
	  tau_eta_match.push_back( floats_tau_eta_nt_obj.at(it2) );
	  tau_dm_match.push_back( floats_tau_decayMode_nt_obj.at(it2) );
	  tau_dmf_match.push_back( floats_tau_decayModeFinding_nt_obj.at(it2) );
	  tau_dmf_new_match.push_back( floats_tau_decayModeFindingNewDMs_nt_obj.at(it2) );
	  tau_chargedIsoPtSum_match.push_back( floats_tau_chargedIsoPtSum_nt_obj.at(it2) );
	  tau_neutralIsoPtSum_match.push_back( floats_tau_neutralIsoPtSum_nt_obj.at(it2) );
	  tau_DBIso_match.push_back( floats_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits_nt_obj.at(it2) );
	  if ( gen_tauvis_pt.at(it)>20 ){         gen_tauvis_eta20_reco.push_back( gen_tauvis_eta.at(it) ); h_gen_tauvis_eta20_reco->Fill( gen_tauvis_eta.at(it) ); }
	  if ( fabs(gen_tauvis_eta.at(it))<4.0 ){ gen_tauvis_pt23_reco.push_back( gen_tauvis_pt.at(it) );   h_gen_tauvis_pt23_reco->Fill(  gen_tauvis_pt.at(it) );}
	  if ( gen_tauvis_pt.at(it)>EFF_MIN_PT && fabs(gen_tauvis_eta.at(it))<EFF_MAX_ETA ){ c_reco++; }

	  if ( !this->is_tau_id(it2,isocut) ) continue;
	  gen_tauvis_eta_id.push_back(  gen_tauvis_eta.at(it) );
	  gen_tauvis_pt_id.push_back(   gen_tauvis_pt.at(it)  );
	  if ( gen_tauvis_pt.at(it)>20 ){         gen_tauvis_eta20_id.push_back( gen_tauvis_eta.at(it) ); h_gen_tauvis_eta20_id->Fill( gen_tauvis_eta.at(it) ); }
	  if ( fabs(gen_tauvis_eta.at(it))<4.0 ){ gen_tauvis_pt23_id.push_back( gen_tauvis_pt.at(it) );   h_gen_tauvis_pt23_id->Fill(  gen_tauvis_pt.at(it) );}
	  if ( gen_tauvis_pt.at(it)>EFF_MIN_PT && fabs(gen_tauvis_eta.at(it))<EFF_MAX_ETA ){ c_id++; }
	
    found_match=1;
  }
      }
	if (found_match) break;
      }
      }
      }

    //  for (unsigned it=0; it<floats_tau_pdgId_nt_obj.size(); it++){
    //    if ( fabs( floats_tau_pdgId_nt_obj.at(it) ) == 15 && floats_tau_genSize_nt_obj.at(it)>0 ){
    for (unsigned it=0; it<floats_tau_eta_nt_obj.size(); it++){
      int found_match=0;
      for (unsigned it2=0; it2<gen_tauvis_eta.size(); it2++){
	if ( deltaR2(gen_tauvis_eta.at(it2),gen_tauvis_phi.at(it2),floats_tau_eta_nt_obj.at(it),floats_tau_phi_nt_obj.at(it)) <MAX_DR2 ){

	//  tau_pt_match.push_back( floats_tau_pt_nt_obj.at(it) );
	//  tau_eta_match.push_back( floats_tau_eta_nt_obj.at(it) );
	//  tau_dm_match.push_back( floats_tau_decayMode_nt_obj.at(it) );
	//  tau_dmf_match.push_back( floats_tau_decayModeFinding_nt_obj.at(it) );
	//  tau_dmf_new_match.push_back( floats_tau_decayModeFindingNewDMs_nt_obj.at(it) );
	//  tau_chargedIsoPtSum_match.push_back( floats_tau_chargedIsoPtSum_nt_obj.at(it) );
	//  tau_neutralIsoPtSum_match.push_back( floats_tau_neutralIsoPtSum_nt_obj.at(it) );
	//  tau_DBIso_match.push_back( floats_tau_byCombinedIsolationDeltaBetaCorrRaw3Hits_nt_obj.at(it) );
	  if (DO62X){
	    tau_DBIso_m_match.push_back( floats_tau_byLooseCombinedIsolationDeltaBetaCorr3Hits_nt_obj.at(it) );
	    //	  tau_againstElMVA_match.push_back( floats_tau_againstElectronMVA5Raw_nt_obj.at(it) );
	    //	  tau_againstElMVA_l_match.push_back( floats_tau_againstElectronLooseMVA5_nt_obj.at(it) );
	  }
	  else{
	    tau_DBIso_m_match.push_back( floats_tau_byMediumCombinedIsolationDeltaBetaCorr3Hits_nt_obj.at(it) );
	    tau_againstElMVA_match.push_back( floats_tau_againstElectronMVA6Raw_nt_obj.at(it) );
	    tau_againstElMVA_l_match.push_back( floats_tau_againstElectronLooseMVA6_nt_obj.at(it) );
	    tau_MVAIso_match.push_back( floats_tau_byIsolationMVArun2v1DBnewDMwLTraw_nt_obj.at(it) );
	    tau_MVAIso_m_match.push_back( floats_tau_byMediumIsolationMVArun2v1DBnewDMwLT_nt_obj.at(it) );
	  }
	  found_match=1;
	}
	if (found_match) break;
      }
      if ( !found_match ){ //potential fake!
	if ( fabs(floats_tau_vz_nt_obj.at(it)) > vzmax/10. ) continue; //speed up execution by checking this already here
	if ( fabs(floats_tau_vz_nt_obj.at(it)) < vzmin/10. ) continue; //speed up execution by checking this already here
        if ( fabs( floats_tau_vz_nt_obj.at(it) - recoVertexs_allprimvtx__nt_obj_position__fCoordinates_fZ[0] )>0.2 ) continue; //pu rej

	for (unsigned it2=0; it2<gen_tauvis_eta_allvz.size(); it2++){ //check if match to gen tau outside of vz limit
	  if ( deltaR2(gen_tauvis_eta_allvz.at(it2),gen_tauvis_phi_allvz.at(it2),floats_tau_eta_nt_obj.at(it),floats_tau_phi_nt_obj.at(it)) <MAX_DR2 ){ found_match=1; break; }
	}
	if ( found_match ) continue;
	for (unsigned it2=0; it2<gen_lep_eta.size(); it2++){ //check if match to deltaR to gen tau
	  if ( deltaR2(gen_lep_eta.at(it2),gen_lep_phi.at(it2),floats_tau_eta_nt_obj.at(it),floats_tau_phi_nt_obj.at(it)) <MAX_DR2 ){ found_match=1; break; }
	}
	if ( found_match ) continue;

	//eta plots
        int isTauId=-1;
        if ( floats_tau_pt_nt_obj.at(it)>20 ){
	  fake_eta20_reco.push_back( floats_tau_eta_nt_obj.at(it) );
	  h_fake_eta20_reco->Fill( floats_tau_eta_nt_obj.at(it) );
	  if ( floats_tau_pt_nt_obj.at(it)>EFF_MIN_PT && fabs(floats_tau_eta_nt_obj.at(it))<EFF_MAX_ETA ){ c_fake_reco++; }
	  isTauId=0;
	  if ( this->is_tau_id(it,isocut) ){
	    fake_eta20_id.push_back( floats_tau_eta_nt_obj.at(it) );
	    h_fake_eta20_id->Fill( floats_tau_eta_nt_obj.at(it) );
	    if ( floats_tau_pt_nt_obj.at(it)>EFF_MIN_PT && fabs(floats_tau_eta_nt_obj.at(it))<EFF_MAX_ETA ){ c_fake_id++; }
	    isTauId=1;
	  }
	}
	//pt plots
	if ( fabs(floats_tau_eta_nt_obj.at(it))<4.0 ){
	  fake_pt23_reco.push_back( floats_tau_pt_nt_obj.at(it) );
	  h_fake_pt23_reco->Fill( floats_tau_pt_nt_obj.at(it) );
          if ( isTauId==1 || ( isTauId==-1 && this->is_tau_id(it,isocut) ) ){ //avoid checking the same tau twice
	    fake_pt23_id.push_back( floats_tau_pt_nt_obj.at(it) );
	    h_fake_pt23_id->Fill( floats_tau_pt_nt_obj.at(it) );
	  }
	}
      }
    }//for loop over reco taus

    tree->Fill();
    //    tree->FlushBaskets();
  }

  tree->Write();

  h_gen_tauvis_eta20_reco->Write();
  h_gen_tauvis_pt23_reco->Write();
  h_gen_tauvis_eta20_id->Write();
  h_gen_tauvis_pt23_id->Write();
  h_gen_tauvis_eta20->Write();
  h_gen_tauvis_pt23->Write();

  h_fake_pt23_reco->Write();
  h_fake_eta20_reco->Write();
  h_fake_pt23_id->Write();
  h_fake_eta20_id->Write();

  TGraphAsymmErrors *g_eff_eta20_reco=new TGraphAsymmErrors(h_gen_tauvis_eta20_reco, h_gen_tauvis_eta20);
  TGraphAsymmErrors *g_eff_pt23_reco=new TGraphAsymmErrors(h_gen_tauvis_pt23_reco, h_gen_tauvis_pt23);
  TGraphAsymmErrors *g_eff_eta20_id=new TGraphAsymmErrors(h_gen_tauvis_eta20_id, h_gen_tauvis_eta20);
  TGraphAsymmErrors *g_eff_pt23_id=new TGraphAsymmErrors(h_gen_tauvis_pt23_id, h_gen_tauvis_pt23);

  g_eff_eta20_reco->SetName("g_eff_eta20_reco");
  g_eff_pt23_reco->SetName("g_eff_pt23_reco");
  g_eff_eta20_id->SetName("g_eff_eta20_id");
  g_eff_pt23_id->SetName("g_eff_pt23_id");

  g_eff_eta20_id->Write();
  g_eff_pt23_id->Write();
  g_eff_eta20_reco->Write();
  g_eff_pt23_reco->Write();

  std::cout << c_all << " " << c_reco << " " << c_id << " " << c_fake_reco << " " << c_fake_id << std::endl;
  std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
  std::cout.precision(3);
  if ( c_all>0 ){
    std::cout << "Reco    eff for pt>" << EFF_MIN_PT << " and |eta|<"<< EFF_MAX_ETA << ":\t" << c_reco*1.0/c_all << " +/- " << sqrt(c_reco)*1.0/c_all  << std::endl;
    std::cout << "Reco+ID eff for pt>" << EFF_MIN_PT << " and |eta|<"<< EFF_MAX_ETA << ":\t" << c_id*1.0/c_all   << " +/- " << sqrt(c_id)*1.0/c_all  << std::endl;
  }
  if ( c_reco>0 ){
    std::cout << "ID      eff for pt>" << EFF_MIN_PT << " and |eta|<"<< EFF_MAX_ETA << ":\t" << c_id*1.0/c_reco   << " +/- " << sqrt(c_id)*1.0/c_reco  << std::endl;
  }
  if ( c_fake_reco>0 ){
    std::cout << "ID      FR  for pt>" << EFF_MIN_PT << " and |eta|<"<< EFF_MAX_ETA << ":\t" << c_fake_id*1.0/c_fake_reco   << " +/- " << sqrt(c_fake_id)*1.0/c_fake_reco  << std::endl;
  }

  std::cout << "Output written to " << ofilename << std::endl;
  f->Close();

  if ( c_fake_id*1.0/c_fake_reco < WP ) return 0;
  else return 1;

}

int input::is_tau_id(unsigned it, double m_isocut){

  if ( floats_tau_chargedIsoPtSum_nt_obj.at(it) > m_isocut ) return 0;
  if ( floats_tau_decayModeFinding_nt_obj.at(it) < 0.5 ) return 0;

  return 1;
}

void input::get_tau_vec(unsigned it, std::vector<float> &dau_pdgid, std::vector<float> &dau_pt, std::vector<float> &dau_eta, std::vector<float> &dau_phi, std::vector<float> &dau_m){
  dau_pdgid.clear();
  dau_pt.clear();
  dau_eta.clear();
  dau_phi.clear();
  dau_m.clear();

  if ( floats_gen_daughter0pt_nt_obj.at(it) >0 ){
    dau_pdgid.push_back( floats_gen_daughter0pdgId_nt_obj.at(it) );
    dau_pt.push_back( floats_gen_daughter0pt_nt_obj.at(it) );
    dau_eta.push_back( floats_gen_daughter0eta_nt_obj.at(it) );
    dau_phi.push_back( floats_gen_daughter0phi_nt_obj.at(it) );
    dau_m.push_back( floats_gen_daughter0m_nt_obj.at(it) );
  }
  if ( floats_gen_daughter1pt_nt_obj.at(it) >0 ){
    dau_pdgid.push_back( floats_gen_daughter1pdgId_nt_obj.at(it) );
    dau_pt.push_back( floats_gen_daughter1pt_nt_obj.at(it) );
    dau_eta.push_back( floats_gen_daughter1eta_nt_obj.at(it) );
    dau_phi.push_back( floats_gen_daughter1phi_nt_obj.at(it) );
    dau_m.push_back( floats_gen_daughter1m_nt_obj.at(it) );
  }
  if ( floats_gen_daughter2pt_nt_obj.at(it) >0 ){
    dau_pdgid.push_back( floats_gen_daughter2pdgId_nt_obj.at(it) );
    dau_pt.push_back( floats_gen_daughter2pt_nt_obj.at(it) );
    dau_eta.push_back( floats_gen_daughter2eta_nt_obj.at(it) );
    dau_phi.push_back( floats_gen_daughter2phi_nt_obj.at(it) );
    dau_m.push_back( floats_gen_daughter2m_nt_obj.at(it) );
  }
  if ( floats_gen_daughter3pt_nt_obj.at(it) >0 ){
    dau_pdgid.push_back( floats_gen_daughter3pdgId_nt_obj.at(it) );
    dau_pt.push_back( floats_gen_daughter3pt_nt_obj.at(it) );
    dau_eta.push_back( floats_gen_daughter3eta_nt_obj.at(it) );
    dau_phi.push_back( floats_gen_daughter3phi_nt_obj.at(it) );
    dau_m.push_back( floats_gen_daughter3m_nt_obj.at(it) );
  }
  if ( floats_gen_daughter4pt_nt_obj.at(it) >0 ){
    dau_pdgid.push_back( floats_gen_daughter4pdgId_nt_obj.at(it) );
    dau_pt.push_back( floats_gen_daughter4pt_nt_obj.at(it) );
    dau_eta.push_back( floats_gen_daughter4eta_nt_obj.at(it) );
    dau_phi.push_back( floats_gen_daughter4phi_nt_obj.at(it) );
    dau_m.push_back( floats_gen_daughter4m_nt_obj.at(it) );
  }
  if ( floats_gen_daughter5pt_nt_obj.at(it) >0 ){
    dau_pdgid.push_back( floats_gen_daughter5pdgId_nt_obj.at(it) );
    dau_pt.push_back( floats_gen_daughter5pt_nt_obj.at(it) );
    dau_eta.push_back( floats_gen_daughter5eta_nt_obj.at(it) );
    dau_phi.push_back( floats_gen_daughter5phi_nt_obj.at(it) );
    dau_m.push_back( floats_gen_daughter5m_nt_obj.at(it) );
  }
  // if ( floats_gen_daughter6pt_nt_obj.at(it) >0 ){
  //   dau_pdgid.push_back( floats_gen_daughter6pdgId_nt_obj.at(it) );
  //   dau_pt.push_back( floats_gen_daughter6pt_nt_obj.at(it) );
  //   dau_eta.push_back( floats_gen_daughter6eta_nt_obj.at(it) );
  //   dau_phi.push_back( floats_gen_daughter6phi_nt_obj.at(it) );
  //   dau_m.push_back( floats_gen_daughter6m_nt_obj.at(it) );
  // }
  // if ( floats_gen_daughter7pt_nt_obj.at(it) >0 ){
  //   dau_pdgid.push_back( floats_gen_daughter7pdgId_nt_obj.at(it) );
  //   dau_pt.push_back( floats_gen_daughter7pt_nt_obj.at(it) );
  //   dau_eta.push_back( floats_gen_daughter7eta_nt_obj.at(it) );
  //   dau_phi.push_back( floats_gen_daughter7phi_nt_obj.at(it) );
  //   dau_m.push_back( floats_gen_daughter7m_nt_obj.at(it) );
  // }
  // if ( floats_gen_daughter8pt_nt_obj.at(it) >0 ){
  //   dau_pdgid.push_back( floats_gen_daughter8pdgId_nt_obj.at(it) );
  //   dau_pt.push_back( floats_gen_daughter8pt_nt_obj.at(it) );
  //   dau_eta.push_back( floats_gen_daughter8eta_nt_obj.at(it) );
  //   dau_phi.push_back( floats_gen_daughter8phi_nt_obj.at(it) );
  //   dau_m.push_back( floats_gen_daughter8m_nt_obj.at(it) );
  // }
  // if ( floats_gen_daughter9pt_nt_obj.at(it) >0 ){
  //   dau_pdgid.push_back( floats_gen_daughter9pdgId_nt_obj.at(it) );
  //   dau_pt.push_back( floats_gen_daughter9pt_nt_obj.at(it) );
  //   dau_eta.push_back( floats_gen_daughter9eta_nt_obj.at(it) );
  //   dau_phi.push_back( floats_gen_daughter9phi_nt_obj.at(it) );
  //   dau_m.push_back( floats_gen_daughter9m_nt_obj.at(it) );
  // }
  // if ( floats_gen_daughter10pt_nt_obj.at(it) >0 ){
  //   dau_pdgid.push_back( floats_gen_daughter10pdgId_nt_obj.at(it) );
  //   dau_pt.push_back( floats_gen_daughter10pt_nt_obj.at(it) );
  //   dau_eta.push_back( floats_gen_daughter10eta_nt_obj.at(it) );
  //   dau_phi.push_back( floats_gen_daughter10phi_nt_obj.at(it) );
  //   dau_m.push_back( floats_gen_daughter10m_nt_obj.at(it) );
  // }

}

double input::get_local_density(double val){

  if ( !f_pu ){
    int npu=-1;
    if ( this->ofilename.Contains("pu200") ) npu=200;
    else if ( this->ofilename.Contains("pu140") ) npu=140;
    else if ( this->ofilename.Contains("pu35") ) npu=35;
    else if ( this->ofilename.Contains("pu0") ) npu=0;
    else{ npu=200; std::cout << "WARNING: Assuming no pu=200" << std::endl; }

    f_pu = new TF1("fb","gaus(0)");
    f_pu->SetParameters(npu*1/(BS_SIZE*sqrt(2*TMath::Pi())),0,BS_SIZE);
  }

  return f_pu->Eval(val);
}

input::input(TString fname, TTree *tree) : fChain(0)
{
  std::cout << "Trying to open file " << fname << std::endl;
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname);
    if (!f || !f->IsOpen()) {
      f = new TFile(fname);
    }
    if ( !f->IsOpen() ){
      std::cout << "File does not exist" << std::endl;
      return;
    }
    this->ifilename=fname;
    this->ofilename=fname;
    //this->ofilename.Replace(0,fname.Last('/')+1,"");
    //this->ofilename.ReplaceAll("ntuple_","../ntuple/step3_");
    //    this->ofilename.ReplaceAll("v3","v3b");                                                                                                                                                                    
    f->GetObject("Events",tree);
  }
  gErrorIgnoreLevel = kBreak;
  Init(tree);
  gErrorIgnoreLevel = kPrint;
}
