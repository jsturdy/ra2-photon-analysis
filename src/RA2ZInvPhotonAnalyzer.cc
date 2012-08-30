// -*- C++ -*-
//
// Package:    RA2ZInvPhotonAnalyzer
// Class:      RA2ZInvPhotonAnalyzer
// 
/**\class RA2ZInvPhotonAnalyzer RA2ZInvPhotonAnalyzer.cc SusyAnalysis/RA2ZInvPhotonAnalyzer/src/RA2ZInvPhotonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seema Sharma
//         Created:  Mon Jun 20 12:58:08 CDT 2011
// $Id: RA2ZInvPhotonAnalyzer.cc,v 1.1 2012/07/20 11:35:34 sturdy Exp $
//
//


// system include files
#include <memory>
#include <iomanip>
#include <cmath>

#include "ZInvisibleBkgds/Photons/interface/RA2ZInvPhotonAnalyzer.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include <DataFormats/METReco/interface/MET.h>


RA2ZInvPhotonAnalyzer::RA2ZInvPhotonAnalyzer(const edm::ParameterSet& pset) {

  // read parameters from config file
  debug_          = pset.getParameter<bool>("Debug");
  doGenAnalysis_  = pset.getParameter<bool>("DoGenAnalysis");
  genParticles_   = pset.getParameter<edm::InputTag>("GenParticles");  
  data_           = pset.getParameter<bool>("Data");
  photonSrc_      = pset.getParameter<edm::InputTag>("PhotonSrc");
  photonIsoSrc_   = pset.getParameter<edm::InputTag>("PhotonIsoSrc");
  jetSrc_         = pset.getParameter<edm::InputTag>("JetSrc");
  jetNoPhotSrc_   = pset.getParameter<edm::InputTag>("JetNoPhotSrc");
  jetHTSrc_       = pset.getParameter<edm::InputTag>("JetHTSource");
  doPUReWeight_   = pset.getParameter<bool>("DoPUReweight");
  puWeightSrc_    = pset.getParameter<edm::InputTag>("PUWeightSource");
  pfRhoSrc_       = pset.getParameter<edm::InputTag>("PFRhoSource");
  
  //minPhoPt_   = pset.getParameter<double>("minPhoPt");
  //ebMaxEta_   = pset.getParameter<double>("ebMaxEta");
  //eeMinEta_   = pset.getParameter<double>("eeMinEta");
  //eeMaxEta_   = pset.getParameter<double>("eeMaxEta");
  //ebSigIeIe_  = pset.getParameter<double>("ebSigIeIe");
  //eeSigIeIe_  = pset.getParameter<double>("eeSigIeIe");
  //recoHoverE_ = pset.getParameter<double>("recoHoverE");

  ra2NJets_         = pset.getParameter<unsigned int>("RA2NJets");
  ra2HT_            = pset.getParameter<double>("RA2HT");
  ra2MHT_           = pset.getParameter<double>("RA2MHT");
  ra2ApplyDphiCuts_ = pset.getParameter<bool>("RA2ApplyDphiCuts");
}

RA2ZInvPhotonAnalyzer::~RA2ZInvPhotonAnalyzer() {

}

void RA2ZInvPhotonAnalyzer::analyze(const edm::Event& ev, const edm::EventSetup& es) {

  using namespace edm;

  // get event-id
  unsigned int event = (ev.id()).event();
  unsigned int run   = (ev.id()).run();
  unsigned int lumi  =  ev.luminosityBlock();

  // reject event if there an isolated electron or muon present
  edm::Handle< std::vector<pat::Electron> > patElectrons;
  ev.getByLabel("patElectronsIDIso", patElectrons);

  edm::Handle< std::vector<pat::Muon> > patMuons;
  ev.getByLabel("patMuonsIDIso", patMuons);
  
  if (patElectrons->size() != 0 || patMuons->size() != 0) { 
    std::cout << "Isolated Lepton found : Event Rejected : ( run, event, lumi ) " 
    << run << " " << event << " " << lumi << std::endl;
    return;
  }
  
  // get vertices
  edm::Handle< std::vector<reco::Vertex> > Vertices;
  ev.getByLabel("goodVertices", Vertices);
  int nVertices = Vertices->size();
  
  // get photons 
  edm::Handle< std::vector<pat::Photon> > patPhotons;
  ev.getByLabel(photonSrc_, patPhotons); 

  edm::Handle< edm::View<pat::Photon> > patPhotonsIso;
  ev.getByLabel(photonIsoSrc_, patPhotonsIso); 

  // get jetcollection
  edm::Handle< std::vector<pat::Jet> > recoJets;
  ev.getByLabel(jetSrc_, recoJets); 
  
  edm::Handle<edm::View<pat::Jet> > jets;
  ev.getByLabel(jetSrc_, jets);

  edm::Handle<edm::View<pat::Jet> > jetsNoPhots;
  ev.getByLabel(jetNoPhotSrc_, jetsNoPhots);

  edm::Handle<edm::View<pat::Jet> > jetsHT;
  ev.getByLabel(jetHTSrc_, jetsHT);

  // if MC, do PU reweighting
  double pu_event_wt = 1.0;
  edm::Handle<double> puweight;
  if( doPUReWeight_ ) {
    ev.getByLabel(puWeightSrc_, puweight);
    pu_event_wt = *puweight;
  }
  h_puWeight ->Fill(pu_event_wt);
  h_Vertices ->Fill(nVertices);
  h_VerticesReWeighted ->Fill(nVertices, pu_event_wt);

  // get pfRho for pileup correction to isolation
  edm::Handle<double> pfRho;
  ev.getByLabel(pfRhoSrc_, pfRho);
  double pf_event_rho = *pfRho;
  
  if( nVertices<60 )  h_pfRho[nVertices]->Fill(pf_event_rho);
  
  std::vector<const pat::Photon*> IsoPhotons;
  std::vector<const pat::Photon*> RA2Photons;
  double photon_he     = 0.05;
  double photon_ptcut  = 50.0;
  double photon_etacut = 2.50;
  double sigma_barrel  = 0.011;
  double sigma_endcap  = 0.030;
  double iso_comb      = 5.0; // combined isolation
  double areaR03       = 0.3 * 0.3 * 4.0*atan(1.0);
  double areaR04       = 0.4 * 0.4 * 4.0*atan(1.0);
  
  //std::cout << "areaR04 " << areaR04 << std::endl;
  h_PhotonsCand_NPhotons->Fill(patPhotons->size(), pu_event_wt); 
  for(unsigned int iphot=0; iphot<patPhotons->size(); ++iphot) {
    
    const pat::Photon *p1 = &((*patPhotons)[iphot]);
    
    double pt           = p1->pt();
    double et           = p1->et();
    double eta          = p1->eta();
    double phi          = p1->phi();
    double scEta        = fabs(p1->superCluster()->eta());
    double sigEtaEta    = p1->sigmaIetaIeta();
    double sigPhiPhi    = 0.0;
    double trackerIso   = p1->trkSumPtHollowConeDR03();
    double ecalIso      = p1->ecalRecHitSumEtConeDR03();
    double hcalIso      = p1->hcalTowerSumEtConeDR03();
    double hcalIso2012  = p1->userFloat("hcalIsoConeDR03_2012");

    double hadOverEm     = p1->hadronicOverEm();
    double hadOverEm2012 = p1->hadTowOverEm();

    double pfChargedIso = p1->userIsolation(pat::IsolationKeys(pat::UserBaseIso+0));
    double pfNeutralIso = p1->userIsolation(pat::IsolationKeys(pat::UserBaseIso+2));
    double pfGammaIso   = p1->userIsolation(pat::IsolationKeys(pat::UserBaseIso+3));

    double pfChargedIsoPU = p1->userFloat("pfChargedPU");
    double pfNeutralIsoPU = p1->userFloat("pfNeutralPU");
    double pfGammaIsoPU   = p1->userFloat("pfGammaPU");

    double pfChargedIsoRel = p1->userFloat("pfChargedRel");
    double pfNeutralIsoRel = p1->userFloat("pfNeutralRel");
    double pfGammaIsoRel   = p1->userFloat("pfGammaRel");
    
    double pfChargedIsoPURel = p1->userFloat("pfChargedPURel");
    double pfNeutralIsoPURel = p1->userFloat("pfNeutralPURel");
    double pfGammaIsoPURel   = p1->userFloat("pfGammaPURel");
    
    bool   isPixel     = p1->hasPixelSeed();
    bool   passElecVeto = p1->userInt("passElectronConvVeto");
    bool   hOverE      = ( hadOverEm < photon_he);
    //bool   hOverE2012  = ( hadOverEm2012 < photon_he);
    bool   hOverE2012  = ( hadOverEm2012 < p1->userFloat("hadTowOverEmTightCut"));
    //bool   kineAcc     = ( p1->et()>photon_ptcut && ((p1->isEE() && std::fabs(p1->eta())<2.5) || p1->isEB()) ) ;
    bool   kineAcc     = ( p1->et() > photon_ptcut && 
			 ((p1->isEE() && std::fabs(eta) > 1.566 && std::fabs(eta) < 2.5) || 
			  (p1->isEB() && std::fabs(eta) < 1.4442) ) ) ;
    //bool   showerShape = ( (std::fabs(eta)<1.4442 && sigEtaEta<sigma_barrel) || (std::fabs(eta)>1.566 && sigEtaEta<sigma_endcap) );
    bool   showerShape = sigEtaEta < p1->userFloat("showerShapeTightCut");

    bool   isPhotonIso       = ( (trackerIso+ecalIso+hcalIso     - areaR04*pf_event_rho)<5.0 );
    bool   isPhotonIso2012   = ( (trackerIso+ecalIso+hcalIso2012 - areaR04*pf_event_rho)<5.0 );
    
    bool   isPhotonPFChargedIso = (pfChargedIsoPU < p1->userFloat("pfChargedTightCut"));
    bool   isPhotonPFNeutralIso = (pfNeutralIsoPU < p1->userFloat("pfNeutralTightCut"));
    bool   isPhotonPFGammaIso   = (pfGammaIsoPU   < p1->userFloat("pfGammaTightCut"));
    bool   isPhotonPFChargedIsoRel = (pfChargedIsoPURel < p1->userFloat("pfChargedRelTightCut"));
    bool   isPhotonPFNeutralIsoRel = (pfNeutralIsoPURel < p1->userFloat("pfNeutralRelTightCut"));
    bool   isPhotonPFGammaIsoRel   = (pfGammaIsoPURel   < p1->userFloat("pfGammaRelTightCut"));

    bool   isPhotonPFIso = ( isPhotonPFChargedIso && isPhotonPFNeutralIso && isPhotonPFGammaIso );

    bool   isPhoton2011    = ( kineAcc && !isPixel && hOverE     && showerShape && isPhotonIso );
    bool   isPhoton2012    = ( kineAcc && !isPixel && hOverE2012 && showerShape && isPhotonIso2012 );
    //bool   isPhoton2012PF  = ( kineAcc && !isPixel && hOverE2012 && showerShape && isPhotonPFIso );
    bool   isPhoton2012PF  = ( kineAcc && passElecVeto && hOverE2012 && showerShape && isPhotonPFIso );

   
    if( iphot<2 ) {
      h_PhotonsCand_Pt [iphot]      ->Fill(pt,  pu_event_wt);
      h_PhotonsCand_Eta[iphot]      ->Fill(eta, pu_event_wt);
      h_PhotonsCand_Phi[iphot]      ->Fill(phi, pu_event_wt);

      h_PhotonsCand_TrkIso[iphot]     ->Fill(trackerIso,  pu_event_wt);
      h_PhotonsCand_EcalIso[iphot]    ->Fill(ecalIso,     pu_event_wt);
      h_PhotonsCand_HcalIso[iphot]    ->Fill(hcalIso,     pu_event_wt);
      h_PhotonsCand_HcalIso2012[iphot]->Fill(hcalIso2012, pu_event_wt);
      
      h_PhotonsCand_PFCharged[iphot] ->Fill(pfChargedIso, pu_event_wt);
      h_PhotonsCand_PFNeutral[iphot] ->Fill(pfNeutralIso, pu_event_wt);
      h_PhotonsCand_PFGamma[iphot]   ->Fill(pfGammaIso,   pu_event_wt);
      
      h_PhotonsCand_PFChargedRel[iphot] ->Fill(pfChargedIsoRel, pu_event_wt);
      h_PhotonsCand_PFNeutralRel[iphot] ->Fill(pfNeutralIsoRel, pu_event_wt);
      h_PhotonsCand_PFGammaRel[iphot]   ->Fill(pfGammaIsoRel,   pu_event_wt);
      
      h_PhotonsCand_PFChargedRelPU[iphot] ->Fill(pfChargedIsoPURel, pu_event_wt);
      h_PhotonsCand_PFNeutralRelPU[iphot] ->Fill(pfNeutralIsoPURel, pu_event_wt);
      h_PhotonsCand_PFGammaRelPU[iphot]   ->Fill(pfGammaIsoPURel,   pu_event_wt);
      
      h_PhotonsCand_hOverE[iphot]    ->Fill(hadOverEm,     pu_event_wt);
      h_PhotonsCand_hOverE2012[iphot]->Fill(hadOverEm2012, pu_event_wt);
      
      if(std::fabs(eta)<1.4442)     h_PhotonsCand_SigEtaEta_EB[iphot]->Fill(sigEtaEta, pu_event_wt);
      else if(std::fabs(eta)>1.566) h_PhotonsCand_SigEtaEta_EE[iphot]->Fill(sigEtaEta, pu_event_wt);
    }

    if( isPhoton2012PF ) IsoPhotons.push_back(p1);
  } // loop over patPhotons

  //compare the sizes of the collections computed within the code an with the external module
  if (IsoPhotons.size() != patPhotonsIso->size()) {
    std::cout<<"isolated photon collection sizes do not match!"<<std::endl;
    std::cout<<"IsoPhotons::"<<IsoPhotons.size()<<"  --  patPhotonsIso::"<<patPhotonsIso->size()<<std::endl;

    if (debug_) {
      bool isPhoton2012PF_tmp(false), 
	kineAcc_tmp(false), 
	isPixel_tmp(false), 
	passElecVeto_tmp(false), 
	hOverE_tmp(false), 
	hOverE2012_tmp(false), 
	showerShape_tmp(false),
	isPhotonPFChargedIso_tmp(false),
	isPhotonPFNeutralIso_tmp(false),
	isPhotonPFGammaIso_tmp(false),
	isPhotonPFIso_tmp(false);
    
      isPixel_tmp      = (*patPhotons)[0].hasPixelSeed();
      passElecVeto_tmp = (*patPhotons)[0].userInt("passElectronConvVeto");
      hOverE_tmp      = ( (*patPhotons)[0].hadronicOverEm() < 0.5);
      hOverE2012_tmp  = ( (*patPhotons)[0].hadTowOverEm() < (*patPhotons)[0].userFloat("hadTowOverEmTightCut"));
      kineAcc_tmp     = ( (*patPhotons)[0].et() > photon_ptcut && 
			  (((*patPhotons)[0].isEE() && std::fabs((*patPhotons)[0].eta()) > 1.566 && std::fabs((*patPhotons)[0].eta()) < 2.5) || 
			   ((*patPhotons)[0].isEB() && std::fabs((*patPhotons)[0].eta()) < 1.4442) ) ) ;
      showerShape_tmp = (*patPhotons)[0].sigmaIetaIeta() < (*patPhotons)[0].userFloat("showerShapeTightCut");
    
      isPhotonPFChargedIso_tmp = ((*patPhotons)[0].userFloat("pfChargedPURel") < (*patPhotons)[0].userFloat("pfChargedTightCut"));
      isPhotonPFNeutralIso_tmp = ((*patPhotons)[0].userFloat("pfNeutralPURel") < (*patPhotons)[0].userFloat("pfNeutralTightCut"));
      isPhotonPFGammaIso_tmp   = ((*patPhotons)[0].userFloat("pfGammaPURel")   < (*patPhotons)[0].userFloat("pfGammaTightCut"));
    
      isPhotonPFIso_tmp = ( isPhotonPFChargedIso_tmp && isPhotonPFNeutralIso_tmp && isPhotonPFGammaIso_tmp );
    
      //isPhoton2012PF_tmp  = ( kineAcc_tmp && !isPixel_tmp && hOverE2012_tmp && showerShape_tmp && isPhotonPFIso_tmp );
      isPhoton2012PF_tmp  = ( kineAcc_tmp && passElecVeto_tmp && hOverE2012_tmp && showerShape_tmp && isPhotonPFIso_tmp );
    
      if (patPhotonsIso->size() > 0 && debug_) {
	std::cout<<
	  //"isPhoton2012PF pass      ("<<isPhoton2012PF_tmp      <<")  -  "<< 
	  "kineAcc pass("     <<kineAcc_tmp             <<")  -  "<< 
	  "!isPixel pass("    <<!isPixel_tmp            <<")  -  "<< 
	  "passElecVeto pass("<<passElecVeto_tmp        <<")  -  "<< 
	  "hOverE pass("      <<hOverE_tmp              <<")  -  "<< 
	  "hOverE2012 pass("  <<hOverE2012_tmp          <<")  -  "<< 
	  "showerShape pass(" <<showerShape_tmp         <<")  -  "<<
	  "pfChargedIso pass("<<isPhotonPFChargedIso_tmp<<")  -  "<<
	  "pfNeutralIso pass("<<isPhotonPFNeutralIso_tmp<<")  -  "<<
	  "pfGammaIso pass("  <<isPhotonPFGammaIso_tmp  <<")"<<//  -  "<<
	  //"isPhotonPFIso pass       ("<<isPhotonPFIso_tmp       <<")  -  "<<
	  std::endl;
      
	std::cout<<"          ::   pt      eta    pixel    hOverE    hOverE2012   sieie    pfCHRel    pfNURel   pfGARel    pfCHPURel    pfNUPURel   pfGAPURel    pfCHPUSub    pfNUPUSub   pfGAPUSub"<<std::endl;
	std::cout<<"tight cuts::"<<
	  20.0<<"    "<<
	  2.5 <<"    "<<
	  0   <<"    "<<
	  0.5 <<"    "<<
	  (*patPhotons)[0].userFloat("hadTowOverEmTightCut")<<"    "<<
	  (*patPhotons)[0].userFloat("showerShapeTightCut")<<"    "<<
	  (*patPhotons)[0].userFloat("pfChargedTightCut")<<"    "<<
	  (*patPhotons)[0].userFloat("pfNeutralTightCut")<<"    "<<
	  (*patPhotons)[0].userFloat("pfGammaTightCut")<<"    "<<
	  std::endl;
	std::cout<<"patPhotons::"<<
	  (*patPhotons)[0].pt()<<"    "<<
	  (*patPhotons)[0].eta()<<"    "<<
	  (*patPhotons)[0].hasPixelSeed()<<"    "<<
	  (*patPhotons)[0].hadronicOverEm()<<"    "<<
	  (*patPhotons)[0].hadTowOverEm()<<"    "<<
	  (*patPhotons)[0].sigmaIetaIeta()<<"    "<<
	  (*patPhotons)[0].userFloat("pfChargedRel")<<"    "<<
	  (*patPhotons)[0].userFloat("pfNeutralRel")<<"    "<<
	  (*patPhotons)[0].userFloat("pfGammaRel")<<"    "<<
	  (*patPhotons)[0].userFloat("pfChargedPURel")<<"    "<<
	  (*patPhotons)[0].userFloat("pfNeutralPURel")<<"    "<<
	  (*patPhotons)[0].userFloat("pfGammaPURel")<<"    "<<
	  (*patPhotons)[0].userFloat("pfChargedPUSub")<<"    "<<
	  (*patPhotons)[0].userFloat("pfNeutralPUSub")<<"    "<<
	  (*patPhotons)[0].userFloat("pfGammaPUSub")<<"    "<<
	  std::endl;
	std::cout<<"patPhotonsIso::"<<
	  (*patPhotonsIso)[0].pt()<<"    "<<
	  (*patPhotonsIso)[0].eta()<<"    "<<
	  (*patPhotonsIso)[0].hasPixelSeed()<<"    "<<
	  (*patPhotonsIso)[0].hadronicOverEm()<<"    "<<
	  (*patPhotonsIso)[0].hadTowOverEm()<<"    "<<
	  (*patPhotonsIso)[0].sigmaIetaIeta()<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfChargedRel")<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfNeutralRel")<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfGammaRel")<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfChargedPURel")<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfNeutralPURel")<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfGammaPURel")<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfChargedPUSub")<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfNeutralPUSub")<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfGammaPUSub")<<"    "<<
	  std::endl;
      }
    }
  }
  // do the gen analysis if MC events
  //if(! ev.isRealData()) {
  //  edm::Handle<reco::GenParticleCollection> genParticles;
  //  ev.getByLabel(genParticles_, genParticles);
  //  doGenAnalysis(genParticles, jetsHT, jets, IsoPhotons, pu_event_wt);
  //}

  // require atleast one isolated photon else return
  if (IsoPhotons.size()<1) return;

  h_PhotonsIso_NPhotons->Fill(IsoPhotons.size(), pu_event_wt); 
  for( unsigned int iphot=0; iphot<IsoPhotons.size(); iphot++) {
    if( iphot<2 ) {
      //double hcalIso2012Spec = IsoPhotons[iphot]->hcalTowerSumEtConeDR04() + (IsoPhotons[iphot]->hadronicOverEm() - IsoPhotons[iphot]->hadTowOverEm())*IsoPhotons[iphot]->superCluster->energy()/cosh(IsoPhotons[iphot]->superCluster()->Eta();

      h_PhotonsIso_Pt [iphot]       ->Fill(IsoPhotons[iphot]->pt(),  pu_event_wt);
      h_PhotonsIso_Eta[iphot]       ->Fill(IsoPhotons[iphot]->eta(), pu_event_wt);
      h_PhotonsIso_Phi[iphot]       ->Fill(IsoPhotons[iphot]->phi(), pu_event_wt);

      h_PhotonsIso_TrkIso[iphot]     ->Fill(IsoPhotons[iphot]->trkSumPtHollowConeDR03(),  pu_event_wt);
      h_PhotonsIso_EcalIso[iphot]    ->Fill(IsoPhotons[iphot]->ecalRecHitSumEtConeDR03(), pu_event_wt);
      h_PhotonsIso_HcalIso[iphot]    ->Fill(IsoPhotons[iphot]->hcalTowerSumEtConeDR03(),  pu_event_wt);
      h_PhotonsIso_HcalIso2012[iphot]->Fill(IsoPhotons[iphot]->userFloat("hcalIsoConeDR03_2012"), pu_event_wt);

      h_PhotonsIso_PFCharged[iphot] ->Fill(IsoPhotons[iphot]->userIsolation(pat::IsolationKeys(pat::UserBaseIso+0)), pu_event_wt);
      h_PhotonsIso_PFNeutral[iphot] ->Fill(IsoPhotons[iphot]->userIsolation(pat::IsolationKeys(pat::UserBaseIso+2)), pu_event_wt);
      h_PhotonsIso_PFGamma[iphot]   ->Fill(IsoPhotons[iphot]->userIsolation(pat::IsolationKeys(pat::UserBaseIso+3)), pu_event_wt);

      double pfChargedIsoRel = IsoPhotons[iphot]->userFloat("pfChargedRel");
      double pfNeutralIsoRel = IsoPhotons[iphot]->userFloat("pfNeutralRel");
      double pfGammaIsoRel   = IsoPhotons[iphot]->userFloat("pfGammaRel");
      
      double pfChargedIsoPURel = IsoPhotons[iphot]->userFloat("pfChargedPURel");
      double pfNeutralIsoPURel = IsoPhotons[iphot]->userFloat("pfNeutralPURel");
      double pfGammaIsoPURel   = IsoPhotons[iphot]->userFloat("pfGammaPURel");
      
      h_PhotonsIso_PFChargedRel[iphot] ->Fill(pfChargedIsoRel, pu_event_wt);
      h_PhotonsIso_PFNeutralRel[iphot] ->Fill(pfNeutralIsoRel, pu_event_wt);
      h_PhotonsIso_PFGammaRel[iphot]   ->Fill(pfGammaIsoRel,   pu_event_wt);
      
      h_PhotonsIso_PFChargedRelPU[iphot] ->Fill(pfChargedIsoPURel, pu_event_wt);
      h_PhotonsIso_PFNeutralRelPU[iphot] ->Fill(pfNeutralIsoPURel, pu_event_wt);
      h_PhotonsIso_PFGammaRelPU[iphot]   ->Fill(pfGammaIsoPURel,   pu_event_wt);
      
      h_PhotonsIso_hOverE[iphot]    ->Fill(IsoPhotons[iphot]->hadronicOverEm(), pu_event_wt);
      h_PhotonsIso_hOverE2012[iphot]->Fill(IsoPhotons[iphot]->hadTowOverEm(), pu_event_wt);
      h_PhotonsIso_R9[iphot]        ->Fill(IsoPhotons[iphot]->r9(), pu_event_wt);

      double eta = IsoPhotons[iphot]->eta();
      double sigEtaEta = IsoPhotons[iphot]->sigmaIetaIeta();
      if(std::fabs(eta)<1.4442)     h_PhotonsIso_SigEtaEta_EB[iphot]->Fill(sigEtaEta, pu_event_wt);
      else if(std::fabs(eta)>1.566) h_PhotonsIso_SigEtaEta_EE[iphot]->Fill(sigEtaEta, pu_event_wt);
    }
  }

  double photon_pt  = IsoPhotons[0]->pt();
  double photon_eta = IsoPhotons[0]->eta();
  double photon_phi = IsoPhotons[0]->phi();

  double photon2_pt  = -1;
  double photon2_eta = -999.;
  double photon2_phi = -999.;
  if (IsoPhotons.size() > 1) {
    photon2_pt  = IsoPhotons[1]->pt();
    photon2_eta = IsoPhotons[1]->eta();
    photon2_phi = IsoPhotons[1]->phi();
  }
  for(int ipt=0; ipt<NPhotPtBins-1; ipt++) {
    for(int ieta=0; ieta<NPhotEtaBins-1; ieta++) {
      if( photon_pt>PhotPtVal[ipt] && photon_pt<PhotPtVal[ipt+1] ) {
	if( std::fabs(photon_eta)>PhotEtaVal[ieta] && std::fabs(photon_eta)<PhotEtaVal[ieta+1] ) {
	  NPhotonPtEta_Iso[ipt][ieta] += pu_event_wt;
	  
	}
      }
    }
  }
  
  // clean jet collection

  std::vector<const pat::Jet*> Jets; // clean jet collection

  int    bestDRPhot = -1;
  int    bestDRJet  = -1;
  double bestDRMin  =999.0;

  double dRMin=999.0;
  int    dRMinJetIdx = -1;
  double dRMin2=999.0;
  int    dRMinJetIdx2 = -1;
  for(unsigned int i=0; i<recoJets->size(); i++) {
    const pat::Jet *r = &((*recoJets)[i]);
   
    double jet_pt  = r->pt();
    double jet_eta = r->eta();
    double jet_phi = r->phi();
    
    if( jet_pt<30.0 ) continue;
    
    double dR = reco::deltaR(photon_eta, photon_phi, jet_eta, jet_phi);
    if(dR<dRMin) { 
      dRMin       = dR;
      dRMinJetIdx = i;
    }

    if (IsoPhotons.size() > 1) {
      double dR2 = reco::deltaR(photon2_eta, photon2_phi, jet_eta, jet_phi);
      if(dR2<dRMin2) { 
	dRMin2       = dR2;
	dRMinJetIdx2 = i;
      }
    }
  }
  
  //original
  //if(dRMinJetIdx>=0) h_drMin_LeadPhotPFJet->Fill( dRMin );
  if (dRMin < dRMin2) {
    bestDRMin = dRMin;
    bestDRPhot = 0;
    bestDRJet = dRMinJetIdx;
  }
  else if (dRMin2 < dRMin) {
    bestDRMin = dRMin2;
    bestDRPhot = 1;
    bestDRJet = dRMinJetIdx2;
  }
  if(bestDRJet>=0) h_drMin_LeadPhotPFJet->Fill( bestDRMin );


  // if no matching jet found in dRMin<0.1 then reject this event
  //if(dRMinJetIdx<0 || dRMin>0.1 ) return;
  if(bestDRJet<0 || bestDRMin>0.1 ) return;

  if ( bestDRPhot == 0 ) {
    RA2Photons.push_back(IsoPhotons[0]);
    if (IsoPhotons.size() > 1)
      RA2Photons.push_back(IsoPhotons[1]);
  }
  if ( bestDRPhot == 1 ) {
    RA2Photons.push_back(IsoPhotons[1]);
    RA2Photons.push_back(IsoPhotons[0]);
    }

  for(unsigned int i=0; i<recoJets->size(); ++i) {
    const pat::Jet *r = &((*recoJets)[i]);
    int jj = (int) i;
    if( jj != bestDRJet ) Jets.push_back(r);
  }

  if(bestDRJet<0) {
    std::cout << "Isolated photons : " << std::endl;
    for( unsigned int iphot=0; iphot<IsoPhotons.size(); iphot++) {
      std::cout << iphot << " " <<IsoPhotons[iphot]->pt() << " " << IsoPhotons[iphot]->eta()
		<< IsoPhotons[iphot]->phi() 
		<< std::endl;
    }
    std::cout << "Jets : " << std::endl;
    for(unsigned int i=0; i<recoJets->size(); i++) {
      const pat::Jet *r = &((*recoJets)[i]);
      std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
    }
    std::cout << "Cleaned Jets : " << std::endl;
    for(unsigned int i=0; i<Jets.size(); ++i) {
      const pat::Jet *r = Jets[i];
      std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
    }
  }
  
  //compare the sizes of the collections computed within the code an with the external module
  if (Jets.size() != jetsNoPhots->size()) {
    std::cout<<"cleaned jet collection sizes do not match!"<<std::endl;
    std::cout<<"Jets::"<<Jets.size()<<"  --  jetsNoPhots::"<<jetsNoPhots->size()<<std::endl;
    std::cout<<"IsoPhotons::"<<IsoPhotons.size()<<"  --  patPhotonsIso::"<<patPhotonsIso->size()<<std::endl;
    if (debug_) {
      bool isPhoton2012PF_tmp(false), 
	kineAcc_tmp(false), 
	isPixel_tmp(false), 
	hOverE_tmp(false), 
	hOverE2012_tmp(false), 
	showerShape_tmp(false),
	isPhotonPFChargedIso_tmp(false),
	isPhotonPFNeutralIso_tmp(false),
	isPhotonPFGammaIso_tmp(false),
	isPhotonPFIso_tmp(false);
      
      isPixel_tmp     = (*patPhotons)[1].hasPixelSeed();
      hOverE_tmp      = ((*patPhotons)[1].hadronicOverEm() < 0.5);
      hOverE2012_tmp  = ((*patPhotons)[1].hadTowOverEm() < (*patPhotons)[1].userFloat("hadTowOverEmTightCut"));
      kineAcc_tmp     = ((*patPhotons)[1].et() > photon_ptcut && 
			 (((*patPhotons)[1].isEE() && std::fabs((*patPhotons)[1].eta()) > 1.566 && std::fabs((*patPhotons)[1].eta()) < 2.5) || 
			  ((*patPhotons)[1].isEB() && std::fabs((*patPhotons)[1].eta()) < 1.4442) ) ) ;
      showerShape_tmp = (*patPhotons)[1].sigmaIetaIeta() < (*patPhotons)[1].userFloat("showerShapeTightCut");
      
      isPhotonPFChargedIso_tmp = ((*patPhotons)[1].userFloat("pfChargedPURel") < (*patPhotons)[1].userFloat("pfChargedTightCut"));
      isPhotonPFNeutralIso_tmp = ((*patPhotons)[1].userFloat("pfNeutralPURel") < (*patPhotons)[1].userFloat("pfNeutralTightCut"));
      isPhotonPFGammaIso_tmp   = ((*patPhotons)[1].userFloat("pfGammaPURel")   < (*patPhotons)[1].userFloat("pfGammaTightCut"));
      
      isPhotonPFIso_tmp = ( isPhotonPFChargedIso_tmp && isPhotonPFNeutralIso_tmp && isPhotonPFGammaIso_tmp );
      
      isPhoton2012PF_tmp  = ( kineAcc_tmp && !isPixel_tmp && hOverE2012_tmp && showerShape_tmp && isPhotonPFIso_tmp );
      
      if (patPhotonsIso->size() > 1) {
	std::cout<<
	  "isPhoton2012PF pass      ("<<isPhoton2012PF_tmp      <<")  -  "<< 
	  "kineAcc pass("     <<kineAcc_tmp             <<")  -  "<< 
	  "!isPixel pass("    <<!isPixel_tmp            <<")  -  "<< 
	  "hOverE pass("      <<hOverE_tmp              <<")  -  "<< 
	  "hOverE2012 pass("  <<hOverE2012_tmp          <<")  -  "<< 
	  "showerShape pass(" <<showerShape_tmp         <<")  -  "<<
	  "pfChargedIso pass("<<isPhotonPFChargedIso_tmp<<")  -  "<<
	  "pfNeutralIso pass("<<isPhotonPFNeutralIso_tmp<<")  -  "<<
	  "pfGammaIso pass("  <<isPhotonPFGammaIso_tmp  <<")"<<//  -  "<<
	  //"isPhotonPFIso pass       ("<<isPhotonPFIso_tmp       <<")  -  "<<
	  std::endl;
      
	std::cout<<"lead      ::index    dRMin    pt(raw)    eta(raw)"<<std::endl;
	std::cout<<dRMinJetIdx<<"    "<<dRMin<<"    "<<
	  (*jets)[dRMinJetIdx].pt()<<"("<<
	  (*jets)[dRMinJetIdx].correctedJet("Uncorrected").pt()<<")    "<<
	  (*jets)[dRMinJetIdx].eta()<<"("<<
	  (*jets)[dRMinJetIdx].correctedJet("Uncorrected").eta()<<")    "<<
	  std::endl;
	std::cout<<"second    ::index2    dRMin2    pt(raw)    eta(raw)"<<std::endl;
	std::cout<<dRMinJetIdx2<<"    "<<dRMin2<<"    "<<
	  (*jets)[dRMinJetIdx2].pt()<<"("<<
	  (*jets)[dRMinJetIdx2].correctedJet("Uncorrected").pt()<<")    "<<
	  (*jets)[dRMinJetIdx2].eta()<<"("<<
	  (*jets)[dRMinJetIdx2].correctedJet("Uncorrected").eta()<<")    "<<
	  std::endl;
	std::cout<<"          ::   pt      eta    pixel    hOverE    hOverE2012   sieie    pfCHPURel    pfNUPURel   pfGAPURel"<<std::endl;
	std::cout<<"tight cuts::"<<
	  50.0<<"    "<<
	  2.5 <<"    "<<
	  0   <<"    "<<
	  0.5 <<"    "<<
	  (*patPhotons)[1].userFloat("hadTowOverEmTightCut")<<"    "<<
	  (*patPhotons)[1].userFloat("showerShapeTightCut")<<"    "<<
	  (*patPhotons)[1].userFloat("pfChargedTightCut")<<"    "<<
	  (*patPhotons)[1].userFloat("pfNeutralTightCut")<<"    "<<
	  (*patPhotons)[1].userFloat("pfGammaTightCut")<<"    "<<
	  "pdgid    status    mother    status"<<
	  std::endl;
	std::cout<<"patPhotons[0]::"<<
	  (*patPhotons)[0].pt()<<"    "<<
	  (*patPhotons)[0].eta()<<"    "<<
	  (*patPhotons)[0].hasPixelSeed()<<"    "<<
	  (*patPhotons)[0].hadronicOverEm()<<"    "<<
	  (*patPhotons)[0].hadTowOverEm()<<"    "<<
	  (*patPhotons)[0].sigmaIetaIeta()<<"    "<<
	  (*patPhotons)[0].userFloat("pfChargedPURel")<<"    "<<
	  (*patPhotons)[0].userFloat("pfNeutralPURel")<<"    "<<
	  (*patPhotons)[0].userFloat("pfGammaPURel")<<"    ";
	if (!data_) {
	  std::cout<<
	    (*patPhotons)[0].genPhoton()->pdgId()<<"    "<<
	    (*patPhotons)[0].genPhoton()->status()<<"    "<<
	    (*patPhotons)[0].genPhoton()->mother()->pdgId()<<"    "<<
	    (*patPhotons)[0].genPhoton()->mother()->status();
	}
	else
	  std::cout<<"0    0    0    0";
	std::cout<<std::endl;
	std::cout<<"patPhotons[1]::"<<
	  (*patPhotons)[1].pt()<<"    "<<
	  (*patPhotons)[1].eta()<<"    "<<
	  (*patPhotons)[1].hasPixelSeed()<<"    "<<
	  (*patPhotons)[1].hadronicOverEm()<<"    "<<
	  (*patPhotons)[1].hadTowOverEm()<<"    "<<
	  (*patPhotons)[1].sigmaIetaIeta()<<"    "<<
	  (*patPhotons)[1].userFloat("pfChargedPURel")<<"    "<<
	  (*patPhotons)[1].userFloat("pfNeutralPURel")<<"    "<<
	  (*patPhotons)[1].userFloat("pfGammaPURel")<<"    ";
	if (!data_) {
	  std::cout<<
	    (*patPhotons)[1].genPhoton()->pdgId()<<"    "<<
	    (*patPhotons)[1].genPhoton()->status()<<"    "<<
	    (*patPhotons)[1].genPhoton()->mother()->pdgId()<<"    "<<
	    (*patPhotons)[1].genPhoton()->mother()->status();
	}
	else
	  std::cout<<"0    0    0    0";
	std::cout<<std::endl;

	std::cout<<"patPhotonsIso[0]::"<<
	  (*patPhotonsIso)[0].pt()<<"    "<<
	  (*patPhotonsIso)[0].eta()<<"    "<<
	  (*patPhotonsIso)[0].hasPixelSeed()<<"    "<<
	  (*patPhotonsIso)[0].hadronicOverEm()<<"    "<<
	  (*patPhotonsIso)[0].hadTowOverEm()<<"    "<<
	  (*patPhotonsIso)[0].sigmaIetaIeta()<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfChargedPURel")<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfNeutralPURel")<<"    "<<
	  (*patPhotonsIso)[0].userFloat("pfGammaPURel")<<"    ";
	if (!data_) {
	  std::cout<<
	    (*patPhotonsIso)[0].genPhoton()->pdgId()<<"    "<<
	    (*patPhotonsIso)[0].genPhoton()->status()<<"    "<<
	    (*patPhotonsIso)[0].genPhoton()->mother()->pdgId()<<"    "<<
	    (*patPhotonsIso)[0].genPhoton()->mother()->status();
	}
	else
	  std::cout<<"0    0    0    0";
	std::cout<<std::endl;
	std::cout<<"patPhotonsIso[1]::"<<
	  (*patPhotonsIso)[1].pt()            <<"    "<<
	  (*patPhotonsIso)[1].eta()           <<"    "<<
	  (*patPhotonsIso)[1].hasPixelSeed()  <<"    "<<
	  (*patPhotonsIso)[1].hadronicOverEm()<<"    "<<
	  (*patPhotonsIso)[1].hadTowOverEm()  <<"    "<<
	  (*patPhotonsIso)[1].sigmaIetaIeta() <<"    "<<
	  (*patPhotonsIso)[1].userFloat("pfChargedPURel")<<"    "<<
	  (*patPhotonsIso)[1].userFloat("pfNeutralPURel")<<"    "<<
	  (*patPhotonsIso)[1].userFloat("pfGammaPURel")  <<"    ";
	if (!data_) {
	  std::cout<<
	    (*patPhotonsIso)[1].genPhoton()->pdgId() <<"    "<<
	    (*patPhotonsIso)[1].genPhoton()->status()<<"    "<<
	    (*patPhotonsIso)[1].genPhoton()->mother()->pdgId()<<"    "<<
	    (*patPhotonsIso)[1].genPhoton()->mother()->status();
	}
	else
	  std::cout<<"0    0    0    0";
	std::cout<<std::endl;
      }
    }
  }

  int n_jets_pt30 = 0, n_jets_pt50eta25 = 0;
  double ht  = 0.0;
  reco::MET::LorentzVector mht(0,0,0,0);
  for(unsigned int i=0; i<Jets.size(); ++i) {
    const pat::Jet *r = Jets[i];
    if(Jets[i]->pt() > 50.0 && fabs(Jets[i]->eta()) < 2.50) {
      n_jets_pt50eta25++;
      ht  += Jets[i]->pt();
    }
    if(Jets[i]->pt() > 30.0 && fabs(Jets[i]->eta()) < 5.0) { 
      n_jets_pt30++;
      mht  -= Jets[i]->p4();
    }
  }
  reco::MET MHT = reco::MET(mht, reco::MET::Point());
  double mht_val =  MHT.pt();

  const pat::Jet  *r, *r1, *r2, *r3;
  double dphi_mht_j1 = 0, dphi_mht_j2 = 0, dphi_mht_j3 = 0;
  if(n_jets_pt30 >= 3) {
    r1 = Jets[0]; r2 = Jets[1]; r3 = Jets[2];
    dphi_mht_j1 = fabs(reco::deltaPhi(r1->phi(),MHT.phi()));
    dphi_mht_j2 = fabs(reco::deltaPhi(r2->phi(),MHT.phi()));
    dphi_mht_j3 = fabs(reco::deltaPhi(r3->phi(),MHT.phi()));
  }

  h_PhotonsIso_NJets_Pt30     ->Fill(n_jets_pt30, pu_event_wt);
  h_PhotonsIso_NJets_Pt50Eta25->Fill(n_jets_pt50eta25, pu_event_wt);
  if(n_jets_pt30 >= 3) {
    h_PhotonsIso_DPhiMHTJet1->Fill(dphi_mht_j1, pu_event_wt);
    h_PhotonsIso_DPhiMHTJet2->Fill(dphi_mht_j2, pu_event_wt);
    h_PhotonsIso_DPhiMHTJet3->Fill(dphi_mht_j3, pu_event_wt);
  }
  h_PhotonsIso_HT ->Fill(ht, pu_event_wt);
  h_PhotonsIso_MHT->Fill(mht_val, pu_event_wt);

  // do RA2 analysis here
  //if(n_jets_pt50eta25 >= 3) {
  if(n_jets_pt50eta25 >= ra2NJets_) {
    
    //if( ht> 500.0) {
    if( ht> ra2HT_) {
      
      //if(mht_val > 200.0) {
      if(mht_val > ra2MHT_) {
	
	h_preRA2_DPhiMHTJet1->Fill(dphi_mht_j1, pu_event_wt);
	h_preRA2_DPhiMHTJet2->Fill(dphi_mht_j2, pu_event_wt);
	h_preRA2_DPhiMHTJet3->Fill(dphi_mht_j3, pu_event_wt);

	if( !ra2ApplyDphiCuts_ || (ra2ApplyDphiCuts_ && (dphi_mht_j1 > 0.5 && dphi_mht_j2 > 0.5 && dphi_mht_j3 > 0.3)) ) {
	  
	  h_RA2_Vertices ->Fill(nVertices);
	  h_RA2_VerticesReWeighted ->Fill(nVertices, pu_event_wt);

	  double pt  = RA2Photons[0]->pt();
	  double eta = RA2Photons[0]->eta();
	  double phi = RA2Photons[0]->phi();	  

	  double totalPhotonIso    = RA2Photons[0]->trkSumPtHollowConeDR03()
	                            +RA2Photons[0]->ecalRecHitSumEtConeDR03()
	                            +RA2Photons[0]->hcalTowerSumEtConeDR03();
	  double totalPhotonIso2012    = RA2Photons[0]->trkSumPtHollowConeDR03()
	                                +RA2Photons[0]->ecalRecHitSumEtConeDR03()
	                                +RA2Photons[0]->userFloat("hcalIsoConeDR03_2012");
	  double totalPhotonRhoIso     = totalPhotonIso     - areaR03*pf_event_rho;
	  double totalPhotonRhoIso2012 = totalPhotonIso2012 - areaR03*pf_event_rho;
	  
	  //if( nVertices<60 ) {
	  //  h_totalIsolation[nVertices]->Fill(totalPhotonIso);
	  //  h_totalIsolationRho[nVertices]->Fill(totalPhotonRhoIso);
	  //}

	  if( debug_) {
	    std::cout << "RA2 " <<std::setw(9)<<run<<":"<<lumi<<":"<<event
		      << std::setw(5)<<" == MHT::"<<mht_val<<std::setw(5)<<" == HT::"<<ht
		      << std::setw(5)<<" == pt::"<<std::setw(9)<<pt<<std::setw(5)<<" == eta::"<<std::setw(8)<<eta
		      << std::setw(3)<<" == R9::"<<std::setw(12)<<RA2Photons[0]->r9()
		      << std::setprecision(4)
		      << std::setw(8) <<" == hOverE::"      <<std::setw(12)<<RA2Photons[0]->hadronicOverEm() 
		      << std::setw(8) <<" == hOverE20::"    <<std::setw(12)<<RA2Photons[0]->hadTowOverEm() 
		      << std::setw(11)<<" == sigEtaEta::"   <<std::setw(10)<<RA2Photons[0]->sigmaIetaIeta()
		      << std::setw(8) <<" == trkIso::"      <<std::setw(10)<<RA2Photons[0]->trkSumPtHollowConeDR03()
		      << std::setw(9) <<" == ecalIso::"     <<std::setw(10)<<RA2Photons[0]->ecalRecHitSumEtConeDR03()
		      << std::setw(9) <<" == hcalIso::"     <<std::setw(10)<<RA2Photons[0]->hcalTowerSumEtConeDR03()
		      << std::setw(9) <<" == hcalIso2012::" <<std::setw(10)<<RA2Photons[0]->userFloat("hcalIsoConeDR03_2012")
		      << std::endl

		      << std::setw(9) <<" == pf_rho::"     <<std::setw(10)<<pf_event_rho
		      << std::setw(9) <<" == pfRho25::"    <<std::setw(10)<<RA2Photons[0]->userFloat("rho25")
		      << std::setw(9) <<" == pfChargedEA::"<<std::setw(10)<<RA2Photons[0]->userFloat("pfChargedEA")
		      << std::setw(9) <<" == pfNeutralEA::"<<std::setw(10)<<RA2Photons[0]->userFloat("pfNeutralEA")
		      << std::setw(9) <<" == pfGammaEA::"  <<std::setw(10)<<RA2Photons[0]->userFloat("pfGammaEA")
		      << std::setw(9) <<" == pfChargedEA*pfRho25::"<<std::setw(10)<<RA2Photons[0]->userFloat("pfChargedEA")*RA2Photons[0]->userFloat("rho25")
		      << std::setw(9) <<" == pfNeutralEA*pfRho25::"<<std::setw(10)<<RA2Photons[0]->userFloat("pfNeutralEA")*RA2Photons[0]->userFloat("rho25")
		      << std::setw(9) <<" == pfGammaEA*pfRho25::"  <<std::setw(10)<<RA2Photons[0]->userFloat("pfGammaEA")*RA2Photons[0]->userFloat("rho25")
		      << std::endl

		      << std::setw(9) <<" == pfCharged::"   <<std::setw(10)<<RA2Photons[0]->userIsolation(pat::IsolationKeys(pat::User1Iso))
		      << std::setw(9) <<" == pfNeutral::"   <<std::setw(10)<<RA2Photons[0]->userIsolation(pat::IsolationKeys(pat::User3Iso))
		      << std::setw(9) <<" == pfGamma::"     <<std::setw(10)<<RA2Photons[0]->userIsolation(pat::IsolationKeys(pat::User4Iso))

		      << std::setw(9) <<" == pfChargedRel::"   <<std::setw(10)<<RA2Photons[0]->userFloat("pfChargedRel")
		      << std::setw(9) <<" == pfNeutralRel::"   <<std::setw(10)<<RA2Photons[0]->userFloat("pfNeutralRel")
		      << std::setw(9) <<" == pfGammaRel::"     <<std::setw(10)<<RA2Photons[0]->userFloat("pfGammaRel")
		      << std::setw(9) <<" == pfChargedPURel::" <<std::setw(10)<<RA2Photons[0]->userFloat("pfChargedPURel")
		      << std::setw(9) <<" == pfNeutralPURel::" <<std::setw(10)<<RA2Photons[0]->userFloat("pfNeutralPURel")
		      << std::setw(9) <<" == pfGammaPURel::"   <<std::setw(10)<<RA2Photons[0]->userFloat("pfGammaPURel")
		      << std::setw(9) <<" == pfChargedPUSub::" <<std::setw(10)<<RA2Photons[0]->userFloat("pfChargedPUSub")
		      << std::setw(9) <<" == pfNeutralPUSub::" <<std::setw(10)<<RA2Photons[0]->userFloat("pfNeutralPUSub")
		      << std::setw(9) <<" == pfGammaPUSub::"   <<std::setw(10)<<RA2Photons[0]->userFloat("pfGammaPUSub")
		      << std::endl

		      << std::setw(9) <<" == totalIso::"   <<std::setw(10)<<totalPhotonIso
		      << std::setw(9) <<" == totalIsoRho::"<<std::setw(10)<<totalPhotonRhoIso
		      << std::endl;
	  }

	  h_RA2_NPhotons->Fill(RA2Photons.size(), pu_event_wt); 
	  for( unsigned int iphot=0; iphot<RA2Photons.size(); iphot++) {
	    if( iphot<2 ) {
	      h_RA2_EtaVsPt[iphot]    ->Fill(RA2Photons[iphot]->eta(), RA2Photons[iphot]->phi(), pu_event_wt);
	      h_RA2_Pt [iphot]        ->Fill(RA2Photons[iphot]->pt() , pu_event_wt);
	      h_RA2_Eta[iphot]        ->Fill(RA2Photons[iphot]->eta(), pu_event_wt);
	      h_RA2_Phi[iphot]        ->Fill(RA2Photons[iphot]->phi(), pu_event_wt);

	      h_RA2_TrkIso[iphot]     ->Fill(RA2Photons[iphot]->trkSumPtHollowConeDR03(),  pu_event_wt);
	      h_RA2_EcalIso[iphot]    ->Fill(RA2Photons[iphot]->ecalRecHitSumEtConeDR03(), pu_event_wt);
	      h_RA2_HcalIso[iphot]    ->Fill(RA2Photons[iphot]->hcalTowerSumEtConeDR03(),  pu_event_wt);
	      h_RA2_HcalIso2012[iphot]->Fill(RA2Photons[iphot]->userFloat("hcalIsoConeDR03_2012"), pu_event_wt);

	      h_RA2_PFCharged[iphot]->Fill(RA2Photons[iphot]->userIsolation(pat::IsolationKeys(pat::UserBaseIso+0)), pu_event_wt);
	      h_RA2_PFNeutral[iphot]->Fill(RA2Photons[iphot]->userIsolation(pat::IsolationKeys(pat::UserBaseIso+2)), pu_event_wt);
	      h_RA2_PFGamma[iphot]  ->Fill(RA2Photons[iphot]->userIsolation(pat::IsolationKeys(pat::UserBaseIso+3)), pu_event_wt);

	      double pfChargedIsoRel = RA2Photons[iphot]->userFloat("pfChargedRel");
	      double pfNeutralIsoRel = RA2Photons[iphot]->userFloat("pfNeutralRel");
	      double pfGammaIsoRel   = RA2Photons[iphot]->userFloat("pfGammaRel");
	      
	      double pfChargedIsoPURel = RA2Photons[iphot]->userFloat("pfChargedPURel");
	      double pfNeutralIsoPURel = RA2Photons[iphot]->userFloat("pfNeutralPURel");
	      double pfGammaIsoPURel   = RA2Photons[iphot]->userFloat("pfGammaPURel");
	      
	      h_RA2_PFChargedRel[iphot] ->Fill(pfChargedIsoRel, pu_event_wt);
	      h_RA2_PFNeutralRel[iphot] ->Fill(pfNeutralIsoRel, pu_event_wt);
	      h_RA2_PFGammaRel[iphot]   ->Fill(pfGammaIsoRel,   pu_event_wt);
	      
	      h_RA2_PFChargedRelPU[iphot] ->Fill(pfChargedIsoPURel, pu_event_wt);
	      h_RA2_PFNeutralRelPU[iphot] ->Fill(pfNeutralIsoPURel, pu_event_wt);
	      h_RA2_PFGammaRelPU[iphot]   ->Fill(pfGammaIsoPURel,   pu_event_wt);

	      h_RA2_hOverE[iphot]     ->Fill(RA2Photons[iphot]->hadronicOverEm(), pu_event_wt);
	      h_RA2_hOverE2012[iphot] ->Fill(RA2Photons[iphot]->hadTowOverEm(),   pu_event_wt);
	      h_RA2_R9[iphot]         ->Fill(RA2Photons[iphot]->r9(),             pu_event_wt);

	      double eta       = RA2Photons[iphot]->eta();
	      double sigEtaEta = RA2Photons[iphot]->sigmaIetaIeta();
	      if(std::fabs(eta)<1.4442)     h_RA2_SigEtaEta_EB[iphot]->Fill(sigEtaEta, pu_event_wt);
	      else if(std::fabs(eta)>1.566) h_RA2_SigEtaEta_EE[iphot]->Fill(sigEtaEta, pu_event_wt);
	    }
	  }

	  for(int ipt=0; ipt<NPhotPtBins-1; ipt++) {
	    for(int ieta=0; ieta<NPhotEtaBins-1; ieta++) {
	      if( photon_pt>PhotPtVal[ipt] && photon_pt<PhotPtVal[ipt+1] ) {
		if( std::fabs(photon_eta)>PhotEtaVal[ieta] && std::fabs(photon_eta)<PhotEtaVal[ieta+1] ) {
		  NPhotonPtEta_RA2[ipt][ieta] += pu_event_wt;
		}
	      }
	    }
	  }

	  // get the event yield in inclusive bins
	  if( ht>350.0 && mht_val>200.0 ) {
	    h_RA2_NJets_HT350_MHT200 ->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT350_MHT200    ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT350_MHT200   ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT350_MHT200  ->Fill(ht+mht_val,       pu_event_wt);
	    h_RA2_PhotPtEta_HT350_MHT200->Fill( std::fabs(photon_eta), photon_pt, pu_event_wt);
	  }
	  if( ht>800.0 && mht_val>200.0 ) {
	    h_RA2_NJets_HT800_MHT200 ->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT800_MHT200    ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT800_MHT200   ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT800_MHT200  ->Fill(ht+mht_val,       pu_event_wt);
	    h_RA2_PhotPtEta_HT800_MHT200->Fill( std::fabs(photon_eta), photon_pt, pu_event_wt);
	  }
	  if( ht>800.0 && mht_val>500.0 ) {
	    h_RA2_NJets_HT800_MHT500 ->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT800_MHT500    ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT800_MHT500   ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT800_MHT500  ->Fill(ht+mht_val,       pu_event_wt);
	    h_RA2_PhotPtEta_HT800_MHT500->Fill( std::fabs(photon_eta), photon_pt, pu_event_wt);
	  }
	  if( ht>500.0 && mht_val>350.0 ) {
	    h_RA2_NJets_HT500_MHT350 ->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT500_MHT350    ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT500_MHT350   ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT500_MHT350  ->Fill(ht+mht_val,       pu_event_wt);
	    h_RA2_PhotPtEta_HT500_MHT350->Fill( std::fabs(photon_eta), photon_pt, pu_event_wt);
	  }
	  if( ht>500.0 && mht_val>200.0 ) {
	    h_RA2_NJets_HT500_MHT200 ->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT500_MHT200    ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT500_MHT200   ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT500_MHT200  ->Fill(ht+mht_val,       pu_event_wt);
	    h_RA2_PhotPtEta_HT500_MHT200->Fill( std::fabs(photon_eta), photon_pt, pu_event_wt);
	  }
	  if( ht>1000.0 && mht_val>600.0 ) {
	    h_RA2_NJets_HT1000_MHT600->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT1000_MHT600   ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT1000_MHT600  ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT1000_MHT600 ->Fill(ht+mht_val,       pu_event_wt);
	    h_RA2_PhotPtEta_HT1000_MHT600->Fill( std::fabs(photon_eta), photon_pt, pu_event_wt);
	  }
	  if( ht>1200.0 && mht_val>400.0 ) {
	    h_RA2_NJets_HT1200_MHT400->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT1200_MHT400   ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT1200_MHT400  ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT1200_MHT400 ->Fill(ht+mht_val,       pu_event_wt);
	    h_RA2_PhotPtEta_HT1200_MHT400->Fill( std::fabs(photon_eta), photon_pt, pu_event_wt);
	  }
	  
	  
	  // get event yields in exclusive bins
	  h_RA2_HTMHT_Excl->Fill( ht, mht_val, pu_event_wt);
	  int ibin = h_RA2_HTMHT_Excl->FindBin(ht, mht_val);
	  int ihtBin, imhtBin, iZ;
	  h_RA2_HTMHT_Excl->GetBinXYZ(ibin, ihtBin, imhtBin, iZ);
	  if( ihtBin >h_RA2_HTMHT_Excl->GetNbinsX() ) ihtBin =h_RA2_HTMHT_Excl->GetNbinsX();
	  if( imhtBin>h_RA2_HTMHT_Excl->GetNbinsY() ) imhtBin=h_RA2_HTMHT_Excl->GetNbinsY();
	  //std::cout << "ibin "<< ibin << " " << ihtBin << " " << imhtBin << std::endl;
	  if( ihtBin>0 && imhtBin>0 )
	    h_RA2_PhotPtEta_Excl[ihtBin-1][imhtBin-1]->Fill( std::fabs(photon_eta), photon_pt, pu_event_wt);


	  // get event yield in inclusive bins
	  int ibinInc = h_RA2_HTMHT_Incl->FindBin(ht, mht_val);
	  int ihtBinInc, imhtBinInc, iZInc;
	  h_RA2_HTMHT_Incl->GetBinXYZ(ibinInc, ihtBinInc, imhtBinInc, iZInc);
	  if( ihtBinInc >h_RA2_HTMHT_Incl->GetNbinsX() ) ihtBinInc =h_RA2_HTMHT_Incl->GetNbinsX();
	  if( imhtBinInc>h_RA2_HTMHT_Incl->GetNbinsY() ) imhtBinInc=h_RA2_HTMHT_Incl->GetNbinsY();
	  for(int iht=1; iht<=ihtBinInc; iht++) {
	    for(int imht=1; imht<=imhtBinInc; imht++) {
	      double binHT  = h_RA2_HTMHT_Incl->GetXaxis()->GetBinCenter(iht);
	      double binMHT = h_RA2_HTMHT_Incl->GetYaxis()->GetBinCenter(imht);
	      h_RA2_HTMHT_Incl->Fill(binHT, binMHT, pu_event_wt);
	    }
	  }


	  // just some monitoring plots
	  h_RA2_NJets_Pt30     ->Fill(n_jets_pt30, pu_event_wt);
	  h_RA2_NJets_Pt50Eta25->Fill(n_jets_pt50eta25, pu_event_wt);
	  if(n_jets_pt30 >= 3) {
	    h_RA2_DPhiMHTJet1->Fill(dphi_mht_j1, pu_event_wt);
	    h_RA2_DPhiMHTJet2->Fill(dphi_mht_j2, pu_event_wt);
	    h_RA2_DPhiMHTJet3->Fill(dphi_mht_j3, pu_event_wt);
	  }
	  h_RA2_HT ->Fill(ht , pu_event_wt);
	  h_RA2_MHT->Fill(mht_val, pu_event_wt);
	  h_RA2_MEff->Fill(ht+mht_val, pu_event_wt);
	  
	} 
	
      }

    }
  }
}

void RA2ZInvPhotonAnalyzer::doGenAnalysis(edm::Handle<reco::GenParticleCollection>& genParticles, 
					  edm::Handle<edm::View<pat::Jet> >& jetsHT,
					  edm::Handle<edm::View<pat::Jet> >& jetsMHT,
					  std::vector<const pat::Photon*>    IsoPhotons,
					  double pu_event_wt) {
  pu_event_wt = 1.0;

  bool electronMistag= false;
  unsigned int ipart = 0;

  for(reco::GenParticleCollection::const_iterator part = genParticles->begin(); part<genParticles->end(); ++part){

    const reco::GenParticle *mother = (reco::GenParticle *) part->mother();

    if( part->pt()>10.0 ) {

      // definition of direct photon 
      bool directPhoton    = ( std::fabs(part->pdgId()) == 22 && part->status() == 3 &&
			       (std::fabs(mother->pdgId())<10 || mother->pdgId()==21 )  ) ;

      // definition of secondary photon or decay photons
      bool secondaryPhoton = ( std::fabs(part->pdgId()) == 22 && part->status() == 1 && std::fabs(mother->pdgId())>100 );

      // definition of fragmentation photon
      bool fragmenPhoton   = ( std::fabs(part->pdgId()) == 22 && part->status() == 1 && (std::fabs(mother->pdgId())<10 || mother->pdgId()==21 )  );
      
      // definition of electron mistags
      // to be done
      
      // consider only direct photons here
      if( directPhoton ) {
	
	double phot_pt  = part->pt();
	double phot_eta = part->eta();
	double phot_phi = part->phi();
	
	// check if genPhoton is within acceptance 
	//why these |eta| cuts and why this pT cut?
	bool kineAccPhot = ( phot_pt>10.0 && 
			     ((std::fabs(phot_eta)<1.4442) || (std::fabs(phot_eta)>1.566 && 
							      std::fabs(phot_eta)<2.5)) );
	///these are our nominal cuts
	//bool kineAccPhot = ( phot_pt>100.0 && 
	//		     ((std::fabs(phot_eta)<1.4442) || (std::fabs(phot_eta)>1.566 && 
	//						     std::fabs(phot_eta)<2.5)) );
	
	// find a matching recoPhoton to check efficiency
	bool   foundPhotonMatch = false;
	double dRMinPhot=99.0;
	int    dRMinPhotIdx = -1;
	for( unsigned int iphot=0; iphot<IsoPhotons.size(); iphot++) {	      
	  double rphot_pt  = IsoPhotons[iphot]->pt();
	  double rphot_eta = IsoPhotons[iphot]->eta();
	  double rphot_phi = IsoPhotons[iphot]->phi();
	  double dR = reco::deltaR(rphot_eta, rphot_phi, phot_eta, phot_phi);
	  if( dR<dRMinPhot) {
	    dRMinPhot    = dR;
	    dRMinPhotIdx = iphot;
	  }
	}
	if( dRMinPhotIdx>=0 && dRMinPhot<0.1 ) {
	  foundPhotonMatch = true;
	} else {
	  //std::cout << "not Matching to leading Phot " << " dRMin " << dRMin << " NRecoPhotons " << IsoPhotons.size() << std::endl;	    
	}

	// find & remove closest reconstructed jet
	double dRMinJet    = 99.0;
	int    dRMinJetIdx = -1;
	for(unsigned int ijet=0; ijet<jetsMHT->size(); ijet++) {

	  const pat::Jet *r = &((*jetsMHT)[ijet]);
	  double jet_pt  = r->pt();
	  double jet_eta = r->eta();
	  double jet_phi = r->phi();	  

	  if( jet_pt>30.0 ) {
	    double dR = reco::deltaR(phot_eta, phot_phi, jet_eta, jet_phi);
	    if(dR<dRMinJet) { 
	      dRMinJet    = dR;
	      dRMinJetIdx = ijet;
	    }
	  }
	}
	
	// clean jet collection	
	std::vector<const pat::Jet*> jetsClean; 
	if(dRMinJetIdx>=0 && dRMinJet<0.1 ) {
	  for(unsigned int i=0; i<jetsMHT->size(); ++i) {
	    const pat::Jet *r = &((*jetsMHT)[i]);
	    int jj = (int) i;
	    if( jj != dRMinJetIdx ) jetsClean.push_back(r);
	  }
	}
	
	// calculate HT/MHT here
	int n_jets_pt50eta25 = 0;
	double ht  = 0.0;
	double ht_jet30 = 0.0;
	reco::MET::LorentzVector mhtV(0,0,0,0);

	bool goodJetId = true;
	for(unsigned int i=0; i<jetsClean.size(); ++i) {
	  const pat::Jet *r = jetsClean[i];
	  if(jetsClean[i]->pt() > 50.0 && std::fabs(jetsClean[i]->eta())<2.50) {
	    n_jets_pt50eta25++;
	    ht  += jetsClean[i]->pt();
	  }
	  if(jetsClean[i]->pt() > 30.0 && std::fabs(jetsClean[i]->eta())<5.0) { 
	    if( jetsClean[i]->neutralHadronEnergyFraction()>0.90 || jetsClean[i]->photonEnergyFraction()/jetsClean[i]->jecFactor(0)>0.95 ) goodJetId=false; 
	    mhtV     -= jetsClean[i]->p4();
	    ht_jet30 += jetsClean[i]->pt();
	  }
	}
	double mht    = mhtV.pt();
	double mhtPhi = mhtV.phi();
	
	h_GenBoson_NJets_Pt30     ->Fill( jetsMHT->size() );
	h_GenBoson_NJets_Pt50Eta25->Fill( n_jets_pt50eta25);
	h_GenBoson_PtVsHT         ->Fill( ht_jet30, part->pt() );
	h_GenBoson_MHTVsHT        ->Fill( ht_jet30, mht );
	h_GenBoson_Pt             ->Fill( part->pt() );
	h_GenBoson_Eta            ->Fill( part->eta() );
	h_GenBoson_Phi            ->Fill( part->phi() );
	h_GenBoson_HT             ->Fill( ht );
	h_GenBoson_MHT            ->Fill( mht );



	// do RA2 analysis here
	double dphi_mht_j1 = 0, dphi_mht_j2 = 0, dphi_mht_j3 = 0;
	if( jetsClean.size()> 2 ){
	  dphi_mht_j1 = std::fabs(reco::deltaPhi(jetsClean[0]->phi(),mhtPhi));
	  dphi_mht_j2 = std::fabs(reco::deltaPhi(jetsClean[1]->phi(),mhtPhi));
	  dphi_mht_j3 = std::fabs(reco::deltaPhi(jetsClean[2]->phi(),mhtPhi));

	}
	if(goodJetId && n_jets_pt50eta25 >= ra2NJets_) {
   
	  if( ht> ra2HT_) {
      
	    if(mht > ra2MHT_) {

	      h_GenBoson_preRA2_dPhiMHTJet1->Fill(dphi_mht_j1);
	      h_GenBoson_preRA2_dPhiMHTJet2->Fill(dphi_mht_j2);
	      h_GenBoson_preRA2_dPhiMHTJet3->Fill(dphi_mht_j3);
	
	      if( !ra2ApplyDphiCuts_ || (ra2ApplyDphiCuts_ && (dphi_mht_j1 > 0.5 && dphi_mht_j2 > 0.5 && dphi_mht_j3 > 0.3)) ) {
		h_GenBoson_RA2_Pt             ->Fill( part->pt() );
		h_GenBoson_RA2_Eta            ->Fill( part->eta() );
		h_GenBoson_RA2_Phi            ->Fill( part->phi() );
		h_GenBoson_RA2_HT             ->Fill( ht );
		h_GenBoson_RA2_MHT            ->Fill( mht );
		
		// get event yields in exclusive bins
		h_GenBosonTotal_RA2_HTMHT_Excl->Fill( ht, mht);

		int ibinInc = h_GenBosonTotal_RA2_HTMHT_Incl->FindBin(ht, mht);
		int ihtBinInc, imhtBinInc, iZInc;
		h_GenBosonTotal_RA2_HTMHT_Incl->GetBinXYZ(ibinInc, ihtBinInc, imhtBinInc, iZInc);
		if( ihtBinInc >h_GenBosonTotal_RA2_HTMHT_Incl->GetNbinsX() ) 
		  ihtBinInc =h_GenBosonTotal_RA2_HTMHT_Incl->GetNbinsX();
		if( imhtBinInc>h_GenBosonTotal_RA2_HTMHT_Incl->GetNbinsY() ) 
		  imhtBinInc=h_GenBosonTotal_RA2_HTMHT_Incl->GetNbinsY();
		for(int iht=1; iht<=ihtBinInc; iht++) {
		  for(int imht=1; imht<=imhtBinInc; imht++) {
		    double binHT  = h_GenBosonTotal_RA2_HTMHT_Incl->GetXaxis()->GetBinCenter(iht);
		    double binMHT = h_GenBosonTotal_RA2_HTMHT_Incl->GetYaxis()->GetBinCenter(imht);
		    h_GenBosonTotal_RA2_HTMHT_Incl->Fill(binHT, binMHT, pu_event_wt);
		  }
		}
		if( kineAccPhot ) {
		  h_GenBosonKineAcc_RA2_HTMHT_Excl->Fill( ht, mht);
		  for(int iht=1; iht<=ihtBinInc; iht++) {
		    for(int imht=1; imht<=imhtBinInc; imht++) {
		      double binHT  = h_GenBosonKineAcc_RA2_HTMHT_Incl->GetXaxis()->GetBinCenter(iht);
		      double binMHT = h_GenBosonKineAcc_RA2_HTMHT_Incl->GetYaxis()->GetBinCenter(imht);
		      h_GenBosonKineAcc_RA2_HTMHT_Incl->Fill(binHT, binMHT, pu_event_wt);
		    }
		  }
		}
		if( foundPhotonMatch ) {
		  h_GenBosonMatched_RA2_HTMHT_Excl->Fill( ht, mht);
		  for(int iht=1; iht<=ihtBinInc; iht++) {
		    for(int imht=1; imht<=imhtBinInc; imht++) {
		      double binHT  = h_GenBosonMatched_RA2_HTMHT_Incl->GetXaxis()->GetBinCenter(iht);
		      double binMHT = h_GenBosonMatched_RA2_HTMHT_Incl->GetYaxis()->GetBinCenter(imht);
		      h_GenBosonMatched_RA2_HTMHT_Incl->Fill(binHT, binMHT, pu_event_wt);
		    }
		  }
		}
	      }
	    }
	  }
	}
     


      } // if direct photon

      /*
      if( debug_ ) {
	if( directPhoton )    std::cout << "It's a Direct Photon "        << std::endl;
	if( secondaryPhoton ) std::cout << "It's a Secondary Photon "    << std::endl;
	if( fragmenPhoton )   std::cout << "It's a Fragmentation Photon " << std::endl;
	if( directPhoton || secondaryPhoton || fragmenPhoton ) {
	  std::cout << ipart << " " << part->pdgId() << " status " << part->status() 
		    << " pt " << part->pt() << " eta " << part->eta() << " phi " << part->phi()
		    << " mother pt " << mother->pt() << " mother Status " << mother->status()
		    << " mother pdgId " << mother->pdgId()
		    << std::endl;
	}
      }
      */
    } // if a parton with pt>10
  
    ipart++;
  } // looop over gen particles

}







void RA2ZInvPhotonAnalyzer::beginJob() {

  double RA2MHTValTemp[NRA2MHTBins] = {200, 350, 500, 600,  800};
  double RA2HTValTemp [NRA2HTBins ] = {350, 500, 800, 1000, 1200, 1400, 1600};
  for(int imht=0; imht<NRA2MHTBins; imht++) RA2MHTVal[imht] = RA2MHTValTemp[imht];
  for(int iht=0;  iht<NRA2HTBins;   iht++)  RA2HTVal [iht]  = RA2HTValTemp[iht];
  
  for(int iht=0;  iht<NRA2HTBins;   iht++){
    for(int imht=0; imht<NRA2MHTBins; imht++) {
      std::cout << "(" << RA2HTVal[iht] << "," << RA2MHTVal[imht] << ")" <<"  ";
    }
    std::cout << std::endl;
  }

  
  double MHTValTemp[NMHTBins] = {200, 250, 300, 350, 400,  500,  600,  800};
  double HTValTemp [NHTBins]  = {350, 400, 500, 600, 800, 1000, 1200, 1400};
  for(int imht=0; imht<NMHTBins; imht++) MHTVal[imht] = MHTValTemp[imht];
  for(int iht=0;  iht<NHTBins;   iht++)  HTVal [iht]  = HTValTemp [iht];
  for(int iht=0;  iht<NHTBins;   iht++){
    for(int imht=0; imht<NMHTBins; imht++) {
      std::cout << "(" << HTVal[iht] << "," << MHTVal[imht] << ")" <<"  ";
    }
    std::cout << std::endl;
  }
 

  double tempEtaVal[NPhotEtaBins] = {0.0,  0.9, 1.442, 1.566, 2.1, 2.5, 5.0};
  double tempPtVal [NPhotPtBins]  = {80.0, 100.0, 120.0, 200.0, 300.0, 400.0, 500.0, 1000.0};
  for(int ipt =0; ipt <NPhotPtBins;  ipt++)  PhotPtVal [ipt]  = tempPtVal[ipt];
  for(int ieta=0; ieta<NPhotEtaBins; ieta++) PhotEtaVal[ieta] = tempEtaVal[ieta];

  for(int ipt=0; ipt<NPhotPtBins; ipt++) {
    for(int ieta=0; ieta<NPhotEtaBins; ieta++) {
      NPhotonPtEta_Iso[ipt][ieta] = 0.0;
      NPhotonPtEta_RA2[ipt][ieta] = 0.0;
    }
  }

  //book histograms
  BookHistograms();
}

void RA2ZInvPhotonAnalyzer::endJob() {

  // treat the overflows properly 
  for(int imht=1; imht<=h_RA2_HTMHT_Excl->GetNbinsY(); imht++) {
    int    overbin    = h_RA2_HTMHT_Excl->GetNbinsX()+1;
    double overflow   = h_RA2_HTMHT_Excl->GetBinContent(overbin, imht);
    int    lastBinMHT = h_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(imht);
    int    lastBinHT  = h_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(overbin-1);
    h_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);

    int    thisBin    = h_RA2_HTMHT_Excl->GetBin(overbin, imht);
    h_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
  }

  for(int iht=1; iht<=h_RA2_HTMHT_Excl->GetNbinsX(); iht++){
    int    overbin    = h_RA2_HTMHT_Excl->GetNbinsY()+1;
    double overflow   = h_RA2_HTMHT_Excl->GetBinContent(iht, overbin);
    double lastBinMHT = h_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(overbin-1);
    double lastBinHT  = h_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(iht);
    h_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);

    int    thisBin  = h_RA2_HTMHT_Excl->GetBin(iht, overbin);
    h_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
  }

  for(int ipt=0; ipt<NPhotPtBins-1; ipt++) {
    for(int ieta=0; ieta<NPhotEtaBins-1; ieta++) {
      h_PhotonIso_PhotPtEta->SetBinContent( h_PhotonIso_PhotPtEta->GetBin(ipt+1, ieta+1), NPhotonPtEta_Iso[ipt][ieta] );	
      h_RA2_PhotPtEta      ->SetBinContent( h_RA2_PhotPtEta->GetBin(ipt+1, ieta+1),       NPhotonPtEta_RA2[ipt][ieta] );	
    }
  }

  for(int imht=1; imht<=h_RA2_HTMHT_Excl->GetNbinsY(); imht++) {
    int    overbin    = h_RA2_HTMHT_Excl->GetNbinsX()+1;
    double overflow   = h_RA2_HTMHT_Excl->GetBinContent(overbin, imht);
    int    lastBinMHT = h_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(imht);
    int    lastBinHT  = h_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(overbin-1);
    h_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);

    int    thisBin    = h_RA2_HTMHT_Excl->GetBin(overbin, imht);
    h_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
  }

  if( doGenAnalysis_) { 
    for(int iht=1; iht<=h_GenBosonTotal_RA2_HTMHT_Excl->GetNbinsX(); iht++){
      int    overbin    = h_GenBosonTotal_RA2_HTMHT_Excl->GetNbinsY()+1;
      double overflow   = h_GenBosonTotal_RA2_HTMHT_Excl->GetBinContent(iht, overbin);
      double lastBinMHT = h_GenBosonTotal_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(overbin-1);
      double lastBinHT  = h_GenBosonTotal_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(iht);
      h_GenBosonTotal_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);
      
      int    thisBin  = h_GenBosonTotal_RA2_HTMHT_Excl->GetBin(iht, overbin);
      h_GenBosonTotal_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
    }

    for(int iht=1; iht<=h_GenBosonKineAcc_RA2_HTMHT_Excl->GetNbinsX(); iht++){
      int    overbin    = h_GenBosonKineAcc_RA2_HTMHT_Excl->GetNbinsY()+1;
      double overflow   = h_GenBosonKineAcc_RA2_HTMHT_Excl->GetBinContent(iht, overbin);
      double lastBinMHT = h_GenBosonKineAcc_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(overbin-1);
      double lastBinHT  = h_GenBosonKineAcc_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(iht);
      h_GenBosonKineAcc_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);
      
      int    thisBin  = h_GenBosonKineAcc_RA2_HTMHT_Excl->GetBin(iht, overbin);
      h_GenBosonKineAcc_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
    }

    for(int iht=1; iht<=h_GenBosonMatched_RA2_HTMHT_Excl->GetNbinsX(); iht++){
      int    overbin    = h_GenBosonMatched_RA2_HTMHT_Excl->GetNbinsY()+1;
      double overflow   = h_GenBosonMatched_RA2_HTMHT_Excl->GetBinContent(iht, overbin);
      double lastBinMHT = h_GenBosonMatched_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(overbin-1);
      double lastBinHT  = h_GenBosonMatched_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(iht);
      h_GenBosonMatched_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);
      
      int    thisBin  = h_GenBosonMatched_RA2_HTMHT_Excl->GetBin(iht, overbin);
      h_GenBosonMatched_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
    }

  }

}


void RA2ZInvPhotonAnalyzer::BookHistograms() {

  char hname[200], htit[200];
  h_puWeight = fs->make<TH1F>("h_puWeight", "h_puWeight", 20, 0.0, 2.0);
  h_Vertices = fs->make<TH1F>("h_Vertices", "h_Vertices", 60, 0.0, 60.0);
  h_VerticesReWeighted = fs->make<TH1F>("h_VerticesReWeighted", "h_VerticesReWeighted", 60, 0.0, 60.0);

  for(int i=0; i<60; i++){
    sprintf(hname, "h_pfRho_nVtx_%i", i+1);
    sprintf(htit, "pfRho nVtx==%i", i+1);
    h_pfRho[i] = fs->make<TH1F>(hname, htit, 200,0.0,20.0);

    //sprintf(hname, "h_totalIsolation_nVtx_%i", i+1);
    //sprintf(hname, "totalIsolation nVtx==%i", i+1);
    //h_totalIsolation[i] = fs->make<TH1F>(hname, htit, 100,0.0,30.0);
    //
    //sprintf(hname, "h_totalIsolationRho_nVtx_%i", i+1);
    //sprintf(hname, "totalIsolationRho nVtx==%i", i+1);
    //h_totalIsolationRho[i] = fs->make<TH1F>(hname, htit, 100,0.0,30.0);

  }

  // photon candidates
  h_PhotonsCand_NPhotons   = fs->make<TH1F>("h_PhotonsCand_NPhotons", "h_PhotonsCand_NPhotons", 10, -0.5, 9.5);
  h_PhotonsCand_NJets_Pt30 = fs->make<TH1F>("h_PhotonsIso_NJets_Pt30", "h_PhotonsIso_NJets_Pt30", 20, -0.5, 19.5);
  h_PhotonsCand_NJets_Pt50Eta25  = fs->make<TH1F>("h_PhotonsCand_NJets_Pt50Eta25", "h_PhotonsCand_NJets_Pt50Eta25", 20, -0.5, 19.5);
  for(int iphot=0; iphot<2; iphot++) {
    sprintf(hname, "h_PhotonsCand_Pt_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_Pt_%i", iphot);
    h_PhotonsCand_Pt[iphot]  = fs->make<TH1F>(hname, htit, 200, 0.0, 2000.0);

    sprintf(hname, "h_PhotonsCand_Eta_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_Eta_%i", iphot);
    h_PhotonsCand_Eta[iphot] = fs->make<TH1F>(hname, htit, 60, -3.0, 3.0);

    sprintf(hname, "h_PhotonsCand_Phi_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_Phi_%i", iphot);
    h_PhotonsCand_Phi[iphot] = fs->make<TH1F>(hname, htit, 64, -3.2, 3.2);
    
    sprintf(hname, "h_PhotonsCand_TrkIso_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_TrkIso_%i", iphot);
    h_PhotonsCand_TrkIso[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);
    
    sprintf(hname, "h_PhotonsCand_EcalIso_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_EcalIso_%i", iphot);
    h_PhotonsCand_EcalIso[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsCand_HcalIso_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_HcalIso_%i", iphot);
    h_PhotonsCand_HcalIso[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsCand_HcalIso2012_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_HcalIso2012_%i", iphot);
    h_PhotonsCand_HcalIso2012[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsCand_PFCharged_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_PFCharged_%i", iphot);
    h_PhotonsCand_PFCharged[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsCand_PFNeutral_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_PFNeutral_%i", iphot);
    h_PhotonsCand_PFNeutral[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsCand_PFGamma_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_PFGamma_%i", iphot);
    h_PhotonsCand_PFGamma[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsCand_PFChargedRel_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_PFChargedRel_%i", iphot);
    h_PhotonsCand_PFChargedRel[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsCand_PFNeutralRel_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_PFNeutralRel_%i", iphot);
    h_PhotonsCand_PFNeutralRel[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsCand_PFGammaRel_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_PFGammaRel_%i", iphot);
    h_PhotonsCand_PFGammaRel[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsCand_PFChargedRelPU_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_PFChargedRelPU_%i", iphot);
    h_PhotonsCand_PFChargedRelPU[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsCand_PFNeutralRelPU_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_PFNeutralRelPU_%i", iphot);
    h_PhotonsCand_PFNeutralRelPU[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsCand_PFGammaRelPU_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_PFGammaRelPU_%i", iphot);
    h_PhotonsCand_PFGammaRelPU[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsCand_SigEtaEta_EB_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_SigEtaEta_EB_%i", iphot);
    h_PhotonsCand_SigEtaEta_EB[iphot] = fs->make<TH1F>(hname, htit, 150, 0.0, 0.05);
    sprintf(hname, "h_PhotonsCand_SigEtaEta_EE_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_SigEtaEta_EE_%i", iphot);
    h_PhotonsCand_SigEtaEta_EE[iphot] = fs->make<TH1F>(hname, htit, 250, 0.0, 0.1);

    sprintf(hname, "h_PhotonsCand_hOverE_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_hOverE_%i", iphot);
    h_PhotonsCand_hOverE[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 1.0);

    sprintf(hname, "h_PhotonsCand_hOverE2012_%i", iphot);
    sprintf(htit,  "h_PhotonsCand_hOverE2012_%i", iphot);
    h_PhotonsCand_hOverE2012[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

  }

  // identified & isolated photons 
  h_PhotonsIso_NPhotons   = fs->make<TH1F>("h_PhotonsIso_NPhotons", "h_PhotonsIso_NPhotons", 10, -0.5, 9.5);
  for(int iphot=0; iphot<2; iphot++) {
    sprintf(hname, "h_PhotonsIso_Pt_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_Pt_%i", iphot);
    h_PhotonsIso_Pt[iphot]  = fs->make<TH1F>(hname, htit, 200, 0.0, 2000.0);

    sprintf(hname, "h_PhotonsIso_Eta_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_Eta_%i", iphot);
    h_PhotonsIso_Eta[iphot] = fs->make<TH1F>(hname, htit, 60, -3.0, 3.0);

    sprintf(hname, "h_PhotonsIso_Phi_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_Phi_%i", iphot);
    h_PhotonsIso_Phi[iphot] = fs->make<TH1F>(hname, htit, 64, -3.2, 3.2);
    
    sprintf(hname, "h_PhotonsIso_TrkIso_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_TrkIso_%i", iphot);
    h_PhotonsIso_TrkIso[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);
    
    sprintf(hname, "h_PhotonsIso_EcalIso_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_EcalIso_%i", iphot);
    h_PhotonsIso_EcalIso[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsIso_HcalIso_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_HcalIso_%i", iphot);
    h_PhotonsIso_HcalIso[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsIso_HcalIso2012_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_HcalIso2012_%i", iphot);
    h_PhotonsIso_HcalIso2012[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsIso_PFCharged_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_PFCharged_%i", iphot);
    h_PhotonsIso_PFCharged[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsIso_PFNeutral_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_PFNeutral_%i", iphot);
    h_PhotonsIso_PFNeutral[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsIso_PFGamma_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_PFGamma_%i", iphot);
    h_PhotonsIso_PFGamma[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_PhotonsIso_PFChargedRel_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_PFChargedRel_%i", iphot);
    h_PhotonsIso_PFChargedRel[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsIso_PFNeutralRel_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_PFNeutralRel_%i", iphot);
    h_PhotonsIso_PFNeutralRel[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsIso_PFGammaRel_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_PFGammaRel_%i", iphot);
    h_PhotonsIso_PFGammaRel[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsIso_PFChargedRelPU_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_PFChargedRelPU_%i", iphot);
    h_PhotonsIso_PFChargedRelPU[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsIso_PFNeutralRelPU_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_PFNeutralRelPU_%i", iphot);
    h_PhotonsIso_PFNeutralRelPU[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsIso_PFGammaRelPU_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_PFGammaRelPU_%i", iphot);
    h_PhotonsIso_PFGammaRelPU[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsIso_SigEtaEta_EB_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_SigEtaEta_EB_%i", iphot);
    h_PhotonsIso_SigEtaEta_EB[iphot] = fs->make<TH1F>(hname, htit, 150, 0.0, 0.05);
    sprintf(hname, "h_PhotonsIso_SigEtaEta_EE_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_SigEtaEta_EE_%i", iphot);
    h_PhotonsIso_SigEtaEta_EE[iphot] = fs->make<TH1F>(hname, htit, 150, 0.0, 0.05);

    sprintf(hname, "h_PhotonsIso_hOverE_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_hOverE_%i", iphot);
    h_PhotonsIso_hOverE[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 1.);

    sprintf(hname, "h_PhotonsIso_hOverE2012_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_hOverE2012_%i", iphot);
    h_PhotonsIso_hOverE2012[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_PhotonsIso_R9_%i", iphot);
    sprintf(htit,  "h_PhotonsIso_R9_%i", iphot);
    h_PhotonsIso_R9[iphot] = fs->make<TH1F>(hname, htit, 50, 0, 2.0);
  }

  h_PhotonsIso_NJets_Pt30       = fs->make<TH1F>("h_PhotonsIso_NJets_Pt30", "h_PhotonsIso_NJets_Pt30", 20, -0.5, 19.5);
  h_PhotonsIso_NJets_Pt50Eta25       = fs->make<TH1F>("h_PhotonsIso_NJets_Pt50Eta25", "h_PhotonsIso_NJets_Pt50Eta25", 20, -0.5, 19.5);
  h_PhotonsIso_HT          = fs->make<TH1F>("h_PhotonsIso_HT",  "h_PhotonsIso_HT",  100, 0, 5000.0);
  h_PhotonsIso_MHT         = fs->make<TH1F>("h_PhotonsIso_MHT", "h_PhotonsIso_MHT", 150, 0, 1500.0);
  h_PhotonsIso_DPhiMHTJet1 = fs->make<TH1F>("h_PhotonsIso_DPhiMHTJet1", "h_PhotonsIso_DPhiMHTJet1", 64, -3.2, 3.2);
  h_PhotonsIso_DPhiMHTJet2 = fs->make<TH1F>("h_PhotonsIso_DPhiMHTJet2", "h_PhotonsIso_DPhiMHTJet2", 64, -3.2, 3.2);
  h_PhotonsIso_DPhiMHTJet3 = fs->make<TH1F>("h_PhotonsIso_DPhiMHTJet3", "h_PhotonsIso_DPhiMHTJet3", 64, -3.2, 3.2);
  
  //h_PhotonIso_PhotPtEta    = fs->make<TH2F>("h_PhotonIso_PhotPtEta",    "h_PhotonIso_PhotPtEta",    NPhotPtBins, -0.5, (double)NPhotPtBins-0.5, NPhotEtaBins, -0.5, (double)NPhotEtaBins-0.5);
  h_PhotonIso_PhotPtEta    = fs->make<TH2F>("h_PhotonIso_PhotPtEta",    "h_PhotonIso_PhotPtEta",    8, -0.5, 8-0.5, 7, -0.5, 7-0.5);
  h_PhotonIso_PhotPtEta->GetXaxis()->SetTitle("Photon Pt");
  h_PhotonIso_PhotPtEta->GetYaxis()->SetTitle("Photon Eta");

  // RA2 cuts in addition to identified & isolated photons 
  h_RA2_Vertices = fs->make<TH1F>("h_RA2_Vertices", "h_RA2_Vertices", 60, 0.0, 60.0);
  h_RA2_VerticesReWeighted = fs->make<TH1F>("h_RA2_VerticesReWeighted", "h_RA2_VerticesReWeighted", 60, 0.0, 60.0);
  h_RA2_NPhotons   = fs->make<TH1F>("h_RA2_NPhotons", "h_RA2_NPhotons", 10, -0.5, 9.5);
  for(int iphot=0; iphot<2; iphot++) {

    sprintf(hname, "h_RA2_EtaVsPhi_%i", iphot);
    sprintf(htit,  "h_RA2_EtaVsPhi_%i", iphot);
    h_RA2_EtaVsPt[iphot] = fs->make<TH2F>(hname, htit,60, -3.0, 3.0, 64, -3.2, 3.2);

    sprintf(hname, "h_RA2_Pt_%i", iphot);
    sprintf(htit,  "h_RA2_Pt_%i", iphot);
    h_RA2_Pt[iphot]  = fs->make<TH1F>(hname, htit, 200, 0.0, 2000.0);

    sprintf(hname, "h_RA2_Eta_%i", iphot);
    sprintf(htit,  "h_RA2_Eta_%i", iphot);
    h_RA2_Eta[iphot] = fs->make<TH1F>(hname, htit, 60, -3.0, 3.0);

    sprintf(hname, "h_RA2_Phi_%i", iphot);
    sprintf(htit,  "h_RA2_Phi_%i", iphot);
    h_RA2_Phi[iphot] = fs->make<TH1F>(hname, htit, 64, -3.2, 3.2);

    sprintf(hname, "h_RA2_TrkIso_%i", iphot);
    sprintf(htit,  "h_RA2_TrkIso_%i", iphot);
    h_RA2_TrkIso[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);
    
    sprintf(hname, "h_RA2_EcalIso_%i", iphot);
    sprintf(htit,  "h_RA2_EcalIso_%i", iphot);
    h_RA2_EcalIso[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_RA2_HcalIso_%i", iphot);
    sprintf(htit,  "h_RA2_HcalIso_%i", iphot);
    h_RA2_HcalIso[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_RA2_HcalIso2012_%i", iphot);
    sprintf(htit,  "h_RA2_HcalIso2012_%i", iphot);
    h_RA2_HcalIso2012[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_RA2_PFCharged_%i", iphot);
    sprintf(htit,  "h_RA2_PFCharged_%i", iphot);
    h_RA2_PFCharged[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_RA2_PFNeutral_%i", iphot);
    sprintf(htit,  "h_RA2_PFNeutral_%i", iphot);
    h_RA2_PFNeutral[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_RA2_PFGamma_%i", iphot);
    sprintf(htit,  "h_RA2_PFGamma_%i", iphot);
    h_RA2_PFGamma[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 10.0);

    sprintf(hname, "h_RA2_PFChargedRel_%i", iphot);
    sprintf(htit,  "h_RA2_PFChargedRel_%i", iphot);
    h_RA2_PFChargedRel[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_RA2_PFNeutralRel_%i", iphot);
    sprintf(htit,  "h_RA2_PFNeutralRel_%i", iphot);
    h_RA2_PFNeutralRel[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_RA2_PFGammaRel_%i", iphot);
    sprintf(htit,  "h_RA2_PFGammaRel_%i", iphot);
    h_RA2_PFGammaRel[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_RA2_PFChargedRelPU_%i", iphot);
    sprintf(htit,  "h_RA2_PFChargedRelPU_%i", iphot);
    h_RA2_PFChargedRelPU[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_RA2_PFNeutralRelPU_%i", iphot);
    sprintf(htit,  "h_RA2_PFNeutralRelPU_%i", iphot);
    h_RA2_PFNeutralRelPU[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_RA2_PFGammaRelPU_%i", iphot);
    sprintf(htit,  "h_RA2_PFGammaRelPU_%i", iphot);
    h_RA2_PFGammaRelPU[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_RA2_SigEtaEta_EB_%i", iphot);
    sprintf(htit,  "h_RA2_SigEtaEta_EB_%i", iphot);
    h_RA2_SigEtaEta_EB[iphot] = fs->make<TH1F>(hname, htit, 150, 0.0, 0.05);
    sprintf(hname, "h_RA2_SigEtaEta_EE_%i", iphot);
    sprintf(htit,  "h_RA2_SigEtaEta_EE_%i", iphot);
    h_RA2_SigEtaEta_EE[iphot] = fs->make<TH1F>(hname, htit, 150, 0.0, 0.05);

    sprintf(hname, "h_RA2_hOverE_%i", iphot);
    sprintf(htit,  "h_RA2_hOverE_%i", iphot);
    h_RA2_hOverE[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 1.0);

    sprintf(hname, "h_RA2_hOverE2012_%i", iphot);
    sprintf(htit,  "h_RA2_hOverE2012_%i", iphot);
    h_RA2_hOverE2012[iphot] = fs->make<TH1F>(hname, htit, 100, 0, 0.5);

    sprintf(hname, "h_RA2_R9_%i", iphot);
    sprintf(htit,  "h_RA2_R9_%i", iphot);
    h_RA2_R9[iphot] = fs->make<TH1F>(hname, htit, 50, 0, 2.0);
  }

  h_drMin_LeadPhotPFJet = fs->make<TH1F>("h_drMin_LeadPhotPFJet", "h_drMin_LeadPhotPFJet", 200, 0.0, 5.0);

  h_RA2_NJets_Pt30       = fs->make<TH1F>("h_RA2_NJets_Pt30", "h_RA2_NJets_Pt30", 20, -0.5, 19.5);
  h_RA2_NJets_Pt50Eta25  = fs->make<TH1F>("h_RA2_NJets_Pt50Eta25", "h_RA2_NJets_Pt50Eta25", 20, -0.5, 19.5);
  h_RA2_HT          = fs->make<TH1F>("h_RA2_HT",  "h_RA2_HT",  100, 0, 5000.0);
  h_RA2_MHT         = fs->make<TH1F>("h_RA2_MHT", "h_RA2_MHT", 150, 0, 1500.0);
  h_RA2_MEff        = fs->make<TH1F>("h_RA2_MEff","h_RA2_MEff",100, 0, 5000.0);
  h_RA2_DPhiMHTJet1 = fs->make<TH1F>("h_RA2_DPhiMHTJet1", "h_RA2_DPhiMHTJet1", 64, -3.2, 3.2);
  h_RA2_DPhiMHTJet2 = fs->make<TH1F>("h_RA2_DPhiMHTJet2", "h_RA2_DPhiMHTJet2", 64, -3.2, 3.2);
  h_RA2_DPhiMHTJet3 = fs->make<TH1F>("h_RA2_DPhiMHTJet3", "h_RA2_DPhiMHTJet3", 64, -3.2, 3.2);
  //h_RA2_PhotPtEta          = fs->make<TH2F>("h_RA2_PhotPtEta",    "h_RA2_PhotPtEta",    NPhotPtBins, -0.5, (double)NPhotPtBins-0.5, NPhotEtaBins, -0.5, (double)NPhotEtaBins-0.5);
  h_RA2_PhotPtEta          = fs->make<TH2F>("h_RA2_PhotPtEta",    "h_RA2_PhotPtEta",    8, -0.5, 8-0.5, 7, -0.5, 7-0.5);
  h_RA2_PhotPtEta->GetXaxis()->SetTitle("Photon Pt");
  h_RA2_PhotPtEta->GetYaxis()->SetTitle("Photon Eta");

  h_preRA2_DPhiMHTJet1 = fs->make<TH1F>("h_preRA2_DPhiMHTJet1", "h_preRA2_DPhiMHTJet1(all ra2 except dphi)", 64, -3.2, 3.2);
  h_preRA2_DPhiMHTJet2 = fs->make<TH1F>("h_preRA2_DPhiMHTJet2", "h_preRA2_DPhiMHTJet2(all ra2 except dphi)", 64, -3.2, 3.2);
  h_preRA2_DPhiMHTJet3 = fs->make<TH1F>("h_preRA2_DPhiMHTJet3", "h_preRA2_DPhiMHTJet3(all ra2 except dphi)", 64, -3.2, 3.2);

  // event yield in exclusive bins
  h_RA2_HTMHT_Excl = fs->make<TH2F>("h_RA2_HTMHT_Excl", "h_RA2_HTMHT_Excl", NRA2HTBins-1, RA2HTVal, NRA2MHTBins-1, RA2MHTVal);
  h_RA2_HTMHT_Excl ->Sumw2();
  h_RA2_HTMHT_Excl->GetXaxis()->SetTitle("HT");
  h_RA2_HTMHT_Excl->GetYaxis()->SetTitle("MHT");

  for(int iht=1; iht<=h_RA2_HTMHT_Excl->GetNbinsX(); iht++){
    for(int imht=1; imht<=h_RA2_HTMHT_Excl->GetNbinsY(); imht++) {
      int htxL  = (int)h_RA2_HTMHT_Excl->GetXaxis()->GetBinLowEdge(iht);
      int htxU  = (int)h_RA2_HTMHT_Excl->GetXaxis()->GetBinUpEdge(iht);

      int mhtyL = (int)h_RA2_HTMHT_Excl->GetYaxis()->GetBinLowEdge(imht);     
      int mhtyU = (int)h_RA2_HTMHT_Excl->GetYaxis()->GetBinUpEdge(imht);     

      sprintf(hname, "h_RA2_PhotPtEta_Excl_%iHT%i_%iMHT%i", htxL,htxU, mhtyL,mhtyU);
      //std::cout <<iht <<" "<< imht << " "<<  hname << std::endl;
      sprintf(htit,  "h_RA2_PhotPtEta_Excl_%iHT%i_%iMHT%i", htxL,htxU, mhtyL,mhtyU);
      h_RA2_PhotPtEta_Excl[iht-1][imht-1]  = fs->make<TH2F>(hname, htit, NPhotEtaBins-1, PhotEtaVal, NPhotPtBins-1, PhotPtVal);
    }
  }

  // event yield in inclusive bins
  h_RA2_HTMHT_Incl = fs->make<TH2F>("h_RA2_HTMHT_Incl", "h_RA2_HTMHT_Incl", NHTBins-1, HTVal, NMHTBins-1, MHTVal);
  h_RA2_HTMHT_Incl ->Sumw2();
  h_RA2_HTMHT_Incl->GetXaxis()->SetTitle("HT");
  h_RA2_HTMHT_Incl->GetYaxis()->SetTitle("MHT");

  sprintf(hname, "h_RA2_NJets_HT350_MHT200");
  sprintf(htit,  "NJets: RA2 Baseline (HT>350 & MHT>200 GeV)");
  h_RA2_NJets_HT350_MHT200 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT350_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_NJets_HT800_MHT200");
  sprintf(htit,  "NJets: RA2 High HT (HT>800 & MHT>200 GeV)");
  h_RA2_NJets_HT800_MHT200 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT800_MHT200 ->Sumw2();

  sprintf(hname, "h_RA2_NJets_HT800_MHT500");
  sprintf(htit,  "NJets: RA2 High MHT (HT>800 & MHT>500 GeV)");
  h_RA2_NJets_HT800_MHT500 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT800_MHT500 ->Sumw2();

  sprintf(hname, "h_RA2_NJets_HT500_MHT350");
  sprintf(htit,  "NJets: RA2 Medium HT,MHT (HT>500 & MHT>350 GeV)");
  h_RA2_NJets_HT500_MHT350 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT500_MHT350 ->Sumw2();

  sprintf(hname, "h_RA2_NJets_HT500_MHT200");
  sprintf(htit,  "NJets: RA2 New Baseline (HT>500 & MHT>200 GeV)");
  h_RA2_NJets_HT500_MHT200 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT500_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_NJets_HT1000_MHT600");
  sprintf(htit,  "NJets: RA2 New High MHT (HT>1000 & MHT>600 GeV)");
  h_RA2_NJets_HT1000_MHT600 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT1000_MHT600 ->Sumw2();

  sprintf(hname, "h_RA2_NJets_HT1200_MHT400");
  sprintf(htit,  "NJets: RA2 New High HT (HT>1200 & MHT>400 GeV)");
  h_RA2_NJets_HT1200_MHT400 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT1200_MHT400 ->Sumw2();


  sprintf(hname, "h_RA2_HT_HT350_MHT200");
  sprintf(htit,  "HT: RA2 Baseline (HT>350 & MHT>200 GeV)");
  h_RA2_HT_HT350_MHT200 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT350_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_HT_HT800_MHT200");
  sprintf(htit,  "HT: RA2 High HT (HT>800 & MHT>200 GeV)");
  h_RA2_HT_HT800_MHT200 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT800_MHT200 ->Sumw2();

  sprintf(hname, "h_RA2_HT_HT800_MHT500");
  sprintf(htit,  "HT: RA2 High MHT (HT>800 & MHT>500 GeV)");
  h_RA2_HT_HT800_MHT500 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT800_MHT500 ->Sumw2();

  sprintf(hname, "h_RA2_HT_HT500_MHT350");
  sprintf(htit,  "HT: RA2 Medium HT,MHT (HT>500 & MHT>350 GeV)");
  h_RA2_HT_HT500_MHT350 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT500_MHT350 ->Sumw2();

  sprintf(hname, "h_RA2_HT_HT500_MHT200");
  sprintf(htit,  "HT: RA2 New Baseline (HT>500 & MHT>200 GeV)");
  h_RA2_HT_HT500_MHT200 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT500_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_HT_HT1000_MHT600");
  sprintf(htit,  "HT: RA2 New High MHT (HT>1000 & MHT>600 GeV)");
  h_RA2_HT_HT1000_MHT600 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT1000_MHT600 ->Sumw2();

  sprintf(hname, "h_RA2_HT_HT1200_MHT400");
  sprintf(htit,  "HT: RA2 New High HT (HT>1200 & MHT>400 GeV)");
  h_RA2_HT_HT1200_MHT400 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT1200_MHT400 ->Sumw2();

  sprintf(hname, "h_RA2_MHT_HT350_MHT200");
  sprintf(htit,  "MHT: RA2 Baseline (HT>350 & MHT>200 GeV)");
  h_RA2_MHT_HT350_MHT200 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT350_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_MHT_HT800_MHT200");
  sprintf(htit,  "MHT: RA2 High HT (HT>800 & MHT>200 GeV)");
  h_RA2_MHT_HT800_MHT200 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT800_MHT200 ->Sumw2();

  sprintf(hname, "h_RA2_MHT_HT800_MHT500");
  sprintf(htit,  "MHT: RA2 High MHT (HT>800 & MHT>500 GeV)");
  h_RA2_MHT_HT800_MHT500 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT800_MHT500 ->Sumw2();

  sprintf(hname, "h_RA2_MHT_HT500_MHT350");
  sprintf(htit,  "MHT: RA2 Medium HT,MHT (HT>500 & MHT>350 GeV)");
  h_RA2_MHT_HT500_MHT350 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT500_MHT350 ->Sumw2();

  sprintf(hname, "h_RA2_MHT_HT500_MHT200");
  sprintf(htit,  "MHT: RA2 New Baseline (HT>500 & MHT>200 GeV)");
  h_RA2_MHT_HT500_MHT200 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT500_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_MHT_HT1000_MHT600");
  sprintf(htit,  "MHT: RA2 New High MHT (HT>1000 & MHT>600 GeV)");
  h_RA2_MHT_HT1000_MHT600 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT1000_MHT600 ->Sumw2();

  sprintf(hname, "h_RA2_MHT_HT1200_MHT400");
  sprintf(htit,  "MHT: RA2 New High HT (HT>1200 & MHT>400 GeV)");
  h_RA2_MHT_HT1200_MHT400 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT1200_MHT400 ->Sumw2();

  sprintf(hname, "h_RA2_MEff_HT350_MHT200");
  sprintf(htit,  "MEff: RA2 Baseline (HT>350 & MHT>200 GeV)");
  h_RA2_MEff_HT350_MHT200 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT350_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_MEff_HT800_MHT200");
  sprintf(htit,  "MEff: RA2 High HT (HT>800 & MHT>200 GeV)");
  h_RA2_MEff_HT800_MHT200 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT800_MHT200 ->Sumw2();

  sprintf(hname, "h_RA2_MEff_HT800_MHT500");
  sprintf(htit,  "MEff: RA2 High MHT (HT>800 & MHT>500 GeV)");
  h_RA2_MEff_HT800_MHT500 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT800_MHT500 ->Sumw2();

  sprintf(hname, "h_RA2_MEff_HT500_MHT350");
  sprintf(htit,  "MEff: RA2 Medium HT,MHT (HT>500 & MHT>350 GeV)");
  h_RA2_MEff_HT500_MHT350 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT500_MHT350 ->Sumw2();

  sprintf(hname, "h_RA2_MEff_HT500_MHT200");
  sprintf(htit,  "MEff: RA2 New Baseline (HT>500 & MHT>200 GeV)");
  h_RA2_MEff_HT500_MHT200 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT500_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_MEff_HT1000_MHT600");
  sprintf(htit,  "MEff: RA2 New High MHT (HT>1000 & MHT>600 GeV)");
  h_RA2_MEff_HT1000_MHT600 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT1000_MHT600 ->Sumw2();

  sprintf(hname, "h_RA2_MEff_HT1200_MHT400");
  sprintf(htit,  "MEff: RA2 New High HT (HT>1200 & MHT>400 GeV)");
  h_RA2_MEff_HT1200_MHT400 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT1200_MHT400 ->Sumw2();

  sprintf(hname, "h_RA2_PhotPtEta_HT350_MHT200");
  sprintf(htit,  "Photon (Pt,Eta): RA2 Baseline (HT>350 & MHT>200 GeV)");
  h_RA2_PhotPtEta_HT350_MHT200 = fs->make<TH2F>(hname, htit, NPhotEtaBins-1, PhotEtaVal, NPhotPtBins-1, PhotPtVal);
  h_RA2_PhotPtEta_HT350_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_PhotPtEta_HT800_MHT200");
  sprintf(htit,  "Photon (Pt,Eta): RA2 High HT (HT>800 & MHT>200 GeV)");
  h_RA2_PhotPtEta_HT800_MHT200 = fs->make<TH2F>(hname, htit, NPhotEtaBins-1, PhotEtaVal, NPhotPtBins-1, PhotPtVal);
  h_RA2_PhotPtEta_HT800_MHT200 ->Sumw2();

  sprintf(hname, "h_RA2_PhotPtEta_HT800_MHT500");
  sprintf(htit,  "Photon (Pt,Eta): RA2 High MHT (HT>800 & MHT>500 GeV)");
  h_RA2_PhotPtEta_HT800_MHT500 = fs->make<TH2F>(hname, htit, NPhotEtaBins-1, PhotEtaVal, NPhotPtBins-1, PhotPtVal);
  h_RA2_PhotPtEta_HT800_MHT500 ->Sumw2();

  sprintf(hname, "h_RA2_PhotPtEta_HT500_MHT350");
  sprintf(htit,  "Photon (Pt,Eta): RA2 Medium HT,MHT (HT>500 & MHT>350 GeV)");
  h_RA2_PhotPtEta_HT500_MHT350 = fs->make<TH2F>(hname, htit, NPhotEtaBins-1, PhotEtaVal, NPhotPtBins-1, PhotPtVal);
  h_RA2_PhotPtEta_HT500_MHT350 ->Sumw2();

  sprintf(hname, "h_RA2_PhotPtEta_HT500_MHT200");
  sprintf(htit,  "Photon (Pt,Eta): RA2 New Baseline (HT>500 & MHT>200 GeV)");
  h_RA2_PhotPtEta_HT500_MHT200 = fs->make<TH2F>(hname, htit, NPhotEtaBins-1, PhotEtaVal, NPhotPtBins-1, PhotPtVal);
  h_RA2_PhotPtEta_HT500_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_PhotPtEta_HT1000_MHT600");
  sprintf(htit,  "Photon (Pt,Eta): RA2 New High MHT (HT>1000 & MHT>600 GeV)");
  h_RA2_PhotPtEta_HT1000_MHT600 = fs->make<TH2F>(hname, htit, NPhotEtaBins-1, PhotEtaVal, NPhotPtBins-1, PhotPtVal);
  h_RA2_PhotPtEta_HT1000_MHT600 ->Sumw2();

  sprintf(hname, "h_RA2_PhotPtEta_HT1200_MHT400");
  sprintf(htit,  "Photon (Pt,Eta): RA2 New High HT (HT>1200 & MHT>400 GeV)");
  h_RA2_PhotPtEta_HT1200_MHT400 = fs->make<TH2F>(hname, htit, NPhotEtaBins-1, PhotEtaVal, NPhotPtBins-1, PhotPtVal);
  h_RA2_PhotPtEta_HT1200_MHT400 ->Sumw2();


  if(doGenAnalysis_) {
    h_GenBoson_preRA2_dPhiMHTJet1 = fs->make<TH1F>("h_GenBoson_preRA2_dPhiMHTJet1", "h_GenBoson_preRA2_dPhiMHTJet1", 128, -3.2, 3.2);
    h_GenBoson_preRA2_dPhiMHTJet2 = fs->make<TH1F>("h_GenBoson_preRA2_dPhiMHTJet2", "h_GenBoson_preRA2_dPhiMHTJet2", 128, -3.2, 3.2);
    h_GenBoson_preRA2_dPhiMHTJet3 = fs->make<TH1F>("h_GenBoson_preRA2_dPhiMHTJet3", "h_GenBoson_preRA2_dPhiMHTJet3", 128, -3.2, 3.2);

    
    h_GenBoson_PtVsHT  = fs->make<TH2F>("h_GenBoson_PtVsHT",  "GenBoson_Pt  Vs HTJet30",  500, 0.0, 5000.0, 100, 0.0, 1000.0);
    h_GenBoson_PtVsHT->Sumw2();
    h_GenBoson_MHTVsHT = fs->make<TH2F>("h_GenBoson_MHTVsHT", "GenBoson_MHT Vs HTJet30",  500, 0.0, 5000.0, 100, 0.0, 1000.0);
    h_GenBoson_MHTVsHT->Sumw2();
    h_GenBoson_NJets_Pt30       = fs->make<TH1F>("h_GenBoson_NJets_Pt30",      "h_GenBoson_NJets_Pt30", 25, -0.5, 24.5);
    h_GenBoson_NJets_Pt30->Sumw2();
    h_GenBoson_NJets_Pt50Eta25  = fs->make<TH1F>("h_GenBoson_NJets_Pt50Eta25", "h_GenBoson_NJets_Pt50Eta25", 25, -0.5, 24.5);
    h_GenBoson_NJets_Pt50Eta25->Sumw2();
    h_GenBoson_Pt        = fs->make<TH1F>("h_GenBoson_Pt",         "Pt of Photon",            100, 0.0, 1000.0);
    h_GenBoson_Pt->Sumw2();
    h_GenBoson_Eta       = fs->make<TH1F>("h_GenBoson_Eta",        "Eta of Photon",          100, -5.0, 5.0);
    h_GenBoson_Eta->Sumw2();
    h_GenBoson_Phi       = fs->make<TH1F>("h_GenBoson_Phi",        "Phi of Photon",           64, -3.2, 3.2);
    h_GenBoson_Phi->Sumw2();
    h_GenBoson_HT        = fs->make<TH1F>("h_GenBoson_HT",         "HT of Photon",            500, 0.0, 5000.0);
    h_GenBoson_HT->Sumw2();
    h_GenBoson_MHT        = fs->make<TH1F>("h_GenBoson_MHT",       "MHT of Photon",           100, 0.0, 1000.0);
    h_GenBoson_MHT->Sumw2();

    h_GenBoson_RA2_Pt        = fs->make<TH1F>("h_GenBoson_RA2_Pt",         "Pt of Photon",            100, 0.0, 1000.0);
    h_GenBoson_RA2_Pt->Sumw2();
    h_GenBoson_RA2_Eta       = fs->make<TH1F>("h_GenBoson_RA2_Eta",        "Eta of Photon",          100, -5.0, 5.0);
    h_GenBoson_RA2_Eta->Sumw2();
    h_GenBoson_RA2_Phi       = fs->make<TH1F>("h_GenBoson_RA2_Phi",        "Phi of Photon",           64, -3.2, 3.2);
    h_GenBoson_RA2_Phi->Sumw2();
    h_GenBoson_RA2_HT        = fs->make<TH1F>("h_GenBoson_RA2_HT",         "HT of Photon",            500, 0.0, 5000.0);
    h_GenBoson_RA2_HT->Sumw2();
    h_GenBoson_RA2_MHT        = fs->make<TH1F>("h_GenBoson_RA2_MHT",       "MHT of Photon",           100, 0.0, 1000.0);
    h_GenBoson_RA2_MHT->Sumw2();

    // event yield in exclusive bins
    h_GenBosonTotal_RA2_HTMHT_Excl = fs->make<TH2F>("h_GenBosonTotal_RA2_HTMHT_Excl", "h_GenBosonTotal_RA2_HTMHT_Excl", NRA2HTBins-1, RA2HTVal, NRA2MHTBins-1, RA2MHTVal);
    h_GenBosonTotal_RA2_HTMHT_Excl ->Sumw2();
    h_GenBosonTotal_RA2_HTMHT_Excl->GetXaxis()->SetTitle("HT");
    h_GenBosonTotal_RA2_HTMHT_Excl->GetYaxis()->SetTitle("MHT");

    h_GenBosonKineAcc_RA2_HTMHT_Excl = fs->make<TH2F>("h_GenBosonKineAcc_RA2_HTMHT_Excl", "h_GenBosonKineAcc_RA2_HTMHT_Excl", NRA2HTBins-1, RA2HTVal, NRA2MHTBins-1, RA2MHTVal);
    h_GenBosonKineAcc_RA2_HTMHT_Excl ->Sumw2();
    h_GenBosonKineAcc_RA2_HTMHT_Excl->GetXaxis()->SetTitle("HT");
    h_GenBosonKineAcc_RA2_HTMHT_Excl->GetYaxis()->SetTitle("MHT");

    h_GenBosonMatched_RA2_HTMHT_Excl = fs->make<TH2F>("h_GenBosonMatched_RA2_HTMHT_Excl", "h_GenBosonMatched_RA2_HTMHT_Excl", NRA2HTBins-1, RA2HTVal, NRA2MHTBins-1, RA2MHTVal);
    h_GenBosonMatched_RA2_HTMHT_Excl ->Sumw2();
    h_GenBosonMatched_RA2_HTMHT_Excl->GetXaxis()->SetTitle("HT");
    h_GenBosonMatched_RA2_HTMHT_Excl->GetYaxis()->SetTitle("MHT");
    
    // event yield in inclusive bins
    h_GenBosonTotal_RA2_HTMHT_Incl = fs->make<TH2F>("h_GenBosonTotal_RA2_HTMHT_Incl", "h_GenBosonTotal_RA2_HTMHT_Incl", NHTBins-1, HTVal, NMHTBins-1, MHTVal);
    h_GenBosonTotal_RA2_HTMHT_Incl ->Sumw2();
    h_GenBosonTotal_RA2_HTMHT_Incl->GetXaxis()->SetTitle("HT");
    h_GenBosonTotal_RA2_HTMHT_Incl->GetYaxis()->SetTitle("MHT");

    h_GenBosonKineAcc_RA2_HTMHT_Incl = fs->make<TH2F>("h_GenBosonKineAcc_RA2_HTMHT_Incl", "h_GenBosonKineAcc_RA2_HTMHT_Incl", NHTBins-1, HTVal, NMHTBins-1, MHTVal);
    h_GenBosonKineAcc_RA2_HTMHT_Incl ->Sumw2();
    h_GenBosonKineAcc_RA2_HTMHT_Incl->GetXaxis()->SetTitle("HT");
    h_GenBosonKineAcc_RA2_HTMHT_Incl->GetYaxis()->SetTitle("MHT");

    h_GenBosonMatched_RA2_HTMHT_Incl = fs->make<TH2F>("h_GenBosonMatched_RA2_HTMHT_Incl", "h_GenBosonMatched_RA2_HTMHT_Incl", NHTBins-1, HTVal, NMHTBins-1, MHTVal);
    h_GenBosonMatched_RA2_HTMHT_Incl ->Sumw2();
    h_GenBosonMatched_RA2_HTMHT_Incl->GetXaxis()->SetTitle("HT");
    h_GenBosonMatched_RA2_HTMHT_Incl->GetYaxis()->SetTitle("MHT");
  }

}


void  RA2ZInvPhotonAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RA2ZInvPhotonAnalyzer);

