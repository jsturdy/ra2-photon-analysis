// -*- C++ -*-
//
// Package:    RA2ZInvDiMuonTreeMaker
// Class:      RA2ZInvDiMuonTreeMaker
// 
/**\class RA2ZInvDiMuonTreeMaker RA2ZInvDiMuonTreeMaker.cc SusyAnalysis/RA2ZInvDiMuonTreeMaker/src/RA2ZInvDiMuonTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seema Sharma
//         Created:  Mon Jun 20 12:58:08 CDT 2011
// $Id: RA2ZInvDiMuonTreeMaker.cc,v 1.9 2013/01/20 02:41:21 sturdy Exp $
//
//


// system include files
#include <memory>
#include <iomanip>
#include <cmath>

#include "ZInvisibleBkgds/Photons/interface/RA2ZInvDiMuonTreeMaker.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include <DataFormats/METReco/interface/MET.h>


RA2ZInvDiMuonTreeMaker::RA2ZInvDiMuonTreeMaker(const edm::ParameterSet& pset) {

  // read parameters from config file
  debug_          = pset.getParameter<bool>("Debug");
  data_           = pset.getParameter<bool>("Data");
  scale_          = pset.getParameter<double>("ScaleFactor");
  muonSrc_        = pset.getParameter<edm::InputTag>("MuonSrc");
  vertexSrc_      = pset.getParameter<edm::InputTag>("VertexSrc");
  jetSrc_         = pset.getParameter<edm::InputTag>("JetSrc");
  bJetSrc_        = pset.getParameter<edm::InputTag>("bJetSrc");
  htJetSrc_       = pset.getParameter<edm::InputTag>("htJetSrc");
  htSrc_          = pset.getParameter<edm::InputTag>("htSource");
  mhtSrc_         = pset.getParameter<edm::InputTag>("mhtSource");
  metSrc_         = pset.getParameter<edm::InputTag>("metSource");
  ra2HTSrc_       = pset.getParameter<edm::InputTag>("ra2HTSource");
  ra2MHTSrc_      = pset.getParameter<edm::InputTag>("ra2MHTSource");
  ra2METSrc_      = pset.getParameter<edm::InputTag>("ra2METSource");
  triggerResults_ = pset.getParameter<edm::InputTag>("TriggerResults");
  runTopTagger_   = pset.getParameter<bool>("runTopTagger");
  looseTopTaggerSrc_   = pset.getParameter<std::string>("looseTopTaggerSource");
  nominalTopTaggerSrc_ = pset.getParameter<std::string>("nominalTopTaggerSource");
  doPUReWeight_   = pset.getParameter<bool>("DoPUReweight");
  puWeightSrc_    = pset.getParameter<edm::InputTag>("PUWeightSource");
  eventWeightSrc_ = pset.getParameter<edm::InputTag >( "EventWeightSource" );

  storeExtraVetos_ = pset.getParameter<bool>("storeExtraVetos");
  ra2ElectronSrc_  = pset.getParameter<edm::InputTag>("ra2ElectronForVeto");
  ra2MuonSrc_      = pset.getParameter<edm::InputTag>("ra2MuonForVeto");
  muonVetoSrc_     = pset.getParameter<edm::InputTag>("muonVetoSource");
  electronVetoSrc_ = pset.getParameter<edm::InputTag>("electronVetoSource");
  isoTrkVetoSrc_   = pset.getParameter<edm::InputTag>("isoTrkVetoSource");
  
  //key to help getting the hlt process from event provenance
  getHLTfromConfig_ = false;
  checkedProcess_ = false;
  processName_    = "";
  reducedValues = 0; 

}

RA2ZInvDiMuonTreeMaker::~RA2ZInvDiMuonTreeMaker() {
  delete reducedValues;
  reducedValues = 0; 
}

void RA2ZInvDiMuonTreeMaker::analyze(const edm::Event& ev, const edm::EventSetup& es) {

  using namespace edm;

  // get event-id
  unsigned int event = (ev.id()).event();
  unsigned int run   = (ev.id()).run();
  unsigned int lumi  =  ev.luminosityBlock();

  m_event = event;
  m_run   = run;
  m_lumi  = lumi;

  // get muons 
  edm::Handle< std::vector<pat::Muon> > patMuons;
  ev.getByLabel(muonSrc_, patMuons); 
  if (debug_)
    std::cout<<"Muon collection has size "<<patMuons->size()<<std::endl;

  edm::Handle< std::vector<pat::Electron> > ra2PATElectrons;
  ev.getByLabel(ra2ElectronSrc_, ra2PATElectrons); 
  m_passRA2ElVeto = true;
  if (ra2PATElectrons->size()>0)
    m_passRA2ElVeto = false;

  edm::Handle< std::vector<pat::Muon> > ra2PATMuons;
  ev.getByLabel(ra2MuonSrc_, ra2PATMuons); 
  m_passRA2MuVeto = true;
  if (ra2PATMuons->size()>2)
    m_passRA2MuVeto = false;

  edm::Handle<reco::GenParticleCollection> gens;
  if (!data_)
    ev.getByLabel("genParticles",gens);

  // get jetcollection
  edm::Handle<reco::VertexCollection > vertices;
  ev.getByLabel(vertexSrc_, vertices);

  edm::Handle<edm::View<pat::Jet> > jets;
  ev.getByLabel(jetSrc_, jets);

  edm::Handle<edm::View<pat::Jet> > htJets;
  ev.getByLabel(htJetSrc_, htJets);

  edm::Handle<edm::View<pat::Jet> > bJets;
  ev.getByLabel(bJetSrc_, bJets);

  edm::Handle<double > ht;
  ev.getByLabel(htSrc_, ht);

  edm::Handle<edm::View<reco::MET> > mht;
  ev.getByLabel(mhtSrc_, mht);

  edm::Handle<edm::View<reco::MET> > met;
  ev.getByLabel(metSrc_, met);

  edm::Handle<double > ra2HT;
  ev.getByLabel(ra2HTSrc_, ra2HT);

  edm::Handle<edm::View<reco::MET> > ra2MHT;
  ev.getByLabel(ra2MHTSrc_, ra2MHT);

  edm::Handle<edm::View<reco::MET> > ra2MET;
  ev.getByLabel(ra2METSrc_, ra2MET);

  ///top tagger variables
  edm::Handle<double> hLoosebestTopJetMass;
  edm::Handle<double> hLooseMTbJet;
  edm::Handle<double> hLooseMTbestTopJet;
  edm::Handle<double> hLooseMT2;
  edm::Handle<double> hLooseMTbestWJet;
  edm::Handle<double> hLooseMTbestbJet;
  edm::Handle<double> hLooseMTremainingTopJet;
  edm::Handle<double> hLooselinearCombMTbJetPlusMTbestTopJet;

  edm::Handle<bool> hLooseremainPassCSVS;

  edm::Handle<int> hLoosebestTopJetIdx;
  edm::Handle<int> hLoosepickedRemainingCombfatJetIdx;

  edm::Handle<double> hNominalbestTopJetMass;
  edm::Handle<double> hNominalMTbJet;
  edm::Handle<double> hNominalMTbestTopJet;
  edm::Handle<double> hNominalMT2;

  edm::Handle<double> hNominalMTbestWJet;
  edm::Handle<double> hNominalMTbestbJet;
  edm::Handle<double> hNominalMTremainingTopJet;
  edm::Handle<double> hNominallinearCombMTbJetPlusMTbestTopJet;

  edm::Handle<bool> hNominalremainPassCSVS;

  edm::Handle<int> hNominalbestTopJetIdx;
  edm::Handle<int> hNominalpickedRemainingCombfatJetIdx;

  if (runTopTagger_) {
    ev.getByLabel(looseTopTaggerSrc_,"bestTopJetMass", hLoosebestTopJetMass);
    ev.getByLabel(looseTopTaggerSrc_,"mTbJet"        , hLooseMTbJet);
    ev.getByLabel(looseTopTaggerSrc_,"mTbestTopJet"  , hLooseMTbestTopJet);
    ev.getByLabel(looseTopTaggerSrc_,"MT2"           , hLooseMT2);
    ev.getByLabel(looseTopTaggerSrc_,"bestTopJetIdx"                   , hLoosebestTopJetIdx);
    ev.getByLabel(looseTopTaggerSrc_,"remainPassCSVS"                  , hLooseremainPassCSVS);
    ev.getByLabel(looseTopTaggerSrc_,"pickedRemainingCombfatJetIdx"    , hLoosepickedRemainingCombfatJetIdx);
    ev.getByLabel(looseTopTaggerSrc_,"mTbestWJet"                      , hLooseMTbestWJet);
    ev.getByLabel(looseTopTaggerSrc_,"mTbestbJet"                      , hLooseMTbestbJet);
    ev.getByLabel(looseTopTaggerSrc_,"mTremainingTopJet"               , hLooseMTremainingTopJet);
    ev.getByLabel(looseTopTaggerSrc_,"linearCombmTbJetPlusmTbestTopJet", hLooselinearCombMTbJetPlusMTbestTopJet);

    ev.getByLabel(nominalTopTaggerSrc_,"bestTopJetMass", hNominalbestTopJetMass);
    ev.getByLabel(nominalTopTaggerSrc_,"mTbJet"        , hNominalMTbJet);
    ev.getByLabel(nominalTopTaggerSrc_,"mTbestTopJet"  , hNominalMTbestTopJet);
    ev.getByLabel(nominalTopTaggerSrc_,"MT2"           , hNominalMT2);
    ev.getByLabel(nominalTopTaggerSrc_,"bestTopJetIdx"                   , hNominalbestTopJetIdx);
    ev.getByLabel(nominalTopTaggerSrc_,"remainPassCSVS"                  , hNominalremainPassCSVS);
    ev.getByLabel(nominalTopTaggerSrc_,"pickedRemainingCombfatJetIdx"    , hNominalpickedRemainingCombfatJetIdx);
    ev.getByLabel(nominalTopTaggerSrc_,"mTbestWJet"                      , hNominalMTbestWJet);
    ev.getByLabel(nominalTopTaggerSrc_,"mTbestbJet"                      , hNominalMTbestbJet);
    ev.getByLabel(nominalTopTaggerSrc_,"mTremainingTopJet"               , hNominalMTremainingTopJet);
    ev.getByLabel(nominalTopTaggerSrc_,"linearCombmTbJetPlusmTbestTopJet", hNominallinearCombMTbJetPlusMTbestTopJet);
  }

  if (debug_) {
    std::cout<<"vertex collection has size "<<vertices->size()<<std::endl;
    std::cout<<"Jet collection has size "<<jets->size()<<std::endl;
    std::cout<<"HT Jet collection has size "<<htJets->size()<<std::endl;
    std::cout<<"b-Jet collection has size "<<bJets->size()<<std::endl;
    std::cout<<"HT value "<<*ht<<std::endl;
    std::cout<<"MHT value "<<(*mht)[0].pt()<<std::endl;
    std::cout<<"MET value "<<(*met)[0].pt()<<std::endl;

    if (runTopTagger_) {
      std::cout<<"hLoosebestTopJetMass value "<<*hLoosebestTopJetMass<<std::endl;
      std::cout<<"hLooseMTbJet value "        <<*hLooseMTbJet<<std::endl;
      std::cout<<"hLooseMTbestTopJet value "  <<*hLooseMTbestTopJet<<std::endl;
      std::cout<<"hLooseMT2 value "           <<*hLooseMT2<<std::endl;
      std::cout<<"hLoosebestTopJetIdx"                   <<*hLoosebestTopJetIdx<<std::endl;
      std::cout<<"hLooseremainPassCSVS"                  <<*hLooseremainPassCSVS<<std::endl;
      std::cout<<"hLoosepickedRemainingCombfatJetIdx"    <<*hLoosepickedRemainingCombfatJetIdx<<std::endl;
      std::cout<<"hLooseMTbestWJet"                      <<*hLooseMTbestWJet<<std::endl;
      std::cout<<"hLooseMTbestbJet"                      <<*hLooseMTbestbJet<<std::endl;
      std::cout<<"hLooseMTremainingTopJet"               <<*hLooseMTremainingTopJet<<std::endl;
      std::cout<<"hLooselinearCombMTbJetPlusMTbestTopJet"<<*hLooselinearCombMTbJetPlusMTbestTopJet<<std::endl;
      
      std::cout<<"hNominalbestTopJetMass value "<<*hNominalbestTopJetMass<<std::endl;
      std::cout<<"hNominalMTbJet value "        <<*hNominalMTbJet<<std::endl;
      std::cout<<"hNominalMTbestTopJet value "  <<*hNominalMTbestTopJet<<std::endl;
      std::cout<<"hNominalMT2 value "           <<*hNominalMT2<<std::endl;
      std::cout<<"hNominalbestTopJetIdx"                   <<*hNominalbestTopJetIdx<<std::endl;
      std::cout<<"hNominalremainPassCSVS"                  <<*hNominalremainPassCSVS<<std::endl;
      std::cout<<"hNominalpickedRemainingCombfatJetIdx"    <<*hNominalpickedRemainingCombfatJetIdx<<std::endl;
      std::cout<<"hNominalMTbestWJet"                      <<*hNominalMTbestWJet<<std::endl;
      std::cout<<"hNominalMTbestbJet"                      <<*hNominalMTbestbJet<<std::endl;
      std::cout<<"hNominalMTremainingTopJet"               <<*hNominalMTremainingTopJet<<std::endl;
      std::cout<<"hNominallinearCombMTbJetPlusMTbestTopJet"<<*hNominallinearCombMTbJetPlusMTbestTopJet<<std::endl;
    }
  }

  // if MC, do PU reweighting
  double pu_event_wt = 1.0;
  edm::Handle<double> puweight;
  if( doPUReWeight_ ) {
    ev.getByLabel(puWeightSrc_, puweight);
    pu_event_wt = *puweight;
  }

  double event_wt = 1.;
  edm::Handle<double> eventWeight;
  ev.getByLabel(eventWeightSrc_,eventWeight);
  event_wt = *eventWeight;

  //m_EventWt = scale_;
  m_EventWt = event_wt;
  m_PUWt     = pu_event_wt;
  m_Vertices = vertices->size();

  edm::Handle<bool> elVeto;
  edm::Handle<bool> muVeto;
  edm::Handle<bool> isoTrkVeto;

  if (storeExtraVetos_) {
    ev.getByLabel(electronVetoSrc_,elVeto);
    m_passDirIsoElVeto = *elVeto;
    ev.getByLabel(muonVetoSrc_,muVeto);
    m_passDirIsoMuVeto = *muVeto;
    ev.getByLabel(isoTrkVetoSrc_,isoTrkVeto);
    m_passIsoTrkVeto = *isoTrkVeto;
  }
  //////
  if(debug_ && patMuons->size() > 0) {
    std::cout << "Isolated muons : " << std::endl;
    for( unsigned int imuon=0; imuon<patMuons->size(); imuon++) {
      std::cout << imuon << " " <<(*patMuons)[imuon].pt() 
		<< " " << (*patMuons)[imuon].eta()
		<< " " << (*patMuons)[imuon].phi() 
		<< std::endl;
    }
  }
  //////
  if (patMuons->size() < 2)
    return;

  m_nMuonsIso = patMuons->size();
  m_Muon1Pt   = (*patMuons)[0].pt();
  m_Muon1Eta  = (*patMuons)[0].eta();
  m_Muon1M    = (*patMuons)[0].mass();
  m_Muon1MinDR  = 10.;
  m_Muon1DRJet1 = 10.;
  m_Muon2Pt   = (*patMuons)[1].pt();
  m_Muon2Eta  = (*patMuons)[1].eta();
  m_Muon2M    = (*patMuons)[1].mass();
  m_Muon2MinDR  = 10.;
  m_Muon2DRJet1 = 10.;
  
  m_DiMuonInvM = ((*patMuons)[0].p4()+(*patMuons)[1].p4()).mass();
  m_DiMuonPt   = ((*patMuons)[0].p4()+(*patMuons)[1].p4()).pt();
  m_DiMuonEta  = ((*patMuons)[0].p4()+(*patMuons)[1].p4()).eta();
  m_DiMuonMinDR  = 10.;
  m_DiMuonDRJet1 = 10.;

  m_nJetsPt30Eta50 = jets  ->size();
  m_nJetsPt50Eta25 = htJets->size();
  m_nJetsCSVM = bJets ->size();
  m_nJetsCSVT = 0;
  m_HT  = *ht;
  m_MHT = (*mht)[0].pt();
  m_MET = (*met)[0].pt();
  m_ra2_HT  = *ra2HT;
  m_ra2_MHT = (*ra2MHT)[0].pt();
  m_ra2_MET = (*ra2MET)[0].pt();

  m_passLooseTopTagger = true;
  m_passLooseTopJetIdx = true;
  m_passLooseTopMassCut = true;
  m_passLooseCSVCut = true;
  m_passLooseRemainingSystem = true; 
  m_passLooseMT2Cuts = true;

  m_passNominalTopTagger = true;
  m_passNominalTopJetIdx = true;
  m_passNominalTopMassCut = true;
  m_passNominalCSVCut = true;
  m_passNominalRemainingSystem = true; 
  m_passNominalMT2Cuts = true;
  
  if (runTopTagger_) {
    m_passLooseTopTagger = true;
    m_passLooseTopJetIdx = true;
    m_passLooseTopMassCut = true;
    m_passLooseCSVCut = true;
    m_passLooseRemainingSystem = true; 
    m_passLooseMT2Cuts = true;
    
    m_passNominalTopTagger = true;
    m_passNominalTopJetIdx = true;
    m_passNominalTopMassCut = true;
    m_passNominalCSVCut = true;
    m_passNominalRemainingSystem = true; 
    m_passNominalMT2Cuts = true;

    m_loose_bestTopJetMass = *hLoosebestTopJetMass;
    m_loose_MTbJet         = *hLooseMTbJet;
    m_loose_MTbestTopJet   = *hLooseMTbestTopJet;
    m_loose_MT2            = *hLooseMT2;
    m_loose_bestTopJetIdx                    = *hLoosebestTopJetIdx;
    m_loose_remainPassCSVS                   = *hLooseremainPassCSVS;
    m_loose_pickedRemainingCombfatJetIdx     = *hLoosepickedRemainingCombfatJetIdx;
    m_loose_MTbestWJet                       = *hLooseMTbestWJet;
    m_loose_MTbestbJet                       = *hLooseMTbestbJet;
    m_loose_MTremainingTopJet                = *hLooseMTremainingTopJet;
    m_loose_linearCombMTbJetPlusMTbestTopJet = *hLooselinearCombMTbJetPlusMTbestTopJet;
    
    if( m_loose_bestTopJetIdx == -1 ) {
      m_passLooseTopJetIdx = false;
      m_passLooseTopTagger = false;
    }
    if( !(m_loose_bestTopJetMass > 80 && m_loose_bestTopJetMass < 270) ) {
      m_passLooseTopMassCut = false;
      m_passLooseTopTagger = false;
    }
    if( !m_loose_remainPassCSVS ) {
      m_passLooseCSVCut = false;
      m_passLooseTopTagger = false;
    }
    if( m_loose_pickedRemainingCombfatJetIdx == -1 && m_nJetsPt30Eta50>=6 ) {
      m_passLooseRemainingSystem = false; 
      m_passLooseTopTagger = false; 
    }
    if( !(m_loose_MT2 > 300 && (m_loose_MTbJet + 0.5*m_loose_MTbestTopJet) > 500) ) {
      m_passLooseMT2Cuts = false;
      m_passLooseTopTagger = false;
    }
    
    m_nominal_bestTopJetMass = *hNominalbestTopJetMass;
    m_nominal_MTbJet         = *hNominalMTbJet;
    m_nominal_MTbestTopJet   = *hNominalMTbestTopJet;
    m_nominal_MT2            = *hNominalMT2;
    m_nominal_bestTopJetIdx                    = *hNominalbestTopJetIdx;
    m_nominal_remainPassCSVS                   = *hNominalremainPassCSVS;
    m_nominal_pickedRemainingCombfatJetIdx     = *hNominalpickedRemainingCombfatJetIdx;
    m_nominal_MTbestWJet                       = *hNominalMTbestWJet;
    m_nominal_MTbestbJet                       = *hNominalMTbestbJet;
    m_nominal_MTremainingTopJet                = *hNominalMTremainingTopJet;
    m_nominal_linearCombMTbJetPlusMTbestTopJet = *hNominallinearCombMTbJetPlusMTbestTopJet;

    if( m_nominal_bestTopJetIdx == -1 ) {
      m_passNominalTopJetIdx = false;
      m_passNominalTopTagger = false;
    }
    if( !(m_nominal_bestTopJetMass > 80 && m_nominal_bestTopJetMass < 270) ) {
      m_passNominalTopMassCut = false;
      m_passNominalTopTagger = false;
    }
    if( !m_nominal_remainPassCSVS ) {
      m_passNominalCSVCut = false;
      m_passNominalTopTagger = false;
    }
    if( m_nominal_pickedRemainingCombfatJetIdx == -1 && m_nJetsPt30Eta50>=6 ) {
      m_passNominalRemainingSystem = false; 
      m_passNominalTopTagger = false; 
    }
    if( !(m_nominal_MT2 > 300 && (m_nominal_MTbJet + 0.5*m_nominal_MTbestTopJet) > 500) ) {
      m_passNominalMT2Cuts = false;
      m_passNominalTopTagger = false;
    }

  }
  
  const pat::Jet  *r1, *r2, *r3, *r4;
  m_dPhiMHT1 = 10.0;  m_dPhiMET1 = 10.0;
  m_dPhiMHT2 = 10.0;  m_dPhiMET2 = 10.0;
  m_dPhiMHT3 = 10.0;  m_dPhiMET3 = 10.0;
  m_dPhiMHT4 = 10.0;  m_dPhiMET4 = 10.0;
  m_Jet1Pt  = -10.;   m_Jet3Pt  = -10.;
  m_Jet1Eta = -10.;   m_Jet3Eta = -10.;
  m_Jet2Pt  = -10.;   m_Jet4Pt  = -10.;
  m_Jet2Eta = -10.;   m_Jet4Eta = -10.;

  if(m_nJetsPt30Eta50 >= 1) {
    r1 = &((*jets)[0]);
    m_dPhiMHT1 = fabs(reco::deltaPhi(r1->phi(),(*mht)[0].phi()));
    m_dPhiMET1 = fabs(reco::deltaPhi(r1->phi(),(*met)[0].phi()));
    m_Jet1Pt = r1->pt();
    m_Jet1Eta = r1->eta();
    if(m_nJetsPt30Eta50 >= 2) {
      r2 = &((*jets)[1]);
      m_dPhiMHT2 = fabs(reco::deltaPhi(r2->phi(),(*mht)[0].phi())); 
      m_dPhiMET2 = fabs(reco::deltaPhi(r2->phi(),(*met)[0].phi())); 
      m_Jet2Pt = r2->pt();
      m_Jet2Eta = r2->eta();
     
      if(m_nJetsPt30Eta50 >= 3) {
	r3 = &((*jets)[2]);
	m_dPhiMHT3 = fabs(reco::deltaPhi(r3->phi(),(*mht)[0].phi()));
	m_dPhiMET3 = fabs(reco::deltaPhi(r3->phi(),(*met)[0].phi()));
	m_Jet3Pt = r3->pt();
	m_Jet3Eta = r3->eta();

	if(m_nJetsPt30Eta50 >= 4) {
	  r4 = &((*jets)[3]);
	  m_dPhiMHT4 = fabs(reco::deltaPhi(r4->phi(),(*mht)[0].phi()));
	  m_dPhiMET4 = fabs(reco::deltaPhi(r4->phi(),(*met)[0].phi()));
	  m_Jet4Pt = r4->pt();
	  m_Jet4Eta = r4->eta();
	}
      }
    }
  }

  m_dPhiMHTMin  = 10.;
  m_dPhiMETMin  = 10.;
  m_nJetsPt30Eta24 = 0;
  m_nJetsPt50Eta24 = 0;
  m_nJetsPt70Eta24 = 0;
  m_nJetsPt50Eta25MInv = 0;
  m_HTMInv = 0.;

  if (jets->size() > 0) {
    double dR = reco::deltaR((*patMuons)[0].eta(),(*patMuons)[0].phi(),(*jets)[0].eta(), (*jets)[0].phi());
    if (m_Muon1DRJet1 > dR)
      m_Muon1DRJet1 = dR;
    dR = reco::deltaR((*patMuons)[1].eta(),(*patMuons)[1].phi(),(*jets)[0].eta(), (*jets)[0].phi());
    if (m_Muon2DRJet1 > dR)
      m_Muon2DRJet1 = dR;
    dR = reco::deltaR(((*patMuons)[0].p4()+(*patMuons)[1].p4()).eta(),
		      ((*patMuons)[0].p4()+(*patMuons)[1].p4()).phi(),
		      (*jets)[0].eta(), (*jets)[0].phi());
    if (m_DiMuonDRJet1 > dR)
      m_DiMuonDRJet1 = dR;
  }

  edm::View<pat::Jet>::const_iterator jet = jets->begin();
  for (; jet!= jets->end(); ++jet) {
    if (fabs(jet->eta() < 2.4)) {
      if (jet->pt() > 30)
	++m_nJetsPt30Eta24;
      if (jet->pt() > 50)
	++m_nJetsPt50Eta24;
      if (jet->pt() > 70)
	++m_nJetsPt70Eta24;
    }
    if (jet->pt() > 50 && fabs(jet->eta() < 2.5)) 
      if (patMuons->size()>1)
	if ((jet->p4()+((*patMuons)[0].p4()+(*patMuons)[1].p4())).mass() > 90.0) {
	  ++m_nJetsPt50Eta25MInv;
	  m_HTMInv += jet->pt();
	}
    
    double dR = reco::deltaR((*patMuons)[0].eta(),(*patMuons)[0].phi(),jet->eta(), jet->phi());
    if (m_Muon1MinDR > dR)
      m_Muon1MinDR = dR;

    dR = reco::deltaR((*patMuons)[1].eta(),(*patMuons)[1].phi(),jet->eta(), jet->phi());
    if (m_Muon2MinDR > dR)
      m_Muon2MinDR = dR;
    
    dR = reco::deltaR(((*patMuons)[0].p4()+(*patMuons)[1].p4()).eta(),
		      ((*patMuons)[0].p4()+(*patMuons)[1].p4()).phi(),
		      jet->eta(), jet->phi());
    if (m_DiMuonMinDR > dR)
      m_DiMuonMinDR = dR;
    
    double tmpDPhi = fabs(reco::deltaPhi(jet->phi(),(*mht)[0].phi()));
    if (tmpDPhi < m_dPhiMHTMin)
      m_dPhiMHTMin = tmpDPhi;

    tmpDPhi = fabs(reco::deltaPhi(jet->phi(),(*met)[0].phi()));
    if (tmpDPhi < m_dPhiMETMin)
      m_dPhiMETMin = tmpDPhi;
  }

  m_dPhiMHTMinBCSVM = 10.;
  m_dPhiMETMinBCSVM = 10.;
  m_dPhiMHTMinBCSVT = 10.;
  m_dPhiMETMinBCSVT = 10.;
  const pat::Jet  *b1, *b2;
  m_JetCSVM1Pt  = -10.;   m_JetCSVT1Pt  = -10.;
  m_JetCSVM1Eta = -10.;   m_JetCSVT1Eta = -10.;
  m_JetCSVM2Pt  = -10.;   m_JetCSVT2Pt  = -10.;
  m_JetCSVM2Eta = -10.;   m_JetCSVT2Eta = -10.;
  edm::View<pat::Jet>::const_iterator bjet = bJets->begin();
  for (; bjet!= bJets->end(); ++bjet) {
    if (bJets->size()>0){
      b1 = &((*bJets)[0]);
      m_JetCSVM1Pt  = b1->pt();
      m_JetCSVM1Eta = b1->eta();
      if (b1->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898) {
	m_JetCSVT1Pt  = b1->pt();
	m_JetCSVT1Eta = b1->eta();
      }
      if (bJets->size()>1){
	b2 = &((*bJets)[1]);
	m_JetCSVM2Pt  = b2->pt();
	m_JetCSVM2Eta = b2->eta();
	if (b2->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898) {
	  m_JetCSVT2Pt  = b2->pt();
	  m_JetCSVT2Eta = b2->eta();
	}
      }
    }
    double tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),(*mht)[0].phi()));
    if (tmpDPhiB < m_dPhiMHTMinBCSVM)
      m_dPhiMHTMinBCSVM = tmpDPhiB;

    tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),(*met)[0].phi()));
    if (tmpDPhiB < m_dPhiMETMinBCSVM)
      m_dPhiMETMinBCSVM = tmpDPhiB;
    if (bjet->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898) {
      ++m_nJetsCSVT;
      tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),(*mht)[0].phi()));
      if (tmpDPhiB < m_dPhiMHTMinBCSVT)
	m_dPhiMHTMinBCSVT = tmpDPhiB;
      
      tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),(*met)[0].phi()));
      if (tmpDPhiB < m_dPhiMETMinBCSVT)
	m_dPhiMETMinBCSVT = tmpDPhiB;
    }
  }

  /******************************************************************
   * Here we do all the HLT related trigger stuff
   *
   *
   ******************************************************************/

  /////Trigger information
  m_Mu13_Mu8 = true;
  m_Mu17_Mu8 = true;
  m_DoubleMu5_IsoMu5 = true;

  // Get the HLT results and check validity for data
  
  if (data_) {
    m_Mu13_Mu8 = false;
    m_Mu17_Mu8 = false;
    m_DoubleMu5_IsoMu5 = false;
    if (!getHLTfromConfig_) 
      if (processName_=="") {
	Handle<trigger::TriggerEvent> hltEventHandle;
	ev.getByLabel(triggerResults_, hltEventHandle);
	processName_ = hltEventHandle.provenance()->processName();
	if (debug_)
	  std::cout<<processName_<<std::endl;
      }
    hlTriggerResults_ = InputTag("TriggerResults","",processName_);
    
    edm::LogInfo("HLTEventSelector") << "Using trigger results for InputTag " << hlTriggerResults_;
    if (debug_)
      std::cout<<"Using trigger results for InputTag " << hlTriggerResults_<<std::endl;
    
    edm::Handle<edm::TriggerResults> hltHandle;
    ev.getByLabel(hlTriggerResults_, hltHandle);
    
    if ( !hltHandle.isValid() ) {
      edm::LogWarning("HLTEventSelector") << "No trigger results for InputTag " << hlTriggerResults_;
      if (debug_)
	std::cout<<"HLT results not valid"<<std::endl;
      return;
    }
    
    const edm::TriggerNames& trgNames = ev.triggerNames(*hltHandle);
    
    int          prescaleSet = hltConfig.prescaleSet(ev,es);
    
    if (debug_)
      std::cout<<"Prescale set is: "<<prescaleSet<<std::endl;
    for (unsigned int hltnum = 0; hltnum < trgNames.size(); ++hltnum) {
      std::string  tmpName     = trgNames.triggerName(hltnum);
      unsigned int trgIndex    = trgNames.triggerIndex(tmpName);
      int          trgResult   = hltHandle->accept(trgIndex);
      
      if (trgResult > 0) {
	if (tmpName.rfind("HLT_Mu13_Mu8_v") != std::string::npos)
	  m_Mu13_Mu8 = true;
	else if (tmpName.rfind("HLT_Mu17_Mu8_v") != std::string::npos)
	  m_Mu17_Mu8 = true;
	else if (tmpName.rfind("HLT_DoubleMu5_IsoMu5_v") != std::string::npos)
	  m_DoubleMu5_IsoMu5 = true;
      }
    }
  }
  //if (reducedValues)
  reducedValues->Fill();
 
}

void RA2ZInvDiMuonTreeMaker::beginJob() {
  //book trees
  BookTree();
}

void RA2ZInvDiMuonTreeMaker::endJob() {

}


void RA2ZInvDiMuonTreeMaker::BookTree() {

  edm::Service<TFileService> fs;
  reducedValues = fs->make<TTree>( "RA2Values", "Variables for reduced studies" );

  reducedValues->Branch("ra2_HT",       &m_HT,       "m_HT/D" );
  reducedValues->Branch("ra2_HTMInv",   &m_HTMInv,   "m_HTMInv/D" );
  reducedValues->Branch("ra2_MHT",      &m_MHT,      "m_MHT/D");
  reducedValues->Branch("ra2_MET",      &m_MET,      "m_MET/D");
  reducedValues->Branch("ra2_ra2HT",    &m_ra2_HT,   "m_ra2_HT/D" );
  reducedValues->Branch("ra2_ra2MHT",   &m_ra2_MHT,  "m_ra2_MHT/D");
  reducedValues->Branch("ra2_ra2MET",   &m_ra2_MET,  "m_ra2_MET/D");
  reducedValues->Branch("ra2_Vertices", &m_Vertices, "m_Vertices/I");
  reducedValues->Branch("ra2_Event",    &m_event,    "m_event/I");
  reducedValues->Branch("ra2_Run",      &m_run,      "m_run/I");
  reducedValues->Branch("ra2_Lumi",     &m_lumi,     "m_lumi/I");

  if (runTopTagger_) {
    reducedValues->Branch("ra2_loose_bestTopJetIdx",               &m_loose_bestTopJetIdx,               "m_loose_bestTopJetIdx/I");
    reducedValues->Branch("ra2_loose_pickedRemainingCombfatJetIdx",&m_loose_pickedRemainingCombfatJetIdx,"m_loose_pickedRemainingCombfatJetIdx/I");

    reducedValues->Branch("ra2_loose_remainPassCSVS",    &m_loose_remainPassCSVS,    "m_loose_remainPassCSVS/O");
    reducedValues->Branch("ra2_passLooseTopTagger",      &m_passLooseTopTagger,      "m_passLooseTopTagger/O");
    reducedValues->Branch("ra2_passLooseTopJetIdx",      &m_passLooseTopJetIdx,      "m_passLooseTopJetIdx/O");
    reducedValues->Branch("ra2_passLooseTopMassCut",     &m_passLooseTopMassCut,     "m_passLooseTopMassCut/O");
    reducedValues->Branch("ra2_passLooseCSVCut",         &m_passLooseCSVCut,         "m_passLooseCSVCut/O");
    reducedValues->Branch("ra2_passLooseRemainingSystem",&m_passLooseRemainingSystem,"m_passLooseRemainingSystem/O");
    reducedValues->Branch("ra2_passLooseMT2Cuts",        &m_passLooseMT2Cuts,        "m_passLooseMT2Cuts/O");

    reducedValues->Branch("ra2_loose_bestTopJetMass", &m_loose_bestTopJetMass, "m_loose_bestTopJetMass/D");
    reducedValues->Branch("ra2_loose_MTbJet",         &m_loose_MTbJet,         "m_loose_MTbJet/D");
    reducedValues->Branch("ra2_loose_MTbestTopJet",   &m_loose_MTbestTopJet,   "m_loose_MTbestTopJet/D");
    reducedValues->Branch("ra2_loose_MT2",            &m_loose_MT2,            "m_loose_MT2/D");
    reducedValues->Branch("ra2_loose_MTbestWJet",     &m_loose_MTbestWJet,     "m_loose_MTbestWJet/D");
    reducedValues->Branch("ra2_loose_MTbestbJet",     &m_loose_MTbestbJet,     "m_loose_MTbestbJet/D");
    reducedValues->Branch("ra2_loose_MTremainingTopJet",               &m_loose_MTremainingTopJet,               "m_loose_MTremainingTopJet/D");
    reducedValues->Branch("ra2_loose_linearCombMTbJetPlusMTbestTopJet",&m_loose_linearCombMTbJetPlusMTbestTopJet,"m_loose_linearCombMTbJetPlusMTbestTopJet/D");
    
    //////
    reducedValues->Branch("ra2_nominal_bestTopJetIdx",               &m_nominal_bestTopJetIdx,               "m_nominal_bestTopJetIdx/I");
    reducedValues->Branch("ra2_nominal_pickedRemainingCombfatJetIdx",&m_nominal_pickedRemainingCombfatJetIdx,"m_nominal_pickedRemainingCombfatJetIdx/I");

    reducedValues->Branch("ra2_nominal_remainPassCSVS",    &m_nominal_remainPassCSVS,    "m_nominal_remainPassCSVS/O");
    reducedValues->Branch("ra2_passNominalTopTagger",      &m_passNominalTopTagger,      "m_passNominalTopTagger/O");
    reducedValues->Branch("ra2_passNominalTopJetIdx",      &m_passNominalTopJetIdx,      "m_passNominalTopJetIdx/O");
    reducedValues->Branch("ra2_passNominalTopMassCut",     &m_passNominalTopMassCut,     "m_passNominalTopMassCut/O");
    reducedValues->Branch("ra2_passNominalCSVCut",         &m_passNominalCSVCut,         "m_passNominalCSVCut/O");
    reducedValues->Branch("ra2_passNominalRemainingSystem",&m_passNominalRemainingSystem,"m_passNominalRemainingSystem/O");
    reducedValues->Branch("ra2_passNominalMT2Cuts",        &m_passNominalMT2Cuts,        "m_passNominalMT2Cuts/O");

    reducedValues->Branch("ra2_nominal_bestTopJetMass", &m_nominal_bestTopJetMass, "m_nominal_bestTopJetMass/D");
    reducedValues->Branch("ra2_nominal_MTbJet",         &m_nominal_MTbJet,         "m_nominal_MTbJet/D");
    reducedValues->Branch("ra2_nominal_MTbestTopJet",   &m_nominal_MTbestTopJet,   "m_nominal_MTbestTopJet/D");
    reducedValues->Branch("ra2_nominal_MT2",            &m_nominal_MT2,            "m_nominal_MT2/D");
    reducedValues->Branch("ra2_nominal_MTbestWJet",     &m_nominal_MTbestWJet,     "m_nominal_MTbestWJet/D");
    reducedValues->Branch("ra2_nominal_MTbestbJet",     &m_nominal_MTbestbJet,     "m_nominal_MTbestbJet/D");
    reducedValues->Branch("ra2_nominal_MTremainingTopJet",               &m_nominal_MTremainingTopJet,               "m_nominal_MTremainingTopJet/D");
    reducedValues->Branch("ra2_nominal_linearCombMTbJetPlusMTbestTopJet",&m_nominal_linearCombMTbJetPlusMTbestTopJet,"m_nominal_linearCombMTbJetPlusMTbestTopJet/D");
  }

  reducedValues->Branch("ra2_dPhiMHT1", &m_dPhiMHT1, "m_dPhiMHT1/D");
  reducedValues->Branch("ra2_dPhiMHT2", &m_dPhiMHT2, "m_dPhiMHT2/D");
  reducedValues->Branch("ra2_dPhiMHT3", &m_dPhiMHT3, "m_dPhiMHT3/D");
  reducedValues->Branch("ra2_dPhiMHT4", &m_dPhiMHT4, "m_dPhiMHT4/D");
  reducedValues->Branch("ra2_dPhiMHTMin", &m_dPhiMHTMin, "m_dPhiMHTMin/D");
  reducedValues->Branch("ra2_dPhiMHTMinBCSVM", &m_dPhiMHTMinBCSVM, "m_dPhiMHTMinBCSVM/D");
  reducedValues->Branch("ra2_dPhiMHTMinBCSVT", &m_dPhiMHTMinBCSVT, "m_dPhiMHTMinBCSVT/D");

  reducedValues->Branch("ra2_dPhiMET1", &m_dPhiMET1, "m_dPhiMET1/D");
  reducedValues->Branch("ra2_dPhiMET2", &m_dPhiMET2, "m_dPhiMET2/D");
  reducedValues->Branch("ra2_dPhiMET3", &m_dPhiMET3, "m_dPhiMET3/D");
  reducedValues->Branch("ra2_dPhiMET4", &m_dPhiMET4, "m_dPhiMET4/D");
  reducedValues->Branch("ra2_dPhiMETMin", &m_dPhiMETMin, "m_dPhiMETMin/D");
  reducedValues->Branch("ra2_dPhiMETMinBCSVM", &m_dPhiMETMinBCSVM, "m_dPhiMETMinBCSVM/D");
  reducedValues->Branch("ra2_dPhiMETMinBCSVT", &m_dPhiMETMinBCSVT, "m_dPhiMETMinBCSVT/D");

  reducedValues->Branch("ra2_Jet1Pt",  &m_Jet1Pt,  "m_Jet1Pt/D");
  reducedValues->Branch("ra2_Jet1Eta", &m_Jet1Eta, "m_Jet1Eta/D");
  reducedValues->Branch("ra2_Jet2Pt",  &m_Jet2Pt,  "m_Jet2Pt/D");
  reducedValues->Branch("ra2_Jet2Eta", &m_Jet2Eta, "m_Jet2Eta/D");
  reducedValues->Branch("ra2_Jet3Pt",  &m_Jet3Pt,  "m_Jet3Pt/D");
  reducedValues->Branch("ra2_Jet3Eta", &m_Jet3Eta, "m_Jet3Eta/D");
  reducedValues->Branch("ra2_Jet4Pt",  &m_Jet4Pt,  "m_Jet4Pt/D");
  reducedValues->Branch("ra2_Jet4Eta", &m_Jet4Eta, "m_Jet4Eta/D");

  reducedValues->Branch("ra2_JetCSVM1Pt",  &m_JetCSVM1Pt,  "m_JetCSVM1Pt/D");
  reducedValues->Branch("ra2_JetCSVM1Eta", &m_JetCSVM1Eta, "m_JetCSVM1Eta/D");
  reducedValues->Branch("ra2_JetCSVM2Pt",  &m_JetCSVM2Pt,  "m_JetCSVM2Pt/D");
  reducedValues->Branch("ra2_JetCSVM2Eta", &m_JetCSVM2Eta, "m_JetCSVM2Eta/D");
  reducedValues->Branch("ra2_JetCSVT1Pt",  &m_JetCSVT1Pt,  "m_JetCSVT1Pt/D");
  reducedValues->Branch("ra2_JetCSVT1Eta", &m_JetCSVT1Eta, "m_JetCSVT1Eta/D");
  reducedValues->Branch("ra2_JetCSVT2Pt",  &m_JetCSVT2Pt,  "m_JetCSVT2Pt/D");
  reducedValues->Branch("ra2_JetCSVT2Eta", &m_JetCSVT2Eta, "m_JetCSVT2Eta/D");

  reducedValues->Branch("ra2_PUWt",    &m_PUWt,    "m_PUWt/D");
  reducedValues->Branch("ra2_EventWt", &m_EventWt, "m_EventWt/D");

  reducedValues->Branch("ra2_nMuonsIso",   &m_nMuonsIso,   "m_nMuonsIso/I");
  reducedValues->Branch("ra2_Muon1Pt",     &m_Muon1Pt,     "m_Muon1Pt/D" );
  reducedValues->Branch("ra2_Muon1Eta",    &m_Muon1Eta,    "m_Muon1Eta/D");
  reducedValues->Branch("ra2_Muon1M",      &m_Muon1M,      "m_Muon1M/D");
  reducedValues->Branch("ra2_Muon1MinDR",  &m_Muon1MinDR,  "m_Muon1MinDR/D" );
  reducedValues->Branch("ra2_Muon1DRJet1", &m_Muon1DRJet1, "m_Muon1DRJet1/D" );
  reducedValues->Branch("ra2_Muon2Pt",     &m_Muon2Pt,     "m_Muon2Pt/D" );
  reducedValues->Branch("ra2_Muon2Eta",    &m_Muon2Eta,    "m_Muon2Eta/D");
  reducedValues->Branch("ra2_Muon2M",      &m_Muon2M,      "m_Muon2M/D");
  reducedValues->Branch("ra2_Muon2MinDR",  &m_Muon2MinDR,  "m_Muon2MinDR/D" );
  reducedValues->Branch("ra2_Muon2DRJet1", &m_Muon2DRJet1, "m_Muon2DRJet1/D" );
  reducedValues->Branch("ra2_DiMuonInvM",  &m_DiMuonInvM,  "m_DiMuonInvM/D");
  reducedValues->Branch("ra2_DiMuonPt",    &m_DiMuonPt,    "m_DiMuonPt/D");
  reducedValues->Branch("ra2_DiMuonEta",   &m_DiMuonEta,   "m_DiMuonEta/D");
  reducedValues->Branch("ra2_DiMuonMinDR", &m_DiMuonMinDR, "m_DiMuonMinDR/D" );
  reducedValues->Branch("ra2_DiMuonDRJet1",&m_DiMuonDRJet1,"m_DiMuonDRJet1/D" );

  reducedValues->Branch("ra2_Mu13_Mu8",         &m_Mu13_Mu8        , "m_Mu13_Mu8/O"        );
  reducedValues->Branch("ra2_Mu17_Mu8",         &m_Mu17_Mu8        , "m_Mu17_Mu8/O"        );
  reducedValues->Branch("ra2_DoubleMu5_IsoMu5", &m_DoubleMu5_IsoMu5, "m_DoubleMu5_IsoMu5/O");

  reducedValues->Branch("ra2_nJetsCSVM", &m_nJetsCSVM, "m_nJetsCSVM/I");
  reducedValues->Branch("ra2_nJetsCSVT", &m_nJetsCSVT, "m_nJetsCSVT/I");
  reducedValues->Branch("ra2_nJetsPt30Eta50", &m_nJetsPt30Eta50, "m_nJetsPt30Eta50/I" );
  reducedValues->Branch("ra2_nJetsPt30Eta24", &m_nJetsPt30Eta24, "m_nJetsPt30Eta24/I");
  reducedValues->Branch("ra2_nJetsPt50Eta24", &m_nJetsPt50Eta24, "m_nJetsPt50Eta24/I");
  reducedValues->Branch("ra2_nJetsPt70Eta24", &m_nJetsPt70Eta24, "m_nJetsPt70Eta24/I");
  reducedValues->Branch("ra2_nJetsPt50Eta25", &m_nJetsPt50Eta25, "m_nJetsPt50Eta25/I" );
  reducedValues->Branch("ra2_nJetsPt50Eta25MInv", &m_nJetsPt50Eta25MInv, "m_nJetsPt50Eta25MInv/I" );

  if (storeExtraVetos_) {
    reducedValues->Branch("ra2_passRA2ElVeto",    &m_passRA2ElVeto   , "m_passRA2ElVeto/O"    );
    reducedValues->Branch("ra2_passRA2MuVeto",    &m_passRA2MuVeto   , "m_passRA2MuVeto/O"    );
    reducedValues->Branch("ra2_passDirIsoElVeto", &m_passDirIsoElVeto, "m_passDirIsoElVeto/O"    );
    reducedValues->Branch("ra2_passDirIsoMuVeto", &m_passDirIsoMuVeto, "m_passDirIsoMuVeto/O"    );
    reducedValues->Branch("ra2_passIsoTrkVeto",   &m_passIsoTrkVeto,   "m_passIsoTrkVeto/O");
  }

  reducedValues->SetAutoSave(1);
}


void  RA2ZInvDiMuonTreeMaker::beginRun(edm::Run const& run, edm::EventSetup const& es) {
  bool changed = false;
  if (data_) {
    if (hltConfig.init(run,es,"HLT",changed)) {
      if (changed) {
	edm::LogWarning("RA2ZInvDiMuonTreeMaker") << "beginRun: The HLT config has changed!";
      }
    }
    else {
      edm::LogError("TriggerEvent") << " HLT config extraction failure";
    }
  }
}

void RA2ZInvDiMuonTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void RA2ZInvDiMuonTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvDiMuonTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvDiMuonTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RA2ZInvDiMuonTreeMaker);

