// -*- C++ -*-
//
// Package:    RA2ZInvTreeMaker
// Class:      RA2ZInvTreeMaker
// 
/**\class RA2ZInvTreeMaker RA2ZInvTreeMaker.cc SusyAnalysis/RA2ZInvTreeMaker/src/RA2ZInvTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seema Sharma
//         Created:  Mon Jun 20 12:58:08 CDT 2011
// $Id: RA2ZInvTreeMaker.cc,v 1.9 2013/02/25 13:51:38 sturdy Exp $
//
//


// system include files
#include <memory>
#include <iomanip>
#include <cmath>

#include "ZInvisibleBkgds/Photons/interface/RA2ZInvTreeMaker.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include <DataFormats/METReco/interface/MET.h>


RA2ZInvTreeMaker::RA2ZInvTreeMaker(const edm::ParameterSet& pset) {

  // read parameters from config file
  debug_          = pset.getParameter<bool>("Debug");
  scale_          = pset.getParameter<double>("ScaleFactor");

  //  genSrc_         = pset.getParameter<edm::InputTag>("genSrc" );
  genJetSrc_      = pset.getParameter<edm::InputTag>("genJetSrc" );
  genMETSrc_      = pset.getParameter<edm::InputTag>("genMETSrc" );
  genLabel_       = pset.getParameter<edm::InputTag>("genLabel");  

  vertexSrc_      = pset.getParameter<edm::InputTag>("VertexSrc");
  jetSrc_         = pset.getParameter<edm::InputTag>("JetSrc");
  bJetSrc_        = pset.getParameter<edm::InputTag>("bJetSrc");
  htJetSrc_       = pset.getParameter<edm::InputTag>("htJetSrc");
  htSrc_          = pset.getParameter<edm::InputTag>("htSource");
  mhtSrc_         = pset.getParameter<edm::InputTag>("mhtSource");
  metSrc_         = pset.getParameter<edm::InputTag>("metSource");

//  runTopTagger_        = pset.getParameter<bool>("runTopTagger");
//  looseTopTaggerSrc_   = pset.getParameter<std::string>("looseTopTaggerSource");
//  nominalTopTaggerSrc_ = pset.getParameter<std::string>("nominalTopTaggerSource");

  doPUReWeight_    = pset.getParameter<bool>("DoPUReweight");
  storeExtraVetos_ = pset.getParameter<bool>("storeExtraVetos");

  ra2ElectronSrc_  = pset.getParameter<edm::InputTag>("ra2ElectronForVeto");
  ra2MuonSrc_      = pset.getParameter<edm::InputTag>("ra2MuonForVeto");
  electronVetoSrc_ = pset.getParameter<edm::InputTag>("electronVetoSource");
  muonVetoSrc_     = pset.getParameter<edm::InputTag>("muonVetoSource");
  isoTrkVetoSrc_   = pset.getParameter<edm::InputTag>("isoTrkVetoSource");
  puWeightSrc_     = pset.getParameter<edm::InputTag>("PUWeightSource");
  eventWeightSrc_  = pset.getParameter<edm::InputTag>("EventWeightSource");

}

RA2ZInvTreeMaker::~RA2ZInvTreeMaker() {
  delete reducedValues;
  reducedValues = 0;
}

void RA2ZInvTreeMaker::analyze(const edm::Event& ev, const edm::EventSetup& es) {

  using namespace edm;

  if (ev.isRealData()) {
    std::cout<<"Trying to run Zvv TreeMaker on data is not going to end well for you..."<<std::endl;
    return;
  }
  // get event-id
  unsigned int event = (ev.id()).event();
  unsigned int run   = (ev.id()).run();
  unsigned int lumi  =  ev.luminosityBlock();

  m_event = event;
  m_run   = run;
  m_lumi  = lumi;

  //gen level information
  edm::Handle<reco::GenParticleCollection> gens;
  ev.getByLabel(genLabel_,gens);
  //edm::Handle<reco::GenJetCollection> genJets;
  edm::Handle<edm::View<reco::GenJet> > genJets;
  ev.getByLabel(genJetSrc_,genJets);
  edm::Handle<edm::View<reco::GenMET> > genMET;
  ev.getByLabel(genMETSrc_,genMET);

  //Reco information
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

  //if (runTopTagger_) {
  //  ev.getByLabel(looseTopTaggerSrc_,"bestTopJetMass", hLoosebestTopJetMass);
  //  ev.getByLabel(looseTopTaggerSrc_,"mTbJet"        , hLooseMTbJet);
  //  ev.getByLabel(looseTopTaggerSrc_,"mTbestTopJet"  , hLooseMTbestTopJet);
  //  ev.getByLabel(looseTopTaggerSrc_,"MT2"           , hLooseMT2);
  //  ev.getByLabel(looseTopTaggerSrc_,"bestTopJetIdx"                   , hLoosebestTopJetIdx);
  //  ev.getByLabel(looseTopTaggerSrc_,"remainPassCSVS"                  , hLooseremainPassCSVS);
  //  ev.getByLabel(looseTopTaggerSrc_,"pickedRemainingCombfatJetIdx"    , hLoosepickedRemainingCombfatJetIdx);
  //  ev.getByLabel(looseTopTaggerSrc_,"mTbestWJet"                      , hLooseMTbestWJet);
  //  ev.getByLabel(looseTopTaggerSrc_,"mTbestbJet"                      , hLooseMTbestbJet);
  //  ev.getByLabel(looseTopTaggerSrc_,"mTremainingTopJet"               , hLooseMTremainingTopJet);
  //  ev.getByLabel(looseTopTaggerSrc_,"linearCombmTbJetPlusmTbestTopJet", hLooselinearCombMTbJetPlusMTbestTopJet);
  //
  //  ev.getByLabel(nominalTopTaggerSrc_,"bestTopJetMass", hNominalbestTopJetMass);
  //  ev.getByLabel(nominalTopTaggerSrc_,"mTbJet"        , hNominalMTbJet);
  //  ev.getByLabel(nominalTopTaggerSrc_,"mTbestTopJet"  , hNominalMTbestTopJet);
  //  ev.getByLabel(nominalTopTaggerSrc_,"MT2"           , hNominalMT2);
  //  ev.getByLabel(nominalTopTaggerSrc_,"bestTopJetIdx"                   , hNominalbestTopJetIdx);
  //  ev.getByLabel(nominalTopTaggerSrc_,"remainPassCSVS"                  , hNominalremainPassCSVS);
  //  ev.getByLabel(nominalTopTaggerSrc_,"pickedRemainingCombfatJetIdx"    , hNominalpickedRemainingCombfatJetIdx);
  //  ev.getByLabel(nominalTopTaggerSrc_,"mTbestWJet"                      , hNominalMTbestWJet);
  //  ev.getByLabel(nominalTopTaggerSrc_,"mTbestbJet"                      , hNominalMTbestbJet);
  //  ev.getByLabel(nominalTopTaggerSrc_,"mTremainingTopJet"               , hNominalMTremainingTopJet);
  //  ev.getByLabel(nominalTopTaggerSrc_,"linearCombmTbJetPlusmTbestTopJet", hNominallinearCombMTbJetPlusMTbestTopJet);
  //}

  if (debug_) {
    std::cout<<"vertex collection has size "<<vertices->size()<<std::endl;
    std::cout<<"Jet collection has size "<<jets->size()<<std::endl;
    std::cout<<"HT Jet collection has size "<<htJets->size()<<std::endl;
    std::cout<<"b-Jet collection has size "<<bJets->size()<<std::endl;
    std::cout<<"HT value "<<*ht<<std::endl;
    std::cout<<"MHT value "<<(*mht)[0].pt()<<std::endl;
    std::cout<<"MET value "<<(*met)[0].pt()<<std::endl;

    //if (runTopTagger_) {
    //  std::cout<<"hLoosebestTopJetMass value "<<*hLoosebestTopJetMass<<std::endl;
    //  std::cout<<"hLooseMTbJet value "        <<*hLooseMTbJet<<std::endl;
    //  std::cout<<"hLooseMTbestTopJet value "  <<*hLooseMTbestTopJet<<std::endl;
    //  std::cout<<"hLooseMT2 value "           <<*hLooseMT2<<std::endl;
    //  std::cout<<"hLoosebestTopJetIdx"                   <<*hLoosebestTopJetIdx<<std::endl;
    //  std::cout<<"hLooseremainPassCSVS"                  <<*hLooseremainPassCSVS<<std::endl;
    //  std::cout<<"hLoosepickedRemainingCombfatJetIdx"    <<*hLoosepickedRemainingCombfatJetIdx<<std::endl;
    //  std::cout<<"hLooseMTbestWJet"                      <<*hLooseMTbestWJet<<std::endl;
    //  std::cout<<"hLooseMTbestbJet"                      <<*hLooseMTbestbJet<<std::endl;
    //  std::cout<<"hLooseMTremainingTopJet"               <<*hLooseMTremainingTopJet<<std::endl;
    //  std::cout<<"hLooselinearCombMTbJetPlusMTbestTopJet"<<*hLooselinearCombMTbJetPlusMTbestTopJet<<std::endl;
    //  
    //  std::cout<<"hNominalbestTopJetMass value "<<*hNominalbestTopJetMass<<std::endl;
    //  std::cout<<"hNominalMTbJet value "        <<*hNominalMTbJet<<std::endl;
    //  std::cout<<"hNominalMTbestTopJet value "  <<*hNominalMTbestTopJet<<std::endl;
    //  std::cout<<"hNominalMT2 value "           <<*hNominalMT2<<std::endl;
    //  std::cout<<"hNominalbestTopJetIdx"                   <<*hNominalbestTopJetIdx<<std::endl;
    //  std::cout<<"hNominalremainPassCSVS"                  <<*hNominalremainPassCSVS<<std::endl;
    //  std::cout<<"hNominalpickedRemainingCombfatJetIdx"    <<*hNominalpickedRemainingCombfatJetIdx<<std::endl;
    //  std::cout<<"hNominalMTbestWJet"                      <<*hNominalMTbestWJet<<std::endl;
    //  std::cout<<"hNominalMTbestbJet"                      <<*hNominalMTbestbJet<<std::endl;
    //  std::cout<<"hNominalMTremainingTopJet"               <<*hNominalMTremainingTopJet<<std::endl;
    //  std::cout<<"hNominallinearCombMTbJetPlusMTbestTopJet"<<*hNominallinearCombMTbJetPlusMTbestTopJet<<std::endl;
    //}
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
  m_PUWt    = pu_event_wt;
  m_Vertices = vertices->size();

  edm::Handle<bool> ra2elVeto;
  edm::Handle<bool> ra2muVeto;
  edm::Handle<bool> elVeto;
  edm::Handle<bool> muVeto;
  edm::Handle<bool> isoTrkVeto;

  //  if (storeExtraVetos_) {
  ev.getByLabel(ra2ElectronSrc_,ra2elVeto);
  m_passRA2ElVeto = *ra2elVeto;
  ev.getByLabel(ra2MuonSrc_,ra2muVeto);
  m_passRA2MuVeto = *ra2muVeto;
  ev.getByLabel(electronVetoSrc_,elVeto);
  m_passDirIsoElVeto = *elVeto;
  ev.getByLabel(muonVetoSrc_,muVeto);
  m_passDirIsoMuVeto = *muVeto;
  ev.getByLabel(isoTrkVetoSrc_,isoTrkVeto);
  m_passIsoTrkVeto = *isoTrkVeto;
  //}
  //////
  if(debug_ && gens->size() > 0) {
    std::cout << "Gen Z-bosons : " << std::endl;
    for( unsigned int igen=0; igen<gens->size(); igen++) {
      std::cout << igen << " " <<(*gens)[igen].pt() 
		<< " " << (*gens)[igen].eta()
		<< " " << (*gens)[igen].phi() 
		<< " " << (*gens)[igen].mass() 
		<< std::endl;
    }
  }
  //////
  m_genBoson1Pt  = -10.;
  m_genBoson1Eta = -10.;
  m_genBoson1M   = -10.;
  m_genBoson1MinDR  = 10.;
  m_genBoson1DRJet1 = 10.;
  m_genBosons  = gens->size();
  if (gens->size() > 0) {
    m_genBoson1Pt  = (*gens)[0].pt();
    m_genBoson1Eta = (*gens)[0].eta();
    m_genBoson1M   = (*gens)[0].mass();
    
    if (jets->size() > 0) {
      double dR = reco::deltaR((*gens)[0].eta(),(*gens)[0].phi(),(*jets)[0].eta(), (*jets)[0].phi());
      if (m_genBoson1DRJet1 > dR)
	m_genBoson1DRJet1 = dR;
    }

    edm::View<pat::Jet>::const_iterator jet = jets->begin();
    for (; jet != jets->end(); ++jet) {
      double dR = reco::deltaR((*gens)[0].eta(),(*gens)[0].phi(),jet->eta(), jet->phi());
      if (m_genBoson1MinDR > dR)
	m_genBoson1MinDR = dR;
    }
  }

  m_nGenJetsPt30Eta50 = 0;
  m_nGenJetsPt50Eta25 = 0;

  m_genMHT = 0;
  m_genHT  = 0;

  m_genDPhiMHT1 = 10.0;   m_genDPhiMET1 = 10.0;
  m_genDPhiMHT2 = 10.0;   m_genDPhiMET2 = 10.0;
  m_genDPhiMHT3 = 10.0;	  m_genDPhiMET3 = 10.0;
  m_genDPhiMHT4 = 10.0;   m_genDPhiMET4 = 10.0;
  m_genJet1Pt  = -10.;   m_genJet3Pt  = -10.;
  m_genJet1Eta = -10.;   m_genJet3Eta = -10.;
  m_genJet2Pt  = -10.;   m_genJet4Pt  = -10.;
  m_genJet2Eta = -10.;   m_genJet4Eta = -10.;
  
  m_genDPhiMHTMin  = 10.;
  m_genDPhiMETMin  = 10.;
  std::vector<double*> genPtVals;
  genPtVals.push_back(&m_genJet1Pt);  genPtVals.push_back(&m_genJet2Pt);  genPtVals.push_back(&m_genJet3Pt);  genPtVals.push_back(&m_genJet4Pt);
  std::vector<double*> genEtaVals;
  genEtaVals.push_back(&m_genJet1Eta);  genEtaVals.push_back(&m_genJet2Eta);  genEtaVals.push_back(&m_genJet3Eta);  genEtaVals.push_back(&m_genJet4Eta);
  std::vector<double*> genDPhiValsMHT;
  genDPhiValsMHT.push_back(&m_genDPhiMHT1);  genDPhiValsMHT.push_back(&m_genDPhiMHT2);  genDPhiValsMHT.push_back(&m_genDPhiMHT3);  genDPhiValsMHT.push_back(&m_genDPhiMHT4);
  std::vector<double*> genDPhiValsMET;
  genDPhiValsMET.push_back(&m_genDPhiMET1);  genDPhiValsMET.push_back(&m_genDPhiMET2);  genDPhiValsMET.push_back(&m_genDPhiMET3);  genDPhiValsMET.push_back(&m_genDPhiMET4);

  RA2ZInvTreeMakerFunctions::computeGenJetVars((*genJets),(*genMET)[0],
					       m_genHT,m_genMHT,
					       &genPtVals,&genEtaVals,
					       &genDPhiValsMHT,&genDPhiValsMET,
					       m_genDPhiMHTMin,m_genDPhiMETMin,
					       m_nGenJetsPt30Eta50,m_nGenJetsPt50Eta25
					       );
  
  ///////////
  m_nJetsPt30Eta50 = jets  ->size();
  m_nJetsPt50Eta25 = htJets->size();
  m_nJetsCSVM = bJets ->size();
  m_nJetsCSVT = 0;
  m_HT  = *ht;
  m_MHT = (*mht)[0].pt();
  m_MET = (*met)[0].pt();

  //m_passLooseTopTagger = true;
  //m_passLooseTopJetIdx = true;
  //m_passLooseTopMassCut = true;
  //m_passLooseCSVCut = true;
  //m_passLooseRemainingSystem = true; 
  //m_passLooseMT2Cuts = true;
  //
  //m_passNominalTopTagger = true;
  //m_passNominalTopJetIdx = true;
  //m_passNominalTopMassCut = true;
  //m_passNominalCSVCut = true;
  //m_passNominalRemainingSystem = true; 
  //m_passNominalMT2Cuts = true;
  //
  //if (runTopTagger_) {
  //  m_passLooseTopTagger = true;
  //  m_passLooseTopJetIdx = true;
  //  m_passLooseTopMassCut = true;
  //  m_passLooseCSVCut = true;
  //  m_passLooseRemainingSystem = true; 
  //  m_passLooseMT2Cuts = true;
  //  
  //  m_passNominalTopTagger = true;
  //  m_passNominalTopJetIdx = true;
  //  m_passNominalTopMassCut = true;
  //  m_passNominalCSVCut = true;
  //  m_passNominalRemainingSystem = true; 
  //  m_passNominalMT2Cuts = true;
  //
  //  m_loose_bestTopJetMass = *hLoosebestTopJetMass;
  //  m_loose_MTbJet         = *hLooseMTbJet;
  //  m_loose_MTbestTopJet   = *hLooseMTbestTopJet;
  //  m_loose_MT2            = *hLooseMT2;
  //  m_loose_bestTopJetIdx                    = *hLoosebestTopJetIdx;
  //  m_loose_remainPassCSVS                   = *hLooseremainPassCSVS;
  //  m_loose_pickedRemainingCombfatJetIdx     = *hLoosepickedRemainingCombfatJetIdx;
  //  m_loose_MTbestWJet                       = *hLooseMTbestWJet;
  //  m_loose_MTbestbJet                       = *hLooseMTbestbJet;
  //  m_loose_MTremainingTopJet                = *hLooseMTremainingTopJet;
  //  m_loose_linearCombMTbJetPlusMTbestTopJet = *hLooselinearCombMTbJetPlusMTbestTopJet;
  //  
  //  if( m_loose_bestTopJetIdx == -1 ) {
  //    m_passLooseTopJetIdx = false;
  //    m_passLooseTopTagger = false;
  //  }
  //  if( !(m_loose_bestTopJetMass > 80 && m_loose_bestTopJetMass < 270) ) {
  //    m_passLooseTopMassCut = false;
  //    m_passLooseTopTagger = false;
  //  }
  //  if( !m_loose_remainPassCSVS ) {
  //    m_passLooseCSVCut = false;
  //    m_passLooseTopTagger = false;
  //  }
  //  if( m_loose_pickedRemainingCombfatJetIdx == -1 && m_nJetsPt30Eta50>=6 ) {
  //    m_passLooseRemainingSystem = false; 
  //    m_passLooseTopTagger = false; 
  //  }
  //  if( !(m_loose_MT2 > 300 && (m_loose_MTbJet + 0.5*m_loose_MTbestTopJet) > 500) ) {
  //    m_passLooseMT2Cuts = false;
  //    m_passLooseTopTagger = false;
  //  }
  //  
  //  m_nominal_bestTopJetMass = *hNominalbestTopJetMass;
  //  m_nominal_MTbJet         = *hNominalMTbJet;
  //  m_nominal_MTbestTopJet   = *hNominalMTbestTopJet;
  //  m_nominal_MT2            = *hNominalMT2;
  //  m_nominal_bestTopJetIdx                    = *hNominalbestTopJetIdx;
  //  m_nominal_remainPassCSVS                   = *hNominalremainPassCSVS;
  //  m_nominal_pickedRemainingCombfatJetIdx     = *hNominalpickedRemainingCombfatJetIdx;
  //  m_nominal_MTbestWJet                       = *hNominalMTbestWJet;
  //  m_nominal_MTbestbJet                       = *hNominalMTbestbJet;
  //  m_nominal_MTremainingTopJet                = *hNominalMTremainingTopJet;
  //  m_nominal_linearCombMTbJetPlusMTbestTopJet = *hNominallinearCombMTbJetPlusMTbestTopJet;
  //
  //  if( m_nominal_bestTopJetIdx == -1 ) {
  //    m_passNominalTopJetIdx = false;
  //    m_passNominalTopTagger = false;
  //  }
  //  if( !(m_nominal_bestTopJetMass > 80 && m_nominal_bestTopJetMass < 270) ) {
  //    m_passNominalTopMassCut = false;
  //    m_passNominalTopTagger = false;
  //  }
  //  if( !m_nominal_remainPassCSVS ) {
  //    m_passNominalCSVCut = false;
  //    m_passNominalTopTagger = false;
  //  }
  //  if( m_nominal_pickedRemainingCombfatJetIdx == -1 && m_nJetsPt30Eta50>=6 ) {
  //    m_passNominalRemainingSystem = false; 
  //    m_passNominalTopTagger = false; 
  //  }
  //  if( !(m_nominal_MT2 > 300 && (m_nominal_MTbJet + 0.5*m_nominal_MTbestTopJet) > 500) ) {
  //    m_passNominalMT2Cuts = false;
  //    m_passNominalTopTagger = false;
  //  }
  //
  //}

  m_dPhiMHT1 = 10.0;  m_dPhiMET1 = 10.0;
  m_dPhiMHT2 = 10.0;  m_dPhiMET2 = 10.0;
  m_dPhiMHT3 = 10.0;  m_dPhiMET3 = 10.0;
  m_dPhiMHT4 = 10.0;  m_dPhiMET4 = 10.0;
  m_Jet1Pt  = -10.;   m_Jet3Pt  = -10.;
  m_Jet1Eta = -10.;   m_Jet3Eta = -10.;
  m_Jet2Pt  = -10.;   m_Jet4Pt  = -10.;
  m_Jet2Eta = -10.;   m_Jet4Eta = -10.;

  m_dPhiMHTMin  = 10.;
  m_dPhiMETMin  = 10.;
  m_nJetsPt30Eta24 = 0;
  m_nJetsPt50Eta24 = 0;
  m_nJetsPt70Eta24 = 0;

  std::vector<double*> ptVals;
  ptVals.push_back(&m_Jet1Pt);  ptVals.push_back(&m_Jet2Pt);  ptVals.push_back(&m_Jet3Pt);  ptVals.push_back(&m_Jet4Pt);
  std::vector<double*> etaVals;
  etaVals.push_back(&m_Jet1Eta);  etaVals.push_back(&m_Jet2Eta);  etaVals.push_back(&m_Jet3Eta);  etaVals.push_back(&m_Jet4Eta);
  std::vector<double*> dphiValsMHT;
  dphiValsMHT.push_back(&m_dPhiMHT1);  dphiValsMHT.push_back(&m_dPhiMHT2);  dphiValsMHT.push_back(&m_dPhiMHT3);  dphiValsMHT.push_back(&m_dPhiMHT4);
  std::vector<double*> dphiValsMET;
  dphiValsMET.push_back(&m_dPhiMET1);  dphiValsMET.push_back(&m_dPhiMET2);  dphiValsMET.push_back(&m_dPhiMET3);  dphiValsMET.push_back(&m_dPhiMET4);
  
  RA2ZInvTreeMakerFunctions::computeAllJetVars((*jets),(*mht)[0],(*met)[0],
					       &ptVals,&etaVals,
					       &dphiValsMHT,&dphiValsMET,
					       m_dPhiMHTMin, m_dPhiMETMin,
					       m_nJetsPt30Eta24, m_nJetsPt50Eta24, m_nJetsPt70Eta24);

  m_nJetsPt50Eta25MInv = 0;
  m_HTMInv = 0.;
  edm::View<pat::Jet>::const_iterator jet = jets->begin();
  for (; jet!= jets->end(); ++jet) {
    if (jet->pt() > 50 && fabs(jet->eta()) < 2.5) 
      if (gens->size())
	if ((jet->p4()+(*gens)[0].p4()).mass() > 90.0) {
	  ++m_nJetsPt50Eta25MInv;
	  m_HTMInv += jet->pt();
	}
    
  }

  m_dPhiMHTMinBCSVM = 10.;
  m_dPhiMETMinBCSVM = 10.;
  m_dPhiMHTMinBCSVT = 10.;
  m_dPhiMETMinBCSVT = 10.;
  m_JetCSVM1Pt  = -10.;   m_JetCSVT1Pt  = -10.;
  m_JetCSVM1Eta = -10.;   m_JetCSVT1Eta = -10.;
  m_JetCSVM2Pt  = -10.;   m_JetCSVT2Pt  = -10.;
  m_JetCSVM2Eta = -10.;   m_JetCSVT2Eta = -10.;
  std::vector<double*> ptValsCSVM;
  std::vector<double*> ptValsCSVT;
  std::vector<double*> etaValsCSVM;
  std::vector<double*> etaValsCSVT;
  ptValsCSVM.push_back(&m_JetCSVM1Pt);
  ptValsCSVM.push_back(&m_JetCSVM2Pt);
  ptValsCSVT.push_back(&m_JetCSVT1Pt);
  ptValsCSVT.push_back(&m_JetCSVT2Pt);
  etaValsCSVM.push_back(&m_JetCSVM1Eta);
  etaValsCSVM.push_back(&m_JetCSVM2Eta);
  etaValsCSVT.push_back(&m_JetCSVT1Eta);
  etaValsCSVT.push_back(&m_JetCSVT2Eta);
  
  RA2ZInvTreeMakerFunctions::computeBJetVars(*bJets,(*mht)[0],(*met)[0],
					     m_nJetsCSVT,
					     &ptValsCSVM,&etaValsCSVM,
					     &ptValsCSVT,&etaValsCSVT,
					     m_dPhiMHTMinBCSVM, m_dPhiMHTMinBCSVT,
					     m_dPhiMETMinBCSVM, m_dPhiMETMinBCSVT
					     );

  //if (reducedValues)
  reducedValues->Fill();
}


void RA2ZInvTreeMaker::beginJob() {
  //book histograms
  BookTree();
}

void RA2ZInvTreeMaker::endJob() {

}


void RA2ZInvTreeMaker::BookTree() {

  edm::Service<TFileService> fs;
  reducedValues = fs->make<TTree>( "RA2Values", "Variables for reduced studies" );

  reducedValues->Branch("ra2_HT",       &m_HT,       "ra2_HT/D" );
  reducedValues->Branch("ra2_HTMInv",   &m_HTMInv,   "ra2_HTMInv/D" );
  reducedValues->Branch("ra2_MHT",      &m_MHT,      "ra2_MHT/D");
  reducedValues->Branch("ra2_MET",      &m_MET,      "ra2_MET/D");
  reducedValues->Branch("ra2_Vertices", &m_Vertices, "ra2_Vertices/I");
  reducedValues->Branch("ra2_Event",    &m_event,    "ra2_event/I");
  reducedValues->Branch("ra2_Run",      &m_run,      "ra2_run/I");
  reducedValues->Branch("ra2_Lumi",     &m_lumi,     "ra2_lumi/I");

  //if (runTopTagger_) {
  //  reducedValues->Branch("ra2_loose_bestTopJetIdx",               &m_loose_bestTopJetIdx,               "ra2_loose_bestTopJetIdx/I");
  //  reducedValues->Branch("ra2_loose_pickedRemainingCombfatJetIdx",&m_loose_pickedRemainingCombfatJetIdx,"ra2_loose_pickedRemainingCombfatJetIdx/I");
  //
  //  reducedValues->Branch("ra2_loose_remainPassCSVS",    &m_loose_remainPassCSVS,    "ra2_loose_remainPassCSVS/O");
  //  reducedValues->Branch("ra2_passLooseTopTagger",      &m_passLooseTopTagger,      "ra2_passLooseTopTagger/O");
  //  reducedValues->Branch("ra2_passLooseTopJetIdx",      &m_passLooseTopJetIdx,      "ra2_passLooseTopJetIdx/O");
  //  reducedValues->Branch("ra2_passLooseTopMassCut",     &m_passLooseTopMassCut,     "ra2_passLooseTopMassCut/O");
  //  reducedValues->Branch("ra2_passLooseCSVCut",         &m_passLooseCSVCut,         "ra2_passLooseCSVCut/O");
  //  reducedValues->Branch("ra2_passLooseRemainingSystem",&m_passLooseRemainingSystem,"ra2_passLooseRemainingSystem/O");
  //  reducedValues->Branch("ra2_passLooseMT2Cuts",        &m_passLooseMT2Cuts,        "ra2_passLooseMT2Cuts/O");
  //
  //  reducedValues->Branch("ra2_loose_bestTopJetMass", &m_loose_bestTopJetMass, "ra2_loose_bestTopJetMass/D");
  //  reducedValues->Branch("ra2_loose_MTbJet",         &m_loose_MTbJet,         "ra2_loose_MTbJet/D");
  //  reducedValues->Branch("ra2_loose_MTbestTopJet",   &m_loose_MTbestTopJet,   "ra2_loose_MTbestTopJet/D");
  //  reducedValues->Branch("ra2_loose_MT2",            &m_loose_MT2,            "ra2_loose_MT2/D");
  //  reducedValues->Branch("ra2_loose_MTbestWJet",     &m_loose_MTbestWJet,     "ra2_loose_MTbestWJet/D");
  //  reducedValues->Branch("ra2_loose_MTbestbJet",     &m_loose_MTbestbJet,     "ra2_loose_MTbestbJet/D");
  //  reducedValues->Branch("ra2_loose_MTremainingTopJet",               &m_loose_MTremainingTopJet,               "ra2_loose_MTremainingTopJet/D");
  //  reducedValues->Branch("ra2_loose_linearCombMTbJetPlusMTbestTopJet",&m_loose_linearCombMTbJetPlusMTbestTopJet,"ra2_loose_linearCombMTbJetPlusMTbestTopJet/D");
  //  
  //  //////
  //  reducedValues->Branch("ra2_nominal_remainPassCSVS",    &m_nominal_remainPassCSVS,    "ra2_nominal_remainPassCSVS/O");
  //  reducedValues->Branch("ra2_passNominalTopTagger",      &m_passNominalTopTagger,      "ra2_passNominalTopTagger/O");
  //  reducedValues->Branch("ra2_passNominalTopJetIdx",      &m_passNominalTopJetIdx,      "ra2_passNominalTopJetIdx/O");
  //  reducedValues->Branch("ra2_passNominalTopMassCut",     &m_passNominalTopMassCut,     "ra2_passNominalTopMassCut/O");
  //  reducedValues->Branch("ra2_passNominalCSVCut",         &m_passNominalCSVCut,         "ra2_passNominalCSVCut/O");
  //  reducedValues->Branch("ra2_passNominalRemainingSystem",&m_passNominalRemainingSystem,"ra2_passNominalRemainingSystem/O");
  //  reducedValues->Branch("ra2_passNominalMT2Cuts",        &m_passNominalMT2Cuts,        "ra2_passNominalMT2Cuts/O");
  //
  //  reducedValues->Branch("ra2_nominal_remainPassCSVS",&m_nominal_remainPassCSVS,"ra2_nominal_remainPassCSVS/O");
  //  reducedValues->Branch("ra2_passNominalTopTagger",  &m_passNominalTopTagger,  "ra2_passNominalTopTagger/O");
  //
  //  reducedValues->Branch("ra2_nominal_bestTopJetMass", &m_nominal_bestTopJetMass, "ra2_nominal_bestTopJetMass/D");
  //  reducedValues->Branch("ra2_nominal_MTbJet",         &m_nominal_MTbJet,         "ra2_nominal_MTbJet/D");
  //  reducedValues->Branch("ra2_nominal_MTbestTopJet",   &m_nominal_MTbestTopJet,   "ra2_nominal_MTbestTopJet/D");
  //  reducedValues->Branch("ra2_nominal_MT2",            &m_nominal_MT2,            "ra2_nominal_MT2/D");
  //  reducedValues->Branch("ra2_nominal_MTbestWJet",     &m_nominal_MTbestWJet,     "ra2_nominal_MTbestWJet/D");
  //  reducedValues->Branch("ra2_nominal_MTbestbJet",     &m_nominal_MTbestbJet,     "ra2_nominal_MTbestbJet/D");
  //  reducedValues->Branch("ra2_nominal_MTremainingTopJet",               &m_nominal_MTremainingTopJet,               "ra2_nominal_MTremainingTopJet/D");
  //  reducedValues->Branch("ra2_nominal_linearCombMTbJetPlusMTbestTopJet",&m_nominal_linearCombMTbJetPlusMTbestTopJet,"ra2_nominal_linearCombMTbJetPlusMTbestTopJet/D");
  //}

  reducedValues->Branch("ra2_genBosons"      ,&m_genBosons      ,"ra2_genBosons/I" );
  reducedValues->Branch("ra2_genBoson1Pt"    ,&m_genBoson1Pt    ,"ra2_genBoson1Pt/D" );
  reducedValues->Branch("ra2_genBoson1Eta"   ,&m_genBoson1Eta   ,"ra2_genBoson1Eta/D" );
  reducedValues->Branch("ra2_genBoson1M"     ,&m_genBoson1M     ,"ra2_genBoson1M/D" );
  reducedValues->Branch("ra2_genBoson1MinDR" ,&m_genBoson1MinDR ,"ra2_genBoson1MinDR/D" );
  reducedValues->Branch("ra2_genBoson1DRJet1",&m_genBoson1DRJet1,"ra2_genBoson1DRJet1/D" );
  reducedValues->Branch("ra2_genHT"          ,&m_genHT          ,"ra2_genHT/D" );
  reducedValues->Branch("ra2_genMHT"         ,&m_genMHT         ,"ra2_genMHT/D");
  reducedValues->Branch("ra2_nGenJetsPt30Eta50",&m_nGenJetsPt30Eta50, "ra2_nGenJetsPt30Eta50/I" );
  reducedValues->Branch("ra2_nGenJetsPt50Eta25",&m_nGenJetsPt50Eta25, "ra2_nGenJetsPt50Eta25/I" );
  reducedValues->Branch("ra2_genDPhiMHT1",  &m_genDPhiMHT1,  "ra2_genDPhiMHT1/D");
  reducedValues->Branch("ra2_genDPhiMHT2",  &m_genDPhiMHT2,  "ra2_genDPhiMHT2/D");
  reducedValues->Branch("ra2_genDPhiMHT3",  &m_genDPhiMHT3,  "ra2_genDPhiMHT3/D");
  reducedValues->Branch("ra2_genDPhiMHT4",  &m_genDPhiMHT4,  "ra2_genDPhiMHT4/D");
  reducedValues->Branch("ra2_genDPhiMHTMin",&m_genDPhiMHTMin,"ra2_genDPhiMHTMin/D");
  //reducedValues->Branch("ra2_boson2Pt", &m_genBoson2Pt, "ra2_genBoson2Pt/D" );
  //reducedValues->Branch("ra2_boson2Eta",&m_genBoson2Eta,"ra2_genBoson2Eta/D" );
  //reducedValues->Branch("ra2_boson2M",  &m_genBoson2M,  "ra2_genBoson2M/D" );

  reducedValues->Branch("ra2_dPhiMHT1", &m_dPhiMHT1, "ra2_dPhiMHT1/D");
  reducedValues->Branch("ra2_dPhiMHT2", &m_dPhiMHT2, "ra2_dPhiMHT2/D");
  reducedValues->Branch("ra2_dPhiMHT3", &m_dPhiMHT3, "ra2_dPhiMHT3/D");
  reducedValues->Branch("ra2_dPhiMHT4", &m_dPhiMHT4, "ra2_dPhiMHT4/D");
  reducedValues->Branch("ra2_dPhiMHTMin", &m_dPhiMHTMin, "ra2_dPhiMHTMin/D");
  reducedValues->Branch("ra2_dPhiMHTMinBCSVM", &m_dPhiMHTMinBCSVM, "ra2_dPhiMHTMinBCSVM/D");
  reducedValues->Branch("ra2_dPhiMHTMinBCSVT", &m_dPhiMHTMinBCSVT, "ra2_dPhiMHTMinBCSVT/D");

  reducedValues->Branch("ra2_dPhiMET1", &m_dPhiMET1, "ra2_dPhiMET1/D");
  reducedValues->Branch("ra2_dPhiMET2", &m_dPhiMET2, "ra2_dPhiMET2/D");
  reducedValues->Branch("ra2_dPhiMET3", &m_dPhiMET3, "ra2_dPhiMET3/D");
  reducedValues->Branch("ra2_dPhiMET4", &m_dPhiMET4, "ra2_dPhiMET4/D");
  reducedValues->Branch("ra2_dPhiMETMin", &m_dPhiMETMin, "ra2_dPhiMETMin/D");
  reducedValues->Branch("ra2_dPhiMETMinBCSVM", &m_dPhiMETMinBCSVM, "ra2_dPhiMETMinBCSVM/D");
  reducedValues->Branch("ra2_dPhiMETMinBCSVT", &m_dPhiMETMinBCSVT, "ra2_dPhiMETMinBCSVT/D");

  reducedValues->Branch("ra2_Jet1Pt",  &m_Jet1Pt,  "ra2_Jet1Pt/D");
  reducedValues->Branch("ra2_Jet1Eta", &m_Jet1Eta, "ra2_Jet1Eta/D");
  reducedValues->Branch("ra2_Jet2Pt",  &m_Jet2Pt,  "ra2_Jet2Pt/D");
  reducedValues->Branch("ra2_Jet2Eta", &m_Jet2Eta, "ra2_Jet2Eta/D");
  reducedValues->Branch("ra2_Jet3Pt",  &m_Jet3Pt,  "ra2_Jet3Pt/D");
  reducedValues->Branch("ra2_Jet3Eta", &m_Jet3Eta, "ra2_Jet3Eta/D");
  reducedValues->Branch("ra2_Jet4Pt",  &m_Jet4Pt,  "ra2_Jet4Pt/D");
  reducedValues->Branch("ra2_Jet4Eta", &m_Jet4Eta, "ra2_Jet4Eta/D");

  reducedValues->Branch("ra2_JetCSVM1Pt",  &m_JetCSVM1Pt,  "ra2_JetCSVM1Pt/D");
  reducedValues->Branch("ra2_JetCSVM1Eta", &m_JetCSVM1Eta, "ra2_JetCSVM1Eta/D");
  reducedValues->Branch("ra2_JetCSVM2Pt",  &m_JetCSVM2Pt,  "ra2_JetCSVM2Pt/D");
  reducedValues->Branch("ra2_JetCSVM2Eta", &m_JetCSVM2Eta, "ra2_JetCSVM2Eta/D");
  reducedValues->Branch("ra2_JetCSVT1Pt",  &m_JetCSVT1Pt,  "ra2_JetCSVT1Pt/D");
  reducedValues->Branch("ra2_JetCSVT1Eta", &m_JetCSVT1Eta, "ra2_JetCSVT1Eta/D");
  reducedValues->Branch("ra2_JetCSVT2Pt",  &m_JetCSVT2Pt,  "ra2_JetCSVT2Pt/D");
  reducedValues->Branch("ra2_JetCSVT2Eta", &m_JetCSVT2Eta, "ra2_JetCSVT2Eta/D");

  reducedValues->Branch("ra2_PUWt",    &m_PUWt,    "ra2_PUWt/D");
  reducedValues->Branch("ra2_EventWt", &m_EventWt, "ra2_EventWt/D");

  reducedValues->Branch("ra2_nJetsCSVM", &m_nJetsCSVM, "ra2_nJetsCSVM/I");
  reducedValues->Branch("ra2_nJetsCSVT", &m_nJetsCSVT, "ra2_nJetsCSVT/I");
  reducedValues->Branch("ra2_nJetsPt30Eta50", &m_nJetsPt30Eta50, "ra2_nJetsPt30Eta50/I" );
  reducedValues->Branch("ra2_nJetsPt30Eta24", &m_nJetsPt30Eta24, "ra2_nJetsPt30Eta24/I");
  reducedValues->Branch("ra2_nJetsPt50Eta24", &m_nJetsPt50Eta24, "ra2_nJetsPt50Eta24/I");
  reducedValues->Branch("ra2_nJetsPt70Eta24", &m_nJetsPt70Eta24, "ra2_nJetsPt70Eta24/I");
  reducedValues->Branch("ra2_nJetsPt50Eta25", &m_nJetsPt50Eta25, "ra2_nJetsPt50Eta25/I" );
  reducedValues->Branch("ra2_nJetsPt50Eta25MInv", &m_nJetsPt50Eta25MInv, "ra2_nJetsPt50Eta25MInv/I" );

  //if (storeExtraVetos_) {
  reducedValues->Branch("ra2_passRA2ElVeto",    &m_passRA2ElVeto   , "ra2_passRA2ElVeto/O"    );
  reducedValues->Branch("ra2_passRA2MuVeto",    &m_passRA2MuVeto   , "ra2_passRA2MuVeto/O"    );
  reducedValues->Branch("ra2_passDirIsoElVeto", &m_passDirIsoElVeto, "ra2_passDirIsoElVeto/O"    );
  reducedValues->Branch("ra2_passDirIsoMuVeto", &m_passDirIsoMuVeto, "ra2_passDirIsoMuVeto/O"    );
  reducedValues->Branch("ra2_passIsoTrkVeto",   &m_passIsoTrkVeto  , "ra2_passIsoTrkVeto/O");
  //}

  reducedValues->SetAutoSave(1);
}


void  RA2ZInvTreeMaker::beginRun(edm::Run const&, edm::EventSetup const&) {
}

void RA2ZInvTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void RA2ZInvTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RA2ZInvTreeMaker);

