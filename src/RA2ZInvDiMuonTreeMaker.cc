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
// $Id: RA2ZInvDiMuonTreeMaker.cc,v 1.4 2012/10/29 10:46:08 sturdy Exp $
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
  topTaggerSrc_   = pset.getParameter<std::string>("topTaggerSource");
  doPUReWeight_   = pset.getParameter<bool>("DoPUReweight");
  puWeightSrc_    = pset.getParameter<edm::InputTag>("PUWeightSource");
  eventWeightSrc_ = pset.getParameter< edm::InputTag >( "EventWeightSource" );
  
  //key to help getting the hlt process from event provenance
  checkedProcess_ = false;
  processName_    = "";

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

  // get muons 
  edm::Handle< std::vector<pat::Muon> > patMuons;
  ev.getByLabel(muonSrc_, patMuons); 
  if (debug_)
    std::cout<<"Muon collection has size "<<patMuons->size()<<std::endl;

  // get jetcollection
  edm::Handle<reco::VertexCollection > vertices;
  ev.getByLabel(vertexSrc_, vertices);
  if (debug_)
    std::cout<<"vertex collection has size "<<vertices->size()<<std::endl;

  edm::Handle<edm::View<pat::Jet> > jets;
  ev.getByLabel(jetSrc_, jets);
  if (debug_)
    std::cout<<"Jet collection has size "<<jets->size()<<std::endl;

  edm::Handle<edm::View<pat::Jet> > htJets;
  ev.getByLabel(htJetSrc_, htJets);
  if (debug_)
    std::cout<<"HT Jet collection has size "<<htJets->size()<<std::endl;

  edm::Handle<edm::View<pat::Jet> > bJets;
  ev.getByLabel(bJetSrc_, bJets);
  if (debug_)
    std::cout<<"b-Jet collection has size "<<bJets->size()<<std::endl;

  edm::Handle<double > ht;
  ev.getByLabel(htSrc_, ht);
  if (debug_)
    std::cout<<"HT value "<<*ht<<std::endl;

  edm::Handle<edm::View<reco::MET> > mht;
  ev.getByLabel(mhtSrc_, mht);
  if (debug_)
    std::cout<<"MHT value "<<(*mht)[0].pt()<<std::endl;

  edm::Handle<edm::View<reco::MET> > met;
  ev.getByLabel(metSrc_, met);
  if (debug_)
    std::cout<<"MET value "<<(*met)[0].pt()<<std::endl;

  ///top tagger variables
  edm::Handle<double > hbestTopJetMass;
  ev.getByLabel(topTaggerSrc_,"bestTopJetMass", hbestTopJetMass);
  if (debug_)
    std::cout<<"hbestTopJetMass value "<<*hbestTopJetMass<<std::endl;

  edm::Handle<double > hTbJet;
  ev.getByLabel(topTaggerSrc_,"mTbJet", hTbJet);
  if (debug_)
    std::cout<<"hTbJet value "<<*hTbJet<<std::endl;

  edm::Handle<double > hTbestTopJet;
  ev.getByLabel(topTaggerSrc_,"mTbestTopJet", hTbestTopJet);
  if (debug_)
    std::cout<<"hTbestTopJet value "<<*hTbestTopJet<<std::endl;

  edm::Handle<double > hMT2;
  ev.getByLabel(topTaggerSrc_,"MT2", hMT2);
  if (debug_)
    std::cout<<"hMT2 value "<<*hMT2<<std::endl;


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

  //std::vector<const pat::Muon*> RA2Muons;
  //double muon_ptcut  = 15.0;

   //////
  if(debug_ && patMuons->size() > 0) {
    std::cout << "Isolated muons : " << std::endl;
    for( unsigned int imuon=0; imuon<patMuons->size(); imuon++) {
      std::cout << imuon << " " <<(*patMuons)[imuon].pt() 
		<< " " << (*patMuons)[imuon].eta()
		<< " " << (*patMuons)[imuon].phi() 
		<< std::endl;
    }
    std::cout << "Jets("<<jets->size()<<") : " << std::endl;
    for(unsigned int i=0; i<jets->size(); i++) {
      const pat::Jet *r = &((*jets)[i]);
      std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
    }
    std::cout << "bJets("<<bJets->size()<<") : " << std::endl;
    for(unsigned int i=0; i<bJets->size(); i++) {
      const pat::Jet *r = &((*bJets)[i]);
      std::cout << i << ": pt::" << r->pt()<<": eta::"<<r->eta()<<": phi::"<<r->phi()<<": csv::"<<r->bDiscriminator("combinedSecondaryVertexBJetTags")<<std::endl;
    }
  }
  //////
  if (patMuons->size() < 2)
    return;

  m_nMuonsIso = patMuons->size();
  m_Muon1Pt   = (*patMuons)[0].pt();
  m_Muon1Eta  = (*patMuons)[0].eta();
  m_Muon2Pt   = (*patMuons)[1].pt();
  m_Muon2Eta  = (*patMuons)[1].eta();
  
  m_DiMuonInvM = ((*patMuons)[0].p4()+(*patMuons)[1].p4()).mass();
  m_DiMuonPt   = ((*patMuons)[0].p4()+(*patMuons)[1].p4()).pt();

  m_nJetsPt30Eta50 = jets  ->size();
  m_nJetsPt50Eta25 = htJets->size();
  m_bJetsPt30Eta24 = bJets ->size();
  m_HT  = *ht;
  m_MHT = (*mht)[0].pt();
  m_MET = (*met)[0].pt();

  m_modMET = ((*met)[0].p4()+(*patMuons)[0].p4()+(*patMuons)[1].p4()).pt();

  m_bestTopJetMass = *hbestTopJetMass;
  m_TbJet          = *hTbJet;
  m_TbestTopJet    = *hTbestTopJet;
  m_MT2            = *hMT2;

  const pat::Jet  *r1, *r2, *r3, *r4;
  m_dPhiMHT1 = 10.0;  m_dPhiMET1 = 10.0;  m_dPhiModMET1 = 10.0;
  m_dPhiMHT2 = 10.0;  m_dPhiMET2 = 10.0;  m_dPhiModMET2 = 10.0;
  m_dPhiMHT3 = 10.0;  m_dPhiMET3 = 10.0;  m_dPhiModMET3 = 10.0;
  m_dPhiMHT4 = 10.0;  m_dPhiMET4 = 10.0;  m_dPhiModMET4 = 10.0;
  m_Jet1Pt  = -10.;   m_Jet3Pt  = -10.;
  m_Jet1Eta = -10.;   m_Jet3Eta = -10.;
  m_Jet2Pt  = -10.;   m_Jet4Pt  = -10.;
  m_Jet2Eta = -10.;   m_Jet4Eta = -10.;

  if(m_nJetsPt30Eta50 >= 1) {
    r1 = &((*jets)[0]);
    m_dPhiMHT1 = fabs(reco::deltaPhi(r1->phi(),(*mht)[0].phi()));
    m_dPhiMET1 = fabs(reco::deltaPhi(r1->phi(),(*met)[0].phi()));
    m_dPhiModMET1 = fabs(reco::deltaPhi(r1->phi(),((*met)[0].p4()+(*patMuons)[0].p4()+(*patMuons)[1].p4()).phi()));
    m_Jet1Pt = r1->pt();
    m_Jet1Eta = r1->eta();
    if(m_nJetsPt30Eta50 >= 2) {
      r2 = &((*jets)[1]);
      m_dPhiMHT2 = fabs(reco::deltaPhi(r2->phi(),(*mht)[0].phi())); 
      m_dPhiMET2 = fabs(reco::deltaPhi(r2->phi(),(*met)[0].phi())); 
      m_dPhiModMET2 = fabs(reco::deltaPhi(r2->phi(),((*met)[0].p4()+(*patMuons)[0].p4()+(*patMuons)[1].p4()).phi()));
      m_Jet2Pt = r2->pt();
      m_Jet2Eta = r2->eta();
     
      if(m_nJetsPt30Eta50 >= 3) {
	r3 = &((*jets)[2]);
	m_dPhiMHT3 = fabs(reco::deltaPhi(r3->phi(),(*mht)[0].phi()));
	m_dPhiMET3 = fabs(reco::deltaPhi(r3->phi(),(*met)[0].phi()));
	m_dPhiModMET3 = fabs(reco::deltaPhi(r3->phi(),((*met)[0].p4()+(*patMuons)[0].p4()+(*patMuons)[1].p4()).phi()));
	m_Jet3Pt = r3->pt();
	m_Jet3Eta = r3->eta();

	if(m_nJetsPt30Eta50 >= 4) {
	  r4 = &((*jets)[3]);
	  m_dPhiMHT4 = fabs(reco::deltaPhi(r4->phi(),(*mht)[0].phi()));
	  m_dPhiMET4 = fabs(reco::deltaPhi(r4->phi(),(*met)[0].phi()));
	  m_dPhiModMET4 = fabs(reco::deltaPhi(r4->phi(),((*met)[0].p4()+(*patMuons)[0].p4()+(*patMuons)[1].p4()).phi()));
	  m_Jet4Pt = r4->pt();
	  m_Jet4Eta = r4->eta();
	}
      }
    }
  }

  m_dPhiMHTMin  = 10.;
  m_dPhiMETMin  = 10.;
  m_dPhiModMETMin  = 10.;
  m_nJetsPt30Eta24 = 0;
  m_nJetsPt50Eta25MInv = 0;
  m_HTMInv = 0.;
  edm::View<pat::Jet>::const_iterator jet = jets->begin();
  for (; jet!= jets->end(); ++jet) {
    if (jet->pt() > 30 && fabs(jet->eta() < 2.4))
      ++m_nJetsPt30Eta24;
    if (jet->pt() > 50 && fabs(jet->eta() < 2.5)) 
      if (patMuons->size()>1)
	if ((jet->p4()+((*patMuons)[0].p4()+(*patMuons)[1].p4())).mass() > 90.0) {
	  ++m_nJetsPt50Eta25MInv;
	  m_HTMInv += jet->pt();
	}
    
    double tmpDPhi = fabs(reco::deltaPhi(jet->phi(),(*mht)[0].phi()));
    if (tmpDPhi < m_dPhiMHTMin)
      m_dPhiMHTMin = tmpDPhi;

    tmpDPhi = fabs(reco::deltaPhi(jet->phi(),(*met)[0].phi()));
    if (tmpDPhi < m_dPhiMETMin)
      m_dPhiMETMin = tmpDPhi;

    tmpDPhi = fabs(reco::deltaPhi(jet->phi(),((*met)[0].p4()+(*patMuons)[0].p4()+(*patMuons)[1].p4()).phi()));
    if (tmpDPhi < m_dPhiModMETMin)
      m_dPhiModMETMin = tmpDPhi;
  }

  m_dPhiMHTMinB = 10.;
  m_dPhiMETMinB = 10.;
  m_dPhiModMETMinB = 10.;
  edm::View<pat::Jet>::const_iterator bjet = bJets->begin();
  for (; bjet!= bJets->end(); ++bjet) {
    double tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),(*mht)[0].phi()));
    if (tmpDPhiB < m_dPhiMHTMinB)
      m_dPhiMHTMinB = tmpDPhiB;

    tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),(*met)[0].phi()));
    if (tmpDPhiB < m_dPhiMETMinB)
      m_dPhiMETMinB = tmpDPhiB;
    tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),((*met)[0].p4()+(*patMuons)[0].p4()+(*patMuons)[1].p4()).phi()));
    if (tmpDPhiB < m_dPhiModMETMinB)
      m_dPhiModMETMinB = tmpDPhiB;
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
	ev.getByLabel("hltTriggerSummaryAOD", hltEventHandle);
	processName_ = hltEventHandle.provenance()->processName();
	if (debug_)
	  std::cout<<processName_<<std::endl;
      }
    hlTriggerResults_ = InputTag("TriggerResults","",processName_);
    
    edm::LogInfo("HLTEventSelector") << "Using trigger results for InputTag " << hlTriggerResults_;
    
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
  reducedValues->Branch("ra2_modMET",   &m_modMET,   "m_modMET/D");
  reducedValues->Branch("ra2_Vertices", &m_Vertices, "m_Vertices/I");

  reducedValues->Branch("ra2_bestTopJetMass", &m_bestTopJetMass, "m_bestTopJetMass/D");
  reducedValues->Branch("ra2_TbJet",          &m_TbJet,          "m_TbJet/D");
  reducedValues->Branch("ra2_TbestTopJet",    &m_TbestTopJet,    "m_TbestTopJet/D");
  reducedValues->Branch("ra2_MT2",            &m_MT2,            "m_MT2/D");

  reducedValues->Branch("ra2_dPhiMHT1", &m_dPhiMHT1, "m_dPhiMHT1/D");
  reducedValues->Branch("ra2_dPhiMHT2", &m_dPhiMHT2, "m_dPhiMHT2/D");
  reducedValues->Branch("ra2_dPhiMHT3", &m_dPhiMHT3, "m_dPhiMHT3/D");
  reducedValues->Branch("ra2_dPhiMHT4", &m_dPhiMHT4, "m_dPhiMHT4/D");
  reducedValues->Branch("ra2_dPhiMHTMin", &m_dPhiMHTMin, "m_dPhiMHTMin/D");
  reducedValues->Branch("ra2_dPhiMHTMinB", &m_dPhiMHTMinB, "m_dPhiMHTMinB/D");

  reducedValues->Branch("ra2_dPhiMET1", &m_dPhiMET1, "m_dPhiMET1/D");
  reducedValues->Branch("ra2_dPhiMET2", &m_dPhiMET2, "m_dPhiMET2/D");
  reducedValues->Branch("ra2_dPhiMET3", &m_dPhiMET3, "m_dPhiMET3/D");
  reducedValues->Branch("ra2_dPhiMET4", &m_dPhiMET4, "m_dPhiMET4/D");
  reducedValues->Branch("ra2_dPhiMETMin", &m_dPhiMETMin, "m_dPhiMETMin/D");
  reducedValues->Branch("ra2_dPhiMETMinB", &m_dPhiMETMinB, "m_dPhiMETMinB/D");

  reducedValues->Branch("ra2_dPhiModMET1", &m_dPhiModMET1, "m_dPhiModMET1/D");
  reducedValues->Branch("ra2_dPhiModMET2", &m_dPhiModMET2, "m_dPhiModMET2/D");
  reducedValues->Branch("ra2_dPhiModMET3", &m_dPhiModMET3, "m_dPhiModMET3/D");
  reducedValues->Branch("ra2_dPhiModMET4", &m_dPhiModMET4, "m_dPhiModMET4/D");
  reducedValues->Branch("ra2_dPhiModMETMin", &m_dPhiModMETMin, "m_dPhiModMETMin/D");
  reducedValues->Branch("ra2_dPhiModMETMinB", &m_dPhiModMETMinB, "m_dPhiModMETMinB/D");

  reducedValues->Branch("ra2_Jet1Pt",  &m_Jet1Pt,  "m_Jet1Pt/D");
  reducedValues->Branch("ra2_Jet1Eta", &m_Jet1Eta, "m_Jet1Eta/D");
  reducedValues->Branch("ra2_Jet2Pt",  &m_Jet2Pt,  "m_Jet2Pt/D");
  reducedValues->Branch("ra2_Jet2Eta", &m_Jet2Eta, "m_Jet2Eta/D");
  reducedValues->Branch("ra2_Jet3Pt",  &m_Jet3Pt,  "m_Jet3Pt/D");
  reducedValues->Branch("ra2_Jet3Eta", &m_Jet3Eta, "m_Jet3Eta/D");
  reducedValues->Branch("ra2_Jet4Pt",  &m_Jet4Pt,  "m_Jet4Pt/D");
  reducedValues->Branch("ra2_Jet4Eta", &m_Jet4Eta, "m_Jet4Eta/D");

  reducedValues->Branch("ra2_PUWt",    &m_PUWt,    "m_PUWt/D");
  reducedValues->Branch("ra2_EventWt", &m_EventWt, "m_EventWt/D");

  reducedValues->Branch("ra2_nMuonsIso",  &m_nMuonsIso,  "m_nMuonsIso/I");
  reducedValues->Branch("ra2_Muon1Pt",    &m_Muon1Pt,    "m_Muon1Pt/D" );
  reducedValues->Branch("ra2_Muon1Eta",   &m_Muon1Eta,   "m_Muon1Eta/D");
  reducedValues->Branch("ra2_Muon2Pt",    &m_Muon2Pt,    "m_Muon2Pt/D" );
  reducedValues->Branch("ra2_Muon2Eta",   &m_Muon2Eta,   "m_Muon2Eta/D");
  reducedValues->Branch("ra2_DiMuonInvM", &m_DiMuonInvM, "m_DiMuonInvM/D");
  reducedValues->Branch("ra2_DiMuonPt",   &m_DiMuonPt,   "m_DiMuonPt/D");

  reducedValues->Branch("ra2_Mu13_Mu8",         &m_Mu13_Mu8        , "m_Mu13_Mu8/O"        );
  reducedValues->Branch("ra2_Mu17_Mu8",         &m_Mu17_Mu8        , "m_Mu17_Mu8/O"        );
  reducedValues->Branch("ra2_DoubleMu5_IsoMu5", &m_DoubleMu5_IsoMu5, "m_DoubleMu5_IsoMu5/O");

  reducedValues->Branch("ra2_nJetsPt30Eta50", &m_nJetsPt30Eta50, "m_nJetsPt30Eta50/I" );
  reducedValues->Branch("ra2_nJetsPt30Eta24", &m_nJetsPt30Eta24, "m_nJetsPt30Eta24/I");
  reducedValues->Branch("ra2_bJetsPt30Eta24", &m_bJetsPt30Eta24, "m_bJetsPt30Eta24/I");
  reducedValues->Branch("ra2_nJetsPt50Eta25", &m_nJetsPt50Eta25, "m_nJetsPt50Eta25/I" );
  reducedValues->Branch("ra2_nJetsPt50Eta25MInv", &m_nJetsPt50Eta25MInv, "m_nJetsPt50Eta25MInv/I" );

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

