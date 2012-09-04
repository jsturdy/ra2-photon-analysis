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
// $Id: RA2ZInvDiMuonTreeMaker.cc,v 1.2 2012/08/31 10:27:22 sturdy Exp $
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
  debug_         = pset.getParameter<bool>("Debug");
  data_          = pset.getParameter<bool>("Data");
  scale_         = pset.getParameter<double>("ScaleFactor");
  muonSrc_       = pset.getParameter<edm::InputTag>("MuonSrc");
  vertexSrc_     = pset.getParameter<edm::InputTag>("VertexSrc");
  jetSrc_        = pset.getParameter<edm::InputTag>("JetSrc");
  bJetSrc_       = pset.getParameter<edm::InputTag>("bJetSrc");
  htJetSrc_      = pset.getParameter<edm::InputTag>("htJetSrc");
  htSrc_         = pset.getParameter<edm::InputTag>("htSource");
  mhtSrc_        = pset.getParameter<edm::InputTag>("mhtSource");
  doPUReWeight_  = pset.getParameter<bool>("DoPUReweight");
  puWeightSrc_   = pset.getParameter<edm::InputTag>("PUWeightSource");
  
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

  // if MC, do PU reweighting
  double pu_event_wt = 1.0;
  edm::Handle<double> puweight;
  if( doPUReWeight_ ) {
    ev.getByLabel(puWeightSrc_, puweight);
    pu_event_wt = *puweight;
  }

  m_EventWt  = scale_;
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
      std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
    }
  }
  //////

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

  const pat::Jet  *r1, *r2, *r3;
  m_dPhi1 = 10.0;
  m_dPhi2 = 10.0;
  m_dPhi3 = 10.0;
  m_Jet1Pt  = -10.;
  m_Jet1Eta = -10.;
  m_Jet2Pt  = -10.;
  m_Jet2Eta = -10.;
  m_Jet3Pt  = -10.;
  m_Jet3Eta = -10.;
  if(m_nJetsPt30Eta50 >= 1) {
    r1 = &((*jets)[0]);
    m_dPhi1 = fabs(reco::deltaPhi(r1->phi(),(*mht)[0].phi()));
    m_Jet1Pt = r1->pt();
    m_Jet1Eta = r1->eta();
    if(m_nJetsPt30Eta50 >= 2) {
      r2 = &((*jets)[1]);
      m_dPhi2 = fabs(reco::deltaPhi(r2->phi(),(*mht)[0].phi())); 
      m_Jet2Pt = r2->pt();
      m_Jet2Eta = r2->eta();
     
      if(m_nJetsPt30Eta50 >= 3) {
	r3 = &((*jets)[2]);
	m_dPhi3 = fabs(reco::deltaPhi(r3->phi(),(*mht)[0].phi()));
	m_Jet3Pt = r3->pt();
	m_Jet3Eta = r3->eta();
      }
    }
  }

  m_dPhiMin  = 10.;
  m_nJetsPt30Eta24 = 0;
  edm::View<pat::Jet>::const_iterator jet = jets->begin();
  for (; jet!= jets->end(); ++jet) {
    if (jet->pt() > 30 && fabs(jet->eta() < 2.4))
      ++m_nJetsPt30Eta24;
    double tmpDPhi = fabs(reco::deltaPhi(jet->phi(),(*mht)[0].phi()));
    if (tmpDPhi < m_dPhiMin)
      m_dPhiMin = tmpDPhi;
  }
  m_dPhiMinB = 10.;
  edm::View<pat::Jet>::const_iterator bjet = bJets->begin();
  for (; bjet!= bJets->end(); ++bjet) {
    double tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),(*mht)[0].phi()));
    if (tmpDPhiB < m_dPhiMinB)
      m_dPhiMinB = tmpDPhiB;
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

  reducedValues->Branch("ra2_HT",  &m_HT,  "m_HT/D" );
  reducedValues->Branch("ra2_MHT", &m_MHT, "m_MHT/D");
  reducedValues->Branch("ra2_Vertices", &m_Vertices, "m_Vertices/I");

  reducedValues->Branch("ra2_dPhi1", &m_dPhi1, "m_dPhi1/D");
  reducedValues->Branch("ra2_dPhi2", &m_dPhi2, "m_dPhi2/D");
  reducedValues->Branch("ra2_dPhi3", &m_dPhi3, "m_dPhi3/D");
  reducedValues->Branch("ra2_dPhiMin", &m_dPhiMin, "m_dPhiMin/D");
  reducedValues->Branch("ra2_dPhiMinB", &m_dPhiMinB, "m_dPhiMinB/D");

  reducedValues->Branch("ra2_Jet1Pt",  &m_Jet1Pt,  "m_Jet1Pt/D");
  reducedValues->Branch("ra2_Jet1Eta", &m_Jet1Eta, "m_Jet1Eta/D");
  reducedValues->Branch("ra2_Jet2Pt",  &m_Jet2Pt,  "m_Jet2Pt/D");
  reducedValues->Branch("ra2_Jet2Eta", &m_Jet2Eta, "m_Jet2Eta/D");
  reducedValues->Branch("ra2_Jet3Pt",  &m_Jet3Pt,  "m_Jet3Pt/D");
  reducedValues->Branch("ra2_Jet3Eta", &m_Jet3Eta, "m_Jet3Eta/D");

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

