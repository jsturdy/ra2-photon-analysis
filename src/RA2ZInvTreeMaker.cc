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
// $Id: RA2ZInvTreeMaker.cc,v 1.1 2012/07/20 11:35:34 sturdy Exp $
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
  data_           = pset.getParameter<bool>("Data");
  scale_          = pset.getParameter<double>("ScaleFactor");
  vertexSrc_      = pset.getParameter<edm::InputTag>("VertexSrc");
  jetSrc_         = pset.getParameter<edm::InputTag>("JetSrc");
  bJetSrc_        = pset.getParameter<edm::InputTag>("bJetSrc");
  jetHTSrc_       = pset.getParameter<edm::InputTag>("JetHTSource");
  doPUReWeight_   = pset.getParameter<bool>("DoPUReweight");
  puWeightSrc_    = pset.getParameter<edm::InputTag>("PUWeightSource");

  if (!data_)
    genLabel_   = pset.getParameter<edm::InputTag>("genLabel");  

}

RA2ZInvTreeMaker::~RA2ZInvTreeMaker() {
  delete reducedValues;
  reducedValues = 0;
}

void RA2ZInvTreeMaker::analyze(const edm::Event& ev, const edm::EventSetup& es) {

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

  edm::Handle<reco::GenParticleCollection> gens;
  m_bosonPt  = 0;
  m_bosonEta = 0;
  m_bosonM   = 0;
  if (!data_) {
    ev.getByLabel(genLabel_,gens);
    m_nBosons  = gens->size();
    if (gens->size() > 0) {
      m_bosonPt  = (*gens)[0].pt();
      m_bosonEta = (*gens)[0].eta();
      m_bosonM   = (*gens)[0].mass();
    }
  }
  edm::Handle<reco::VertexCollection > vertices;
  ev.getByLabel(vertexSrc_, vertices);

  // get jetcollection
  edm::Handle<edm::View<pat::Jet> > jets;
  ev.getByLabel(jetSrc_, jets);

  edm::Handle<edm::View<pat::Jet> > jetsHT;
  ev.getByLabel(jetHTSrc_, jetsHT);

  // get jetcollection
  edm::Handle<edm::View<pat::Jet> > bJets;
  ev.getByLabel(bJetSrc_, bJets);

  // if MC, do PU reweighting
  double pu_event_wt = 1.0;
  edm::Handle<double> puweight;
  if( doPUReWeight_ ) {
    ev.getByLabel(puWeightSrc_, puweight);
    pu_event_wt = *puweight;
  }

  m_EventWt = scale_;
  m_PUWt    = pu_event_wt;
  m_Vertices = vertices->size();

  m_nJetsPt30Eta50 = 0;
  m_nJetsPt50Eta25 = 0;
  m_HT  = 0.0;
  reco::MET::LorentzVector mht(0,0,0,0);
  for(unsigned int i=0; i<jets->size(); ++i) {
    const pat::Jet r = (*jets)[i];
    if((*jets)[i].pt() > 50.0 && fabs((*jets)[i].eta()) < 2.50) {
      m_nJetsPt50Eta25++;
      m_HT  += (*jets)[i].pt();
    }
    if((*jets)[i].pt() > 30.0 && fabs((*jets)[i].eta()) < 5.0) { 
      m_nJetsPt30Eta50++;
      mht  -= (*jets)[i].p4();
    }
  }

  ///Count the b-jets
  m_bJetsPt30Eta24 = 0;
  m_bJetsPt30Eta24 = bJets->size();

  reco::MET MHT = reco::MET(mht, reco::MET::Point());
  m_MHT =  MHT.pt();

  m_dPhi1 = -1.0;
  m_dPhi2 = -1.0;
  m_dPhi3 = -1.0;
  if(m_nJetsPt30Eta50 >= 1) {
    m_dPhi1 = fabs(reco::deltaPhi((*jets)[0].phi(),MHT.phi()));
    if(m_nJetsPt30Eta50 >= 2) {
      m_dPhi2 = fabs(reco::deltaPhi((*jets)[1].phi(),MHT.phi()));
      if(m_nJetsPt30Eta50 >= 3) 
	m_dPhi3 = fabs(reco::deltaPhi((*jets)[2].phi(),MHT.phi()));
    }
  }

  //if(n_jets_pt50Eta25 >= ra2NJets_) {
  //  //if( ht> 500.0) {
  //  if( ht> ra2HT_) {
  //    //if(mht_val > 200.0) {
  //    if(mht_val > ra2MHT_) {
  //	
  //	if( !ra2ApplyDphiCuts_ || (ra2ApplyDphiCuts_ && (dphi_mht_j1 > 0.5 && dphi_mht_j2 > 0.5 && dphi_mht_j3 > 0.3)) ) {
  //	  
  //	} 
  //    }
  //  }
  //}
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

  reducedValues->Branch("ra2_HT",  &m_HT,  "m_HT/D" );
  reducedValues->Branch("ra2_MHT", &m_MHT, "m_MHT/D");
  reducedValues->Branch("ra2_Vertices", &m_Vertices, "m_Vertices/I");

  reducedValues->Branch("ra2_nBosons",  &m_nBosons,  "m_nBosons/I" );
  reducedValues->Branch("ra2_bosonPt",  &m_bosonPt,  "m_bosonPt/D" );
  reducedValues->Branch("ra2_bosonEta", &m_bosonEta, "m_bosonEta/D" );
  reducedValues->Branch("ra2_bosonM",   &m_bosonM,   "m_bosonM/D" );

  reducedValues->Branch("ra2_dPhi1", &m_dPhi1, "m_dPhi1/D");
  reducedValues->Branch("ra2_dPhi2", &m_dPhi2, "m_dPhi2/D");
  reducedValues->Branch("ra2_dPhi3", &m_dPhi3, "m_dPhi3/D");

  reducedValues->Branch("ra2_PUWt",    &m_PUWt,    "m_PUWt/D");
  reducedValues->Branch("ra2_EventWt", &m_EventWt, "m_EventWt/D");

  reducedValues->Branch("ra2_nJetsPt30Eta50", &m_nJetsPt30Eta50, "nJetsPt30Eta50/I" );
  reducedValues->Branch("ra2_bJetsPt30Eta24", &m_bJetsPt30Eta24, "bJetsPt30Eta24/I");
  reducedValues->Branch("ra2_nJetsPt50Eta25", &m_nJetsPt50Eta25, "nJetsPt50Eta25/I" );

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

