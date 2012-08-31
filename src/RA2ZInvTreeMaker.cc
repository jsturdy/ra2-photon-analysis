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
// $Id: RA2ZInvTreeMaker.cc,v 1.1 2012/08/30 09:43:40 sturdy Exp $
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
  genLabel_       = pset.getParameter<edm::InputTag>("genLabel");  
  vertexSrc_      = pset.getParameter<edm::InputTag>("VertexSrc");
  jetSrc_         = pset.getParameter<edm::InputTag>("JetSrc");
  bJetSrc_        = pset.getParameter<edm::InputTag>("bJetSrc");
  htJetSrc_       = pset.getParameter<edm::InputTag>("htJetSrc");
  htSrc_          = pset.getParameter<edm::InputTag>("htSource");
  mhtSrc_         = pset.getParameter<edm::InputTag>("mhtSource");
  doPUReWeight_   = pset.getParameter<bool>("DoPUReweight");
  puWeightSrc_    = pset.getParameter<edm::InputTag>("PUWeightSource");


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

  edm::Handle<reco::GenParticleCollection> gens;
  ev.getByLabel(genLabel_,gens);

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

  m_EventWt = scale_;
  m_PUWt    = pu_event_wt;
  m_Vertices = vertices->size();

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
  m_boson1Pt  = -10.;
  m_boson1Eta = -10.;
  m_boson1M   = -10.;
  m_boson2Pt  = -10.;
  m_boson2Eta = -10.;
  m_boson2M   = -10.;
  m_nBosons  = gens->size();
  if (gens->size() > 0) {
    m_boson1Pt  = (*gens)[0].pt();
    m_boson1Eta = (*gens)[0].eta();
    m_boson1M   = (*gens)[0].mass();
    if (gens->size() > 1) {
      m_boson2Pt  = (*gens)[1].pt();
      m_boson2Eta = (*gens)[1].eta();
      m_boson2M   = (*gens)[1].mass();
    }
  }

  ///////////
  m_nJetsPt30Eta50 = jets  ->size();
  m_nJetsPt50Eta25 = htJets->size();
  m_bJetsPt30Eta24 = bJets ->size();
  m_HT  = *ht;
  m_MHT = (*mht)[0].pt();

  const pat::Jet  *r1, *r2, *r3;
  m_dPhi1 = -1.0;
  m_dPhi2 = -1.0;
  m_dPhi3 = -1.0;
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
  reducedValues->Branch("ra2_boson1Pt", &m_boson1Pt, "m_boson1Pt/D" );
  reducedValues->Branch("ra2_boson1Eta",&m_boson1Eta,"m_boson1Eta/D" );
  reducedValues->Branch("ra2_boson1M",  &m_boson1M,  "m_boson1M/D" );
  reducedValues->Branch("ra2_boson2Pt", &m_boson2Pt, "m_boson2Pt/D" );
  reducedValues->Branch("ra2_boson2Eta",&m_boson2Eta,"m_boson2Eta/D" );
  reducedValues->Branch("ra2_boson2M",  &m_boson2M,  "m_boson2M/D" );

  reducedValues->Branch("ra2_dPhi1", &m_dPhi1, "m_dPhi1/D");
  reducedValues->Branch("ra2_dPhi2", &m_dPhi2, "m_dPhi2/D");
  reducedValues->Branch("ra2_dPhi3", &m_dPhi3, "m_dPhi3/D");

  reducedValues->Branch("ra2_Jet1Pt",  &m_Jet1Pt,  "m_Jet1Pt/D");
  reducedValues->Branch("ra2_Jet1Eta", &m_Jet1Eta, "m_Jet1Eta/D");
  reducedValues->Branch("ra2_Jet2Pt",  &m_Jet2Pt,  "m_Jet2Pt/D");
  reducedValues->Branch("ra2_Jet2Eta", &m_Jet2Eta, "m_Jet2Eta/D");
  reducedValues->Branch("ra2_Jet3Pt",  &m_Jet3Pt,  "m_Jet3Pt/D");
  reducedValues->Branch("ra2_Jet3Eta", &m_Jet3Eta, "m_Jet3Eta/D");

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

