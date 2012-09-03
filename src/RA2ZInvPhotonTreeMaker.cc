// -*- C++ -*-
//
// Package:    RA2ZInvPhotonTreeMaker
// Class:      RA2ZInvPhotonTreeMaker
// 
/**\class RA2ZInvPhotonTreeMaker RA2ZInvPhotonTreeMaker.cc SusyAnalysis/RA2ZInvPhotonTreeMaker/src/RA2ZInvPhotonTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seema Sharma
//         Created:  Mon Jun 20 12:58:08 CDT 2011
// $Id: RA2ZInvPhotonTreeMaker.cc,v 1.2 2012/08/31 10:27:22 sturdy Exp $
//
//


// system include files
#include <memory>
#include <iomanip>
#include <cmath>

#include "ZInvisibleBkgds/Photons/interface/RA2ZInvPhotonTreeMaker.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include <DataFormats/METReco/interface/MET.h>


RA2ZInvPhotonTreeMaker::RA2ZInvPhotonTreeMaker(const edm::ParameterSet& pset) {

  // read parameters from config file
  debug_          = pset.getParameter<bool>("Debug");
  data_           = pset.getParameter<bool>("Data");
  scale_          = pset.getParameter<double>("ScaleFactor");
  photonSrc_      = pset.getParameter<edm::InputTag>("PhotonSrc");
  vertexSrc_      = pset.getParameter<edm::InputTag>("VertexSrc");
  jetSrc_         = pset.getParameter<edm::InputTag>("JetSrc");
  bJetSrc_        = pset.getParameter<edm::InputTag>("bJetSrc");
  htJetSrc_       = pset.getParameter<edm::InputTag>("htJetSrc");
  htSrc_          = pset.getParameter<edm::InputTag>("htSource");
  mhtSrc_         = pset.getParameter<edm::InputTag>("mhtSource");
  doPUReWeight_   = pset.getParameter<bool>("DoPUReweight");
  puWeightSrc_    = pset.getParameter<edm::InputTag>("PUWeightSource");
  
  //key to help getting the hlt process from event provenance
  checkedProcess_ = false;
  processName_    = "";

}

RA2ZInvPhotonTreeMaker::~RA2ZInvPhotonTreeMaker() {
  delete reducedValues;
  reducedValues = 0; 
}

void RA2ZInvPhotonTreeMaker::analyze(const edm::Event& ev, const edm::EventSetup& es) {

  using namespace edm;

  // get event-id
  unsigned int event = (ev.id()).event();
  unsigned int run   = (ev.id()).run();
  unsigned int lumi  =  ev.luminosityBlock();

  // get photons 
  edm::Handle< std::vector<pat::Photon> > patPhotons;
  ev.getByLabel(photonSrc_, patPhotons); 

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
  if(debug_ && patPhotons->size() > 0) {
    std::cout << "Isolated photons : " << std::endl;
    for( unsigned int iphot=0; iphot<patPhotons->size(); iphot++) {
      std::cout << iphot << " " <<(*patPhotons)[iphot].pt() 
		<< " " << (*patPhotons)[iphot].eta()
		<< " " << (*patPhotons)[iphot].phi() 
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
  
  m_nPhotonsIso = patPhotons->size();
  m_Photon1Pt  = (*patPhotons)[0].pt();
  m_Photon1Eta = (*patPhotons)[0].eta();
  m_Photon2Pt  = -10.;
  m_Photon2Eta = -10.;

  if (m_nPhotonsIso > 1) {
    m_Photon2Pt  = (*patPhotons)[1].pt();
    m_Photon2Eta = (*patPhotons)[1].eta();
  }

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
  m_Photon70PFMET100 = true;
  m_Photon70PFHT400 = true;
  m_Photon70PFNoPUHT400 = true;
  m_Photon135 = true;
  m_Photon150 = true;

  // Get the HLT results and check validity
  
  if (data_) {
    m_Photon70PFMET100 = false;
    m_Photon70PFHT400 = false;
    m_Photon70PFNoPUHT400 = false;
    m_Photon135 = false;
    m_Photon150 = false;

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
	if (tmpName.rfind("HLT_Photon70_CaloIdXL_PFHT400_v") != std::string::npos)
	  m_Photon70PFHT400 = true;
	else if (tmpName.rfind("HLT_Photon70_CaloIdXL_PFNoPUHT400_v") != std::string::npos)
	  m_Photon70PFNoPUHT400 = true;
	else if (tmpName.rfind("HLT_Photon70_CaloIdXL_PFMET100_v") != std::string::npos)
	  m_Photon70PFMET100 = true;
	else if (tmpName.rfind("HLT_Photon135_v") != std::string::npos)
	  m_Photon135 = true;
	else if (tmpName.rfind("HLT_Photon150_v") != std::string::npos)
	  m_Photon150 = true;
      }
    }
  }
  //if (reducedValues)
  reducedValues->Fill();
 
}

void RA2ZInvPhotonTreeMaker::beginJob() {
  //book trees
  BookTree();
}

void RA2ZInvPhotonTreeMaker::endJob() {

}


void RA2ZInvPhotonTreeMaker::BookTree() {

  edm::Service<TFileService> fs;

  reducedValues = fs->make<TTree>( "RA2Values", "Variables for reduced studies" );

  reducedValues->Branch("ra2_HT",  &m_HT,  "m_HT/D" );
  reducedValues->Branch("ra2_MHT", &m_MHT, "m_MHT/D");
  reducedValues->Branch("ra2_Vertices", &m_Vertices, "m_Vertices/I");

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

  reducedValues->Branch("ra2_nPhotonsIso", &m_nPhotonsIso, "m_nPhotonsIso/I");
  reducedValues->Branch("ra2_Photon1Pt",   &m_Photon1Pt,   "m_Photon1Pt/D" );
  reducedValues->Branch("ra2_Photon1Eta",  &m_Photon1Eta,  "m_Photon1Eta/D");
  reducedValues->Branch("ra2_Photon2Pt",   &m_Photon2Pt,   "m_Photon2Pt/D" );
  reducedValues->Branch("ra2_Photon2Eta",  &m_Photon2Eta,  "m_Photon2Eta/D");

  reducedValues->Branch("ra2_Photon70PFMET100",    &m_Photon70PFMET100,    "m_Photon70PFMET100/O" );
  reducedValues->Branch("ra2_Photon70PFHT400",     &m_Photon70PFHT400,     "m_Photon70PFHT400/O" );
  reducedValues->Branch("ra2_Photon70PFNoPUHT400", &m_Photon70PFNoPUHT400, "m_Photon70PFNoPUHT400/O" );
  reducedValues->Branch("ra2_Photon135",           &m_Photon135,           "m_Photon135/O" );
  reducedValues->Branch("ra2_Photon150",           &m_Photon150,           "m_Photon150/O" );

  reducedValues->Branch("ra2_nJetsPt30Eta50", &m_nJetsPt30Eta50, "m_nJetsPt30Eta50/I" );
  reducedValues->Branch("ra2_nJetsPt30Eta24", &m_nJetsPt30Eta24, "m_nJetsPt30Eta24/I");
  reducedValues->Branch("ra2_bJetsPt30Eta24", &m_bJetsPt30Eta24, "m_bJetsPt30Eta24/I");
  reducedValues->Branch("ra2_nJetsPt50Eta25", &m_nJetsPt50Eta25, "m_nJetsPt50Eta25/I" );

  reducedValues->SetAutoSave(1);
}


void  RA2ZInvPhotonTreeMaker::beginRun(edm::Run const& run, edm::EventSetup const& es) {
  bool changed = false;
  if (data_) {
    if (hltConfig.init(run,es,"HLT",changed)) {
      if (changed) {
	edm::LogWarning("RA2ZInvPhotonTreeMaker") << "beginRun: The HLT config has changed!";
      }
    }
    else {
      edm::LogError("TriggerEvent") << " HLT config extraction failure";
    }
  }
}

void RA2ZInvPhotonTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RA2ZInvPhotonTreeMaker);

