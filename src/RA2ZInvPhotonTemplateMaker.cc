// -*- C++ -*-
//
// Package:    RA2ZInvPhotonTemplateMaker
// Class:      RA2ZInvPhotonTemplateMaker
// 
/**\class RA2ZInvPhotonTemplateMaker RA2ZInvPhotonTemplateMaker.cc SusyAnalysis/RA2ZInvPhotonTemplateMaker/src/RA2ZInvPhotonTemplateMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seema Sharma
//         Created:  Mon Jun 20 12:58:08 CDT 2011
// $Id: RA2ZInvPhotonTemplateMaker.cc,v 1.14 2013/01/20 10:09:10 sturdy Exp $
//
//


// system include files
#include <memory>
#include <iomanip>
#include <cmath>

#include "ZInvisibleBkgds/Photons/interface/RA2ZInvPhotonTemplateMaker.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include <DataFormats/METReco/interface/MET.h>


RA2ZInvPhotonTemplateMaker::RA2ZInvPhotonTemplateMaker(const edm::ParameterSet& pset) {

  // read parameters from config file
  debug_          = pset.getParameter<bool>("Debug");
  debugString_    = pset.getParameter< std::string >( "DebugString" );
  data_           = pset.getParameter<bool>("Data");
  scale_          = pset.getParameter<double>("ScaleFactor");
  photonSrc_      = pset.getParameter<edm::InputTag>("PhotonSrc");
  tightPhotonSrc_ = pset.getParameter<edm::InputTag>("TightPhotonSrc");
  vertexSrc_      = pset.getParameter<edm::InputTag>("VertexSrc");
  jetSrc_         = pset.getParameter<edm::InputTag>("JetSrc");
  htJetSrc_       = pset.getParameter<edm::InputTag>("htJetSrc");
  htSrc_          = pset.getParameter<edm::InputTag>("htSource");
  mhtSrc_         = pset.getParameter<edm::InputTag>("mhtSource");
  metSrc_         = pset.getParameter<edm::InputTag>("metSource");
  doPUReWeight_   = pset.getParameter<bool>("DoPUReweight");
  puWeightSrc_    = pset.getParameter<edm::InputTag>("PUWeightSource");
  eventWeightSrc_ = pset.getParameter< edm::InputTag >( "EventWeightSource" );

  reducedValues = 0; 

}

RA2ZInvPhotonTemplateMaker::~RA2ZInvPhotonTemplateMaker() {
  delete reducedValues;
  reducedValues = 0; 
}

void RA2ZInvPhotonTemplateMaker::analyze(const edm::Event& ev, const edm::EventSetup& es) {

  using namespace edm;

  // get event-id
  unsigned int event = (ev.id()).event();
  unsigned int run   = (ev.id()).run();
  unsigned int lumi  =  ev.luminosityBlock();

  m_event = event;
  m_run   = run;
  m_lumi  = lumi;

  // get photons 
  edm::Handle< std::vector<pat::Photon> > patPhotons;
  ev.getByLabel(photonSrc_, patPhotons); 

  edm::Handle< std::vector<pat::Photon> > patPhotonsTight;
  ev.getByLabel(tightPhotonSrc_, patPhotonsTight); 

  edm::Handle<reco::GenParticleCollection> gens;
  if (!data_)
    ev.getByLabel("genParticles",gens);

  edm::Handle<reco::VertexCollection > vertices;
  ev.getByLabel(vertexSrc_, vertices);

  edm::Handle<edm::View<pat::Jet> > jets;
  ev.getByLabel(jetSrc_, jets);

  edm::Handle<edm::View<pat::Jet> > htJets;
  ev.getByLabel(htJetSrc_, htJets);

  edm::Handle<double > ht;
  ev.getByLabel(htSrc_, ht);

  edm::Handle<edm::View<reco::MET> > mht;
  ev.getByLabel(mhtSrc_, mht);

  edm::Handle<edm::View<reco::MET> > met;
  ev.getByLabel(metSrc_, met);


  if (debug_) {
    std::cout<<"vertex collection has size "<<vertices->size()<<std::endl;
    std::cout<<"Jet collection has size "   <<jets->size()<<std::endl;
    std::cout<<"HT Jet collection has size "<<htJets->size()<<std::endl;
    std::cout<<"HT value " <<*ht<<std::endl;
    std::cout<<"MHT value "<<(*mht)[0].pt()<<std::endl;
    std::cout<<"MET value "<<(*met)[0].pt()<<std::endl;
   
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
  
  //////
  if(debug_ && patPhotons->size() > 0) {
    std::cout<<debugString_<<std::endl;
    std::cout << "Isolated photons passTightID  passTightISO | pt  eta  phi  conv  !pixel  hadTowOverEm  (cut)  sigieie  (cut) | PU  cut  isoAlt   puSub EA" << std::endl;
    for( unsigned int iphot=0; iphot<patPhotons->size(); iphot++) {
      bool passTightID = false;
      bool passTightISO = false;

      if ((*patPhotons)[iphot].et() > 100 &&
	  (fabs((*patPhotons)[iphot].eta()) < 2.5 &&
	   (
	    fabs((*patPhotons)[iphot].eta()) < 1.4442 ||
	    fabs((*patPhotons)[iphot].eta()) > 1.566 
	    ))&&
	  ((*patPhotons)[iphot].hadTowOverEm() < (*patPhotons)[iphot].userFloat("hadTowOverEmTightCut")) &&
	  ((*patPhotons)[iphot].sigmaIetaIeta() < (*patPhotons)[iphot].userFloat("showerShapeTightCut"))
	  )
	passTightID = true;
      
      if (((*patPhotons)[iphot].userFloat("pfChargedPU") < (*patPhotons)[iphot].userFloat("pfChargedTightCut")) &&
	  ((*patPhotons)[iphot].userFloat("pfNeutralPU") < (*patPhotons)[iphot].userFloat("pfNeutralTightCut")) &&
	  ((*patPhotons)[iphot].userFloat("pfGammaPU")   < (*patPhotons)[iphot].userFloat("pfGammaTightCut"))
	  )
	passTightISO = true;

      std::cout << "ph" << iphot 
		<< " " << passTightID
		<< " " << passTightISO<<std::endl
		<< " | " <<(*patPhotons)[iphot].pt() 
		<< " " << (*patPhotons)[iphot].eta()
		<< " " << (*patPhotons)[iphot].phi();
      if (!data_ && (*patPhotons)[iphot].genPhoton())
	std::cout<< " " << (*patPhotons)[iphot].genPhoton()->pdgId();
      std::cout<< " " << (*patPhotons)[iphot].userFloat("passElectronConvVeto") 
		<< " " << !((*patPhotons)[iphot].hasPixelSeed() )
		<< " " << (*patPhotons)[iphot].hadTowOverEm() 
		<< " " << (*patPhotons)[iphot].userFloat("hadTowOverEmTightCut") 
		<< " " << (*patPhotons)[iphot].sigmaIetaIeta() 
		<< " " << (*patPhotons)[iphot].userFloat("showerShapeTightCut") 
		<< std::endl
		<< " | pfCharged - " << (*patPhotons)[iphot].userFloat("pfChargedPU") 
		<< " " << (*patPhotons)[iphot].userFloat("pfChargedTightCut") 
		<< " " << (*patPhotons)[iphot].userFloat("pfChargedIsoAlt") 
		<< " " << (*patPhotons)[iphot].userFloat("pfChargedPUSub") 
		<< " " << (*patPhotons)[iphot].userFloat("pfChargedEA") 
		<< std::endl
		<< " | pfNeutral - " << (*patPhotons)[iphot].userFloat("pfNeutralPU") 
		<< " " << (*patPhotons)[iphot].userFloat("pfNeutralTightCut") 
		<< " " << (*patPhotons)[iphot].userFloat("pfNeutralIsoAlt") 
		<< " " << (*patPhotons)[iphot].userFloat("pfNeutralPUSub") 
		<< " " << (*patPhotons)[iphot].userFloat("pfNeutralEA")
		<< std::endl
		<< " | pfGamma - " << (*patPhotons)[iphot].userFloat("pfGammaPU") 
		<< " " << (*patPhotons)[iphot].userFloat("pfGammaTightCut") 
		<< " " << (*patPhotons)[iphot].userFloat("pfGammaIsoAlt") 
		<< " " << (*patPhotons)[iphot].userFloat("pfGammaPUSub") 
		<< " " << (*patPhotons)[iphot].userFloat("pfGammaEA") 
		<< std::endl;
    }
  }
  //////

  if (patPhotons->size() < 1)
    return;

  m_nPhotonsID         = patPhotons->size();
  m_nPhotonsTightIso   = patPhotonsTight->size();

  m_Photon1PDGID = 0;
  double photon1MinDRGen = 100;
  int tmpPhoton1PDGID = 0;
  if (!data_) {
    reco::GenParticleCollection::const_iterator genp = gens->begin();
    for (; genp!= gens->end(); ++genp) {
      double dR = reco::deltaR(genp->eta(),genp->phi(),(*patPhotons)[0].eta(), (*patPhotons)[0].phi());
      if (dR < photon1MinDRGen) {
	if (debug_ && dR < 1.) {
	  std::cout<<"DR("<<dR<<"), pt("<<genp->pt()<<"), eta("<<genp->eta()<<"), phi("<<genp->phi()<<")"<<std::endl;
	  std::cout<<"status("<<genp->status()<<"), pdgId("<<genp->pdgId()<<"), numberOfDaughters("<<genp->numberOfDaughters()<<")"<<std::endl;
	}
	photon1MinDRGen = dR;
	tmpPhoton1PDGID = genp->pdgId();
      }
    }
  }
  if (photon1MinDRGen < 0.5)
    m_Photon1PDGID = tmpPhoton1PDGID;
  
  m_Photon1pfCH  = (*patPhotons)[0].userFloat("pfChargedPU");
  m_Photon1pfNU  = (*patPhotons)[0].userFloat("pfNeutralPU");
  m_Photon1pfGA  = (*patPhotons)[0].userFloat("pfGammaPU");
  m_Photon1Pt  = (*patPhotons)[0].pt();
  m_Photon1Eta = (*patPhotons)[0].eta();
  m_Photon1Phi = (*patPhotons)[0].phi();
  m_Photon1SigmaIetaIeta = (*patPhotons)[0].sigmaIetaIeta();
  m_Photon1HadTowOverEm  = (*patPhotons)[0].hadTowOverEm();
  m_Photon1EConvVeto  = (*patPhotons)[0].userFloat("passElectronConvVeto");
  m_Photon1PixelVeto  = !((*patPhotons)[0].hasPixelSeed());

  m_Photon1IsTightID  = (((*patPhotons)[0].hadTowOverEm() < (*patPhotons)[0].userFloat("hadTowOverEmTightCut")) &&
			 ((*patPhotons)[0].sigmaIetaIeta() < (*patPhotons)[0].userFloat("showerShapeTightCut"))
			 );

  m_Photon1IsTightPFIso  = (((*patPhotons)[0].userFloat("pfChargedPU")<(*patPhotons)[0].userFloat("pfChargedTightCut"))&&
			    ((*patPhotons)[0].userFloat("pfNeutralPU")<(*patPhotons)[0].userFloat("pfNeutralTightCut"))&&
			    ((*patPhotons)[0].userFloat("pfGammaPU")<(*patPhotons)[0].userFloat("pfGammaTightCut"))
			    );
  m_Photon1PassPFCh  = ((*patPhotons)[0].userFloat("pfChargedPU")<(*patPhotons)[0].userFloat("pfChargedTightCut"));
  m_Photon1PassPFNu  = ((*patPhotons)[0].userFloat("pfNeutralPU")<(*patPhotons)[0].userFloat("pfNeutralTightCut"));
  m_Photon1PassPFGa  = ((*patPhotons)[0].userFloat("pfGammaPU")<(*patPhotons)[0].userFloat("pfGammaTightCut"));
  
  
  m_Photon1MinDR  = 10.;
  m_Photon1DRJet1 = 10.;

  m_nJetsPt30Eta50 = jets  ->size();
  m_nJetsPt50Eta25 = htJets->size();
  m_HT  = *ht;
  m_MHT = (*mht)[0].pt();
  m_MET = (*met)[0].pt();

  const pat::Jet  *r1, *r2, *r3, *r4;
  m_dPhiMHT1 = 10.0;  m_dPhiMET1 = 10.0;
  m_dPhiMHT2 = 10.0;  m_dPhiMET2 = 10.0;
  m_dPhiMHT3 = 10.0;  m_dPhiMET3 = 10.0;
  m_dPhiMHT4 = 10.0;  m_dPhiMET4 = 10.0;

  if(m_nJetsPt30Eta50 >= 1) {
    r1 = &((*jets)[0]);
    m_dPhiMHT1 = fabs(reco::deltaPhi(r1->phi(),(*mht)[0].phi()));
    m_dPhiMET1 = fabs(reco::deltaPhi(r1->phi(),(*met)[0].phi()));
    if(m_nJetsPt30Eta50 >= 2) {
      r2 = &((*jets)[1]);
      m_dPhiMHT2 = fabs(reco::deltaPhi(r2->phi(),(*mht)[0].phi())); 
      m_dPhiMET2 = fabs(reco::deltaPhi(r2->phi(),(*met)[0].phi())); 
     
      if(m_nJetsPt30Eta50 >= 3) {
	r3 = &((*jets)[2]);
	m_dPhiMHT3 = fabs(reco::deltaPhi(r3->phi(),(*mht)[0].phi()));
	m_dPhiMET3 = fabs(reco::deltaPhi(r3->phi(),(*met)[0].phi()));

	if(m_nJetsPt30Eta50 >= 4) {
	  r4 = &((*jets)[3]);
	  m_dPhiMHT4 = fabs(reco::deltaPhi(r4->phi(),(*mht)[0].phi()));
	  m_dPhiMET4 = fabs(reco::deltaPhi(r4->phi(),(*met)[0].phi()));
	}
      }
    }
  }

  m_dPhiMHTMin  = 10.;
  m_dPhiMETMin  = 10.;

  if (jets->size() > 0) {
    double dR = reco::deltaR((*patPhotons)[0].eta(),(*patPhotons)[0].phi(),(*jets)[0].eta(), (*jets)[0].phi());
    if (m_Photon1DRJet1 > dR)
      m_Photon1DRJet1 = dR;
  }

  edm::View<pat::Jet>::const_iterator jet = jets->begin();
  for (; jet!= jets->end(); ++jet) {
    double dR = reco::deltaR((*patPhotons)[0].eta(),(*patPhotons)[0].phi(),jet->eta(), jet->phi());
    if (m_Photon1MinDR > dR)
      m_Photon1MinDR = dR;
    
    double tmpDPhi = fabs(reco::deltaPhi(jet->phi(),(*mht)[0].phi()));
    if (tmpDPhi < m_dPhiMHTMin)
      m_dPhiMHTMin = tmpDPhi;
    
    tmpDPhi = fabs(reco::deltaPhi(jet->phi(),(*met)[0].phi()));
    if (tmpDPhi < m_dPhiMETMin)
      m_dPhiMETMin = tmpDPhi;
  }


  //if (reducedValues)
  reducedValues->Fill();
 
}

void RA2ZInvPhotonTemplateMaker::beginJob() {
  //book trees
  BookTree();
}

void RA2ZInvPhotonTemplateMaker::endJob() {

}


void RA2ZInvPhotonTemplateMaker::BookTree() {

  edm::Service<TFileService> fs;

  reducedValues = fs->make<TTree>( "RA2Values", "Variables for reduced studies" );

  reducedValues->Branch("ra2_HT",       &m_HT,       "m_HT/D" );
  reducedValues->Branch("ra2_MHT",      &m_MHT,      "m_MHT/D");
  reducedValues->Branch("ra2_MET",      &m_MET,      "m_MET/D");
  reducedValues->Branch("ra2_Vertices", &m_Vertices, "m_Vertices/I");
  reducedValues->Branch("ra2_Event",    &m_event,    "m_event/I");
  reducedValues->Branch("ra2_Run",      &m_run,      "m_run/I");
  reducedValues->Branch("ra2_Lumi",     &m_lumi,     "m_lumi/I");

  
  reducedValues->Branch("ra2_dPhiMHT1", &m_dPhiMHT1, "m_dPhiMHT1/D");
  reducedValues->Branch("ra2_dPhiMHT2", &m_dPhiMHT2, "m_dPhiMHT2/D");
  reducedValues->Branch("ra2_dPhiMHT3", &m_dPhiMHT3, "m_dPhiMHT3/D");
  reducedValues->Branch("ra2_dPhiMHT4", &m_dPhiMHT4, "m_dPhiMHT4/D");
  reducedValues->Branch("ra2_dPhiMHTMin", &m_dPhiMHTMin, "m_dPhiMHTMin/D");

  reducedValues->Branch("ra2_dPhiMET1", &m_dPhiMET1, "m_dPhiMET1/D");
  reducedValues->Branch("ra2_dPhiMET2", &m_dPhiMET2, "m_dPhiMET2/D");
  reducedValues->Branch("ra2_dPhiMET3", &m_dPhiMET3, "m_dPhiMET3/D");
  reducedValues->Branch("ra2_dPhiMET4", &m_dPhiMET4, "m_dPhiMET4/D");
  reducedValues->Branch("ra2_dPhiMETMin", &m_dPhiMETMin, "m_dPhiMETMin/D");

  reducedValues->Branch("ra2_PUWt",    &m_PUWt,    "m_PUWt/D");
  reducedValues->Branch("ra2_EventWt", &m_EventWt, "m_EventWt/D");

  reducedValues->Branch("ra2_nPhotonsID",        &m_nPhotonsID,        "m_nPhotonsID/I");
  reducedValues->Branch("ra2_nPhotonsTightIso",  &m_nPhotonsTightIso,  "m_nPhotonsTightIso/I");

  reducedValues->Branch("ra2_Photon1PDGID",&m_Photon1PDGID,"m_Photon1PDGID/I");

  reducedValues->Branch("ra2_Photon1pfCH", &m_Photon1pfCH, "m_Photon1pfCH/D");
  reducedValues->Branch("ra2_Photon1pfNU", &m_Photon1pfNU, "m_Photon1pfNU/D");
  reducedValues->Branch("ra2_Photon1pfGA", &m_Photon1pfGA, "m_Photon1pfGA/D");
  reducedValues->Branch("ra2_Photon1Pt",   &m_Photon1Pt,   "m_Photon1Pt/D" );
  reducedValues->Branch("ra2_Photon1Eta",  &m_Photon1Eta,  "m_Photon1Eta/D");
  reducedValues->Branch("ra2_Photon1Phi",  &m_Photon1Phi,  "m_Photon1Phi/D");
  reducedValues->Branch("ra2_Photon1SigmaIetaIeta", &m_Photon1SigmaIetaIeta, "m_Photon1SigmaIetaIeta/D");
  reducedValues->Branch("ra2_Photon1HadTowOverEm",  &m_Photon1HadTowOverEm,  "m_Photon1HadTowOverEm/D");
  reducedValues->Branch("ra2_Photon1MinDR", &m_Photon1MinDR, "m_Photon1MinDR/D");
  reducedValues->Branch("ra2_Photon1DRJet1",&m_Photon1DRJet1,"m_Photon1DRJet1/D");
  reducedValues->Branch("ra2_Photon1EConvVeto",   &m_Photon1EConvVeto,   "m_Photon1EConvVeto/O" );
  reducedValues->Branch("ra2_Photon1PixelVeto",   &m_Photon1PixelVeto,   "m_Photon1PixelVeto/O" );
  reducedValues->Branch("ra2_Photon1IsTightID",   &m_Photon1IsTightID,   "m_Photon1IsTightID/O" );
  reducedValues->Branch("ra2_Photon1IsTightPFIso",   &m_Photon1IsTightPFIso,   "m_Photon1IsTightPFIso/O" );
  reducedValues->Branch("ra2_Photon1PassPFCh",   &m_Photon1PassPFCh,   "m_Photon1PassPFCh/O" );
  reducedValues->Branch("ra2_Photon1PassPFNu",   &m_Photon1PassPFNu,   "m_Photon1PassPFNu/O" );
  reducedValues->Branch("ra2_Photon1PassPFGa",   &m_Photon1PassPFGa,   "m_Photon1PassPFGa/O" );

  reducedValues->Branch("ra2_nJetsPt30Eta50", &m_nJetsPt30Eta50, "m_nJetsPt30Eta50/I" );
  reducedValues->Branch("ra2_nJetsPt50Eta25", &m_nJetsPt50Eta25, "m_nJetsPt50Eta25/I" );

  reducedValues->SetAutoSave(1);
}


void  RA2ZInvPhotonTemplateMaker::beginRun(edm::Run const& run, edm::EventSetup const& es) {
}

void RA2ZInvPhotonTemplateMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonTemplateMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonTemplateMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonTemplateMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RA2ZInvPhotonTemplateMaker);


//  LocalWords:  reco
