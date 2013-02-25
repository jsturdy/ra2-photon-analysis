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
// $Id: RA2ZInvPhotonTreeMaker.cc,v 1.14 2013/01/20 10:09:10 sturdy Exp $
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
  debugString_    = pset.getParameter< std::string >( "DebugString" );
  data_           = pset.getParameter<bool>("Data");
  scale_          = pset.getParameter<double>("ScaleFactor");
  photonSrc_      = pset.getParameter<edm::InputTag>("PhotonSrc");
  tightPhotonSrc_ = pset.getParameter<edm::InputTag>("TightPhotonSrc");
  vertexSrc_      = pset.getParameter<edm::InputTag>("VertexSrc");
  jetSrc_         = pset.getParameter<edm::InputTag>("JetSrc");
  bJetSrc_        = pset.getParameter<edm::InputTag>("bJetSrc");
  htJetSrc_       = pset.getParameter<edm::InputTag>("htJetSrc");
  htSrc_          = pset.getParameter<edm::InputTag>("htSource");
  mhtSrc_         = pset.getParameter<edm::InputTag>("mhtSource");
  ra2JetSrc_      = pset.getParameter< edm::InputTag >( "ra2JetSrc" );
  ra2HTSrc_       = pset.getParameter<edm::InputTag>("ra2HTSource");
  ra2MHTSrc_      = pset.getParameter<edm::InputTag>("ra2MHTSource");
  triggerResults_ = pset.getParameter<edm::InputTag>("TriggerResults");
  doPUReWeight_   = pset.getParameter<bool>("DoPUReweight");
  puWeightSrc_    = pset.getParameter<edm::InputTag>("PUWeightSource");
  eventWeightSrc_ = pset.getParameter< edm::InputTag >( "EventWeightSource" );

  storeExtraVetos_ = pset.getParameter<bool>("storeExtraVetos");
  ra2ElectronSrc_  = pset.getParameter<edm::InputTag>("ra2ElectronForVeto");
  ra2MuonSrc_      = pset.getParameter<edm::InputTag>("ra2MuonForVeto");
  electronVetoSrc_ = pset.getParameter<edm::InputTag>("electronVetoSource");
  muonVetoSrc_     = pset.getParameter<edm::InputTag>("muonVetoSource");
  isoTrkVetoSrc_   = pset.getParameter<edm::InputTag>("isoTrkVetoSource");

  computeMET_      = pset.getParameter<bool>("computeMET");
  if (computeMET_) {
    metSrc_         = pset.getParameter<edm::InputTag>("metSource");
    ra2METSrc_      = pset.getParameter<edm::InputTag>("ra2METSource");
  }
  genSrc_    = pset.getParameter< edm::InputTag >( "genSrc" );
  runGenStudy_     = pset.getParameter<bool>("runGenStudy");
  if (runGenStudy_) {
    maxDR_     = pset.getParameter< double >( "maxDR" );
    genJetSrc_ = pset.getParameter< edm::InputTag >( "genJetSrc" );
    if (computeMET_) 
      genMETSrc_ = pset.getParameter< edm::InputTag >( "genMETSrc" );
  }
  //key to help getting the hlt process from event provenance
  getHLTfromConfig_ = false;
  checkedProcess_ = false;
  processName_    = "";
  reducedValues = 0; 

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

  m_event = event;
  m_run   = run;
  m_lumi  = lumi;

  // get photons 
  edm::Handle< edm::View<pat::Photon> > patPhotons;
  ev.getByLabel(photonSrc_, patPhotons); 

  edm::Handle< edm::View<pat::Photon> > patPhotonsTight;
  ev.getByLabel(tightPhotonSrc_, patPhotonsTight); 

  edm::Handle< edm::View<pat::Electron> > ra2PATElectrons;
  ev.getByLabel(ra2ElectronSrc_, ra2PATElectrons); 
  m_passRA2ElVeto = true;
  if (ra2PATElectrons.isValid() && ra2PATElectrons->size()>0)
    m_passRA2ElVeto = false;

  edm::Handle< edm::View<pat::Muon> > ra2PATMuons;
  ev.getByLabel(ra2MuonSrc_, ra2PATMuons); 
  m_passRA2MuVeto = true;
  if (ra2PATMuons.isValid() && ra2PATMuons->size()>0)
    m_passRA2MuVeto = false;

  edm::Handle<reco::GenParticleCollection> gens;
  ev.getByLabel(genSrc_,gens);
  edm::Handle<reco::GenJetCollection> genJets;
  edm::Handle<edm::View<reco::GenMET> > genMET;
  if (runGenStudy_) {
    ev.getByLabel(genJetSrc_,genJets);
    //get the METs
    ev.getByLabel(genMETSrc_,genMET);
  }

  edm::Handle<reco::VertexCollection > vertices;
  ev.getByLabel(vertexSrc_, vertices);

  edm::Handle<edm::View<pat::Jet> > jets;
  ev.getByLabel(jetSrc_, jets);

  edm::Handle<edm::View<pat::Jet> > ra2Jets;
  ev.getByLabel(ra2JetSrc_, ra2Jets);

  edm::Handle<edm::View<pat::Jet> > htJets;
  ev.getByLabel(htJetSrc_, htJets);

  edm::Handle<edm::View<pat::Jet> > bJets;
  ev.getByLabel(bJetSrc_, bJets);

  edm::Handle<double > ht;
  ev.getByLabel(htSrc_, ht);

  edm::Handle<edm::View<reco::MET> > mht;
  ev.getByLabel(mhtSrc_, mht);

  edm::Handle<double > ra2HT;
  ev.getByLabel(ra2HTSrc_, ra2HT);

  edm::Handle<edm::View<reco::MET> > ra2MHT;
  ev.getByLabel(ra2MHTSrc_, ra2MHT);

  edm::Handle<edm::View<reco::MET> > met;
  edm::Handle<edm::View<reco::MET> > ra2MET;

  if (computeMET_) {
    ev.getByLabel(metSrc_, met);
    ev.getByLabel(ra2METSrc_, ra2MET);
  }
  
  if (debug_) {
    std::cout<<debugString_<<std::endl;
    std::cout<<"vertex collection has size "<<vertices->size()<<std::endl;
    std::cout<<"Full Jet collection has size "   <<ra2Jets->size()<<std::endl;
    std::cout<<"Jet collection has size "   <<jets->size()<<std::endl;
    std::cout<<"HT Jet collection has size "<<htJets->size()<<std::endl;
    std::cout<<"b-Jet collection has size " <<bJets->size()<<std::endl;
    std::cout<<"HT  value(uncorr.) "<<*ht           <<"("<<*ra2HT<<")"<<std::endl;
    std::cout<<"MHT value(uncorr.) "<<(*mht)[0].pt()<<"("<<(*ra2MHT)[0].pt()<<")"<<std::endl;
    if (computeMET_) 
      std::cout<<"MET value(uncorr.) "<<(*met)[0].pt()<<"("<<(*ra2MET)[0].pt()<<")"<<std::endl;
    
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
  if(debug_ && patPhotons->size() > 0) {
    std::cout<<debugString_<<std::endl;
    std::cout << "Isolated photons passID  passTightISO | pt  eta  phi  conv  !pixel  hadTowOverEm  (cut)  sigieie  (cut) | PU  cut  isoAlt   puSub EA" << std::endl;
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

  m_nPhotonsID         = patPhotons->size();
  m_nPhotonsTightIso   = patPhotonsTight->size();

  if (patPhotons->size()>0) {
    
    m_Photon1PDGID = 0;
    double photon1MinDRGen = 100;
    int tmpPhoton1PDGID = 0;
    m_Photon2PDGID = 0;
    double photon2MinDRGen = 100;
    int tmpPhoton2PDGID = 0;
    if (!data_) {
      reco::GenParticleCollection::const_iterator genp = gens->begin();
      for (; genp!= gens->end(); ++genp) {
	double dR = reco::deltaR(genp->eta(),genp->phi(),(*patPhotons)[0].eta(), (*patPhotons)[0].phi());
	if (dR < photon1MinDRGen) {
	  if (debug_ && dR < 1.) {
	    std::cout<<debugString_<<std::endl;
	    std::cout<<"DR("<<dR<<"), pt("<<genp->pt()<<"), eta("<<genp->eta()<<"), phi("<<genp->phi()<<")"<<std::endl;
	    std::cout<<"status("<<genp->status()<<"), pdgId("<<genp->pdgId()<<"), numberOfDaughters("<<genp->numberOfDaughters()<<")"<<std::endl;
	  }
	  photon1MinDRGen = dR;
	  tmpPhoton1PDGID = genp->pdgId();
	}
	if (m_nPhotonsID>1) {
	  dR = reco::deltaR(genp->eta(),genp->phi(),(*patPhotons)[1].eta(), (*patPhotons)[1].phi());
	  if (dR < photon2MinDRGen) {
	    if (debug_ && dR < 1.) {
	      std::cout<<debugString_<<std::endl;
	      std::cout<<"DR("<<dR<<"), pt("<<genp->pt()<<"), eta("<<genp->eta()<<"), phi("<<genp->phi()<<")"<<std::endl;
	      std::cout<<"status("<<genp->status()<<"), pdgId("<<genp->pdgId()<<"), numberOfDaughters("<<genp->numberOfDaughters()<<")"<<std::endl;
	    }
	    photon2MinDRGen = dR;
	    tmpPhoton2PDGID = genp->pdgId();
	  }
	}
      }
    }
    if (photon1MinDRGen < 0.5)
      m_Photon1PDGID = tmpPhoton1PDGID;
    if (photon2MinDRGen < 0.5)
      m_Photon2PDGID = tmpPhoton2PDGID;
    
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
    
    m_Photon1IsTightPFIso  = (((*patPhotons)[0].userFloat("pfChargedPU")<(*patPhotons)[0].userFloat("pfChargedTightCut"))&&
			      ((*patPhotons)[0].userFloat("pfNeutralPU")<(*patPhotons)[0].userFloat("pfNeutralTightCut"))&&
			      ((*patPhotons)[0].userFloat("pfGammaPU")<(*patPhotons)[0].userFloat("pfGammaTightCut"))
			      );
    
    m_Photon1MinDR  = 10.;
    m_Photon1DRJet1 = 10.;
    
    m_Photon2MinDR  = 10.;
    m_Photon2DRJet1 = 10.;
    
    m_Photon2Pt  = -10.;
    m_Photon2Eta = -10.;
    if (m_nPhotonsID > 1) {
      m_Photon2pfCH = (*patPhotons)[1].userFloat("pfChargedPU");
      m_Photon2pfNU = (*patPhotons)[1].userFloat("pfNeutralPU");
      m_Photon2pfGA = (*patPhotons)[1].userFloat("pfGammaPU");
      m_Photon2Pt  = (*patPhotons)[1].pt();
      m_Photon2Eta = (*patPhotons)[1].eta();
      m_Photon2Phi = (*patPhotons)[1].phi();
      m_Photon2SigmaIetaIeta = (*patPhotons)[1].sigmaIetaIeta();
      m_Photon2HadTowOverEm  = (*patPhotons)[1].hadTowOverEm();
      m_Photon2MinDR  = 10.;
      m_Photon2DRJet1 = 10.;
      m_Photon2EConvVeto  = (*patPhotons)[1].userFloat("passElectronConvVeto");
      m_Photon2PixelVeto  = !((*patPhotons)[1].hasPixelSeed());
    
      m_Photon2IsTightPFIso  = (((*patPhotons)[1].userFloat("pfChargedPU")<(*patPhotons)[1].userFloat("pfChargedTightCut"))&&
				((*patPhotons)[1].userFloat("pfNeutralPU")<(*patPhotons)[1].userFloat("pfNeutralTightCut"))&&
				((*patPhotons)[1].userFloat("pfGammaPU")<(*patPhotons)[1].userFloat("pfGammaTightCut"))
				);
    
    }
  }//Done with reco photon information
  
  m_nJetsPt30Eta50 = jets  ->size();
  m_nJetsPt50Eta25 = htJets->size();
  m_nJetsCSVM = bJets ->size();
  m_nJetsCSVT = 0;
  m_HT  = *ht;
  m_MHT = (*mht)[0].pt();
  m_ra2_HT  = *ra2HT;
  m_ra2_MHT = (*ra2MHT)[0].pt();

  if (computeMET_) {
    m_MET = (*met)[0].pt();
    m_ra2_MET = (*ra2MET)[0].pt();
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
    if (computeMET_) 
      m_dPhiMET1 = fabs(reco::deltaPhi(r1->phi(),(*met)[0].phi()));
    m_Jet1Pt = r1->pt();
    m_Jet1Eta = r1->eta();
    if(m_nJetsPt30Eta50 >= 2) {
      r2 = &((*jets)[1]);
      m_dPhiMHT2 = fabs(reco::deltaPhi(r2->phi(),(*mht)[0].phi())); 
      if (computeMET_) 
	m_dPhiMET2 = fabs(reco::deltaPhi(r2->phi(),(*met)[0].phi())); 
      m_Jet2Pt = r2->pt();
      m_Jet2Eta = r2->eta();
     
      if(m_nJetsPt30Eta50 >= 3) {
	r3 = &((*jets)[2]);
	m_dPhiMHT3 = fabs(reco::deltaPhi(r3->phi(),(*mht)[0].phi()));
	if (computeMET_) 
	  m_dPhiMET3 = fabs(reco::deltaPhi(r3->phi(),(*met)[0].phi()));
	m_Jet3Pt = r3->pt();
	m_Jet3Eta = r3->eta();

	if(m_nJetsPt30Eta50 >= 4) {
	  r4 = &((*jets)[3]);
	  m_dPhiMHT4 = fabs(reco::deltaPhi(r4->phi(),(*mht)[0].phi()));
	  if (computeMET_) 
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

  if (jets->size() > 0 && patPhotons->size()) {
    double dR = reco::deltaR((*patPhotons)[0].eta(),(*patPhotons)[0].phi(),(*jets)[0].eta(), (*jets)[0].phi());
    if (m_Photon1DRJet1 > dR)
      m_Photon1DRJet1 = dR;
    if (m_nPhotonsID>1) {
      dR = reco::deltaR((*patPhotons)[1].eta(),(*patPhotons)[1].phi(),(*jets)[0].eta(), (*jets)[0].phi());
      if (m_Photon2DRJet1 > dR)
	m_Photon2DRJet1 = dR;
    }
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
      if (patPhotons->size())
	if ((jet->p4()+(*patPhotons)[0].p4()).mass() > 90.0) {
	  ++m_nJetsPt50Eta25MInv;
	  m_HTMInv += jet->pt();
	}
    
    double tmpDPhi = fabs(reco::deltaPhi(jet->phi(),
					 (*mht)[0].phi()));
    if (tmpDPhi < m_dPhiMHTMin)
      m_dPhiMHTMin = tmpDPhi;
    
    if (computeMET_) {
      tmpDPhi = fabs(reco::deltaPhi(jet->phi(),(*met)[0].phi()));
      if (tmpDPhi < m_dPhiMETMin)
	m_dPhiMETMin = tmpDPhi;
    }
    
    if (patPhotons->size()) {
      double dR = reco::deltaR((*patPhotons)[0].eta(),(*patPhotons)[0].phi(),jet->eta(), jet->phi());
      if (m_Photon1MinDR > dR)
	m_Photon1MinDR = dR;
      if (m_nPhotonsID>1) {
	dR = reco::deltaR((*patPhotons)[1].eta(),(*patPhotons)[1].phi(),jet->eta(), jet->phi());
	if (m_Photon2MinDR > dR)
	  m_Photon2MinDR = dR;
      }
    }
  }///end reco photon test
  
  //saving gen information
  if (runGenStudy_ && gens->size()) {
    m_genBosons = gens->size();
    m_gen_nGenJetsPt30Eta50 = 0;
    m_gen_nGenJetsPt50Eta25 = 0;
    //loop over the gen jets, compute genHT, genMHT and nGenJets
    reco::GenJetCollection::const_iterator gjet = genJets->begin();
    double htGenJets(0.);
    reco::MET::LorentzVector mhtGen(0,0,0,0);
    for (; gjet != genJets->end(); ++gjet)
      if (gjet->pt() > 30) 
	if (fabs(gjet->eta()) < 5.0) {
	  ++m_gen_nGenJetsPt30Eta50 ;
	  mhtGen -= gjet->p4();
	  if (gjet->pt() > 50. && fabs(gjet->eta()) < 2.5) {
	    ++m_gen_nGenJetsPt50Eta25;
	    htGenJets += gjet->pt();
	  }
	}
    reco::MET genMHT = reco::MET(mhtGen, reco::MET::Point());
    m_gen_GenMHT = genMHT.pt();
    m_gen_GenHT  = htGenJets;
    
    ///loop over the reco jets and compute quantities
    //std::vector<const pat::Jet*> genRemoved_htJets;
    std::vector<const pat::Jet*> genRemoved_mhtJets;
    int jetIndex = -1;
    int genIndex  = -1;
    double bestDR = 1000.;
    int iJet(0);
    jet = ra2Jets->begin();
    for (; jet!= ra2Jets->end(); ++jet) {
      double dR = reco::deltaR((*gens)[0].eta(),(*gens)[0].phi(),jet->eta(),jet->phi());
      if (dR < bestDR) {
	bestDR = dR;
	jetIndex = iJet;
	genIndex  = 0;
      }
      ++iJet;
    }//presumably found the matched jet
    //now recompute event without that jet
    double htNoGen(0.);//, mhtGenJets(0.);
    reco::MET::LorentzVector mhtNoGen(0,0,0,0);
    m_gen_nJetsPt50Eta25 = 0;
    m_gen_nJetsPt30Eta50 = 0;
    iJet = 0;
    jet = ra2Jets->begin();
    for (; jet!= ra2Jets->end(); ++jet) {
      if ((iJet == jetIndex && bestDR > maxDR_) || iJet != jetIndex) {
	if (jet->pt() > 50 && fabs(jet->eta() < 2.5)) {
	  htNoGen += jet->pt();
	  ++m_gen_nJetsPt50Eta25;
	}
	if (jet->pt() > 30 && fabs(jet->eta() < 5.0)) {
	  mhtNoGen -= jet->p4();
	  genRemoved_mhtJets.push_back(&(*jet));
	  ++m_gen_nJetsPt30Eta50;
	}
      }
      ++iJet;
    }
    reco::MET ra2MHTNoGen = reco::MET(mhtNoGen, reco::MET::Point());
    m_gen_MHT = ra2MHTNoGen.pt();
    m_gen_HT  = htNoGen;
    
    //if (jet->pt() > 30 && fabs(jet->eta() < 5.0)) {
    //  
    //}
    m_gen_dPhiMHT1 = 10.0;  m_gen_dPhiMET1 = 10.0;
    m_gen_dPhiMHT2 = 10.0;  m_gen_dPhiMET2 = 10.0;
    m_gen_dPhiMHT3 = 10.0;  m_gen_dPhiMET3 = 10.0;
    m_gen_dPhiMHT4 = 10.0;  m_gen_dPhiMET4 = 10.0;
    
    if (genRemoved_mhtJets.size() > 0) {
      m_gen_dPhiMHT1 = fabs(reco::deltaPhi(genRemoved_mhtJets.at(0)->phi(), ra2MHTNoGen.phi()));
      //m_gen_dPhiMET1 = fabs(reco::deltaPhi(genRemoved_mhtJets.at(0)->phi(), (*genMET)[0].phi()));
      if (genRemoved_mhtJets.size() > 1) {
	m_gen_dPhiMHT2 = fabs(reco::deltaPhi(genRemoved_mhtJets.at(1)->phi(), ra2MHTNoGen.phi()));
	//m_gen_dPhiMET2 = fabs(reco::deltaPhi(genRemoved_mhtJets.at(1)->phi(), (*genMET)[0].phi()));
	if (genRemoved_mhtJets.size() > 2) {
	  m_gen_dPhiMHT3 = fabs(reco::deltaPhi(genRemoved_mhtJets.at(2)->phi(), ra2MHTNoGen.phi()));
	  //m_gen_dPhiMET3 = fabs(reco::deltaPhi(genRemoved_mhtJets.at(2)->phi(), (*genMET)[0].phi()));
	  if (genRemoved_mhtJets.size() > 3) {
	    m_gen_dPhiMHT4 = fabs(reco::deltaPhi(genRemoved_mhtJets.at(3)->phi(), ra2MHTNoGen.phi()));
	    //m_gen_dPhiMET4 = fabs(reco::deltaPhi(genRemoved_mhtJets.at(3)->phi(), (*genMET)[0].phi()));
	  }
	}
      }
    }//done looping over mht jets with GEN removed
    
    m_genBoson1Pt  = (*gens)[0].pt();
    m_genBoson1Eta = (*gens)[0].eta();
    m_genBoson1M   = (*gens)[0].mass();
    
    if (m_genBosons > 1) {
      m_genBoson2Pt  = (*gens)[1].pt();
      m_genBoson2Eta = (*gens)[1].eta();
      m_genBoson2M   = (*gens)[1].mass();
    }
    
    unsigned int pjet = 0;
    for (;pjet < genRemoved_mhtJets.size();++pjet) {
      double dR = reco::deltaR((*gens)[0].eta(),(*gens)[0].phi(),
			       genRemoved_mhtJets.at(pjet)->eta(),
			       genRemoved_mhtJets.at(pjet)->phi());
      if (m_genBoson1MinDR > dR)
	m_genBoson1MinDR = dR;
      if (gens->size()>1) {
	dR = reco::deltaR((*gens)[1].eta(),(*gens)[1].phi(),
			  genRemoved_mhtJets.at(pjet)->eta(),
			  genRemoved_mhtJets.at(pjet)->phi());
	if (m_genBoson2MinDR > dR)
	  m_genBoson2MinDR = dR;
      }
    }
    
    ///now do the matching to a reco photon
    m_genMatchRecoID        = false;
    m_genMatchRecoIDPixV    = false;
    m_genMatchRecoIDCSEV    = false;
    m_genMatchRecoIDIso     = false;
    m_gen1MatchRecoID       = false;
    m_gen1MatchRecoIDPixV   = false;
    m_gen1MatchRecoIDCSEV   = false;
    m_gen1MatchRecoIDIso    = false;
    m_reco1MatchRecoID      = false;
    m_reco1MatchRecoIDPixV  = false;
    m_reco1MatchRecoIDCSEV  = false;
    m_reco1MatchRecoIDIso   = false;
    
    if (patPhotons->size()) {
      reco::GenParticleCollection::const_iterator genp = gens->begin();
      int gphot = 0;
      if (debug_) {
	std::cout<<debugString_<<"::matching information"<<std::endl;
	printf("gen idx(pt,eta,phi) reco idx(pt,eta,phi) -- dR   pixv  csev iso\n");
      }
      for (; genp != gens->end(); ++genp) {
	int bestDRPhot = 0;
	double bestDRMin = 999.0;
	int phot = 0;
	//here or outside the genp loop?
	bool tmpPassPixV = false;
	bool tmpPassCSEV = false;
	bool tmpPassIso = false;
	edm::View<pat::Photon>::const_iterator recop = patPhotons->begin();
	for (; recop != patPhotons->end(); ++recop){
	  double dR = reco::deltaR(genp->eta(),genp->phi(),recop->eta(), recop->phi());
	  if (debug_) {
	    std::cout<<debugString_<<std::endl;
	    tmpPassPixV  = !(recop->hasPixelSeed());
	    tmpPassCSEV  = recop->userFloat("pfChargedPU");
	    tmpPassIso  = ((recop->userFloat("pfChargedPU")<recop->userFloat("pfChargedTightCut"))&&
			   (recop->userFloat("pfNeutralPU")<recop->userFloat("pfNeutralTightCut"))&&
			   (recop->userFloat("pfGammaPU")<recop->userFloat("pfGammaTightCut")));
	    printf("gen%d(%2.2f,%2.2f,%2.2f) reco%d(%2.2f,%2.2f,%2.2f) -- dR(%2.2f)   pixv(%d)  csev(%d) iso(%d)\n",
		   gphot,genp->pt(),genp->eta(),genp->phi(),
		   phot,recop->pt(),recop->eta(),recop->phi(),
		   dR,tmpPassPixV,tmpPassCSEV,tmpPassIso);
	  }
	  if (dR < bestDRMin) {
	    bestDRPhot = phot;
	    bestDRMin = dR;
	    tmpPassPixV  = !(recop->hasPixelSeed());
	    tmpPassCSEV  = recop->userFloat("pfChargedPU");
	    tmpPassIso  = ((recop->userFloat("pfChargedPU")<recop->userFloat("pfChargedTightCut"))&&
			   (recop->userFloat("pfNeutralPU")<recop->userFloat("pfNeutralTightCut"))&&
			   (recop->userFloat("pfGammaPU")<recop->userFloat("pfGammaTightCut")));
	  }
	  ++phot;
	}
	if (bestDRMin < 0.2) {
	  ////recoMatched = &((*patPhotons)[bestDRPhot]);
	  //tmpPassPixV  = !((*patPhotons)[bestDRPhot].hasPixelSeed());
	  //tmpPassCSEV  = (*patPhotons)[bestDRPhot].userFloat("pfChargedPU");
	  //tmpPassIso  = (((*patPhotons)[bestDRPhot].userFloat("pfChargedPU")<(*patPhotons)[bestDRPhot].userFloat("pfChargedTightCut"))&&
	  //		 ((*patPhotons)[bestDRPhot].userFloat("pfNeutralPU")<(*patPhotons)[bestDRPhot].userFloat("pfNeutralTightCut"))&&
	  //		 ((*patPhotons)[bestDRPhot].userFloat("pfGammaPU")<(*patPhotons)[bestDRPhot].userFloat("pfGammaTightCut")));
	  
	  m_genMatchRecoID = true;
	  
	  if (tmpPassPixV)
	    m_genMatchRecoIDPixV = true;
	  if (tmpPassCSEV)
	    m_genMatchRecoIDCSEV = true;
	  if (tmpPassIso)
	    m_genMatchRecoIDIso = true;
	  
	  if (bestDRPhot==0) {
	    m_reco1MatchRecoID = true;
	    
	    if (tmpPassPixV)
	      m_reco1MatchRecoIDPixV = true;
	    if (tmpPassCSEV)
	      m_reco1MatchRecoIDCSEV = true;
	    if (tmpPassIso)
	      m_reco1MatchRecoIDIso = true;
	  }
	  
	}
	if (m_genMatchRecoID && gphot==0) {
	  m_gen1MatchRecoID = true;
	  
	  if (tmpPassPixV)
	    m_gen1MatchRecoIDPixV = true;
	  if (tmpPassCSEV)
	    m_gen1MatchRecoIDCSEV = true;
	  if (tmpPassIso)
	    m_gen1MatchRecoIDIso = true;
	}
	++gphot;
      }//end loop over gen particles
    }
      
  }//end gen photon test
  
  //}
  
  const pat::Jet  *b1, *b2;
  m_JetCSVM1Pt  = -10.;   m_JetCSVT1Pt  = -10.;
  m_JetCSVM1Eta = -10.;   m_JetCSVT1Eta = -10.;
  m_JetCSVM2Pt  = -10.;   m_JetCSVT2Pt  = -10.;
  m_JetCSVM2Eta = -10.;   m_JetCSVT2Eta = -10.;
  edm::View<pat::Jet>::const_iterator bjet = bJets->begin();
  for (; bjet!= bJets->end(); ++bjet) {
    if (bjet->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898) 
      ++m_nJetsCSVT;
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
	ev.getByLabel(triggerResults_, hltEventHandle);
	processName_ = hltEventHandle.provenance()->processName();
	if (debug_){
	  std::cout<<debugString_<<std::endl;
	  std::cout<<processName_<<std::endl;
	}
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

  reducedValues->Branch("ra2_HT",       &m_HT,       "ra2_HT/D" );
  reducedValues->Branch("ra2_HTMInv",   &m_HTMInv,   "ra2_HTMInv/D" );
  reducedValues->Branch("ra2_MHT",      &m_MHT,      "ra2_MHT/D");
  reducedValues->Branch("ra2_ra2HT",    &m_ra2_HT,   "ra2_ra2HT/D" );
  reducedValues->Branch("ra2_ra2MHT",   &m_ra2_MHT,  "ra2_ra2MHT/D");
  reducedValues->Branch("ra2_Vertices", &m_Vertices, "ra2_Vertices/I");
  reducedValues->Branch("ra2_Event",    &m_event,    "ra2_event/I");
  reducedValues->Branch("ra2_Run",      &m_run,      "ra2_run/I");
  reducedValues->Branch("ra2_Lumi",     &m_lumi,     "ra2_lumi/I");

  reducedValues->Branch("ra2_dPhiMHT1", &m_dPhiMHT1, "ra2_dPhiMHT1/D");
  reducedValues->Branch("ra2_dPhiMHT2", &m_dPhiMHT2, "ra2_dPhiMHT2/D");
  reducedValues->Branch("ra2_dPhiMHT3", &m_dPhiMHT3, "ra2_dPhiMHT3/D");
  reducedValues->Branch("ra2_dPhiMHT4", &m_dPhiMHT4, "ra2_dPhiMHT4/D");
  reducedValues->Branch("ra2_dPhiMHTMin", &m_dPhiMHTMin, "ra2_dPhiMHTMin/D");

  if (computeMET_) {
    reducedValues->Branch("ra2_MET",      &m_MET,      "ra2_MET/D");
    reducedValues->Branch("ra2_ra2MET",   &m_ra2_MET,  "ra2_ra2MET/D");
  
    reducedValues->Branch("ra2_dPhiMET1", &m_dPhiMET1, "ra2_dPhiMET1/D");
    reducedValues->Branch("ra2_dPhiMET2", &m_dPhiMET2, "ra2_dPhiMET2/D");
    reducedValues->Branch("ra2_dPhiMET3", &m_dPhiMET3, "ra2_dPhiMET3/D");
    reducedValues->Branch("ra2_dPhiMET4", &m_dPhiMET4, "ra2_dPhiMET4/D");
    reducedValues->Branch("ra2_dPhiMETMin", &m_dPhiMETMin, "ra2_dPhiMETMin/D");
  }
  if (runGenStudy_) {
    reducedValues->Branch("ra2_gen_HT",    &m_gen_HT,    "ra2_gen_HT/D" );
    reducedValues->Branch("ra2_gen_MHT",   &m_gen_MHT,   "ra2_gen_MHT/D");
    reducedValues->Branch("ra2_gen_GenHT", &m_gen_GenHT, "ra2_gen_GenHT/D" );
    reducedValues->Branch("ra2_gen_GenMHT",&m_gen_GenMHT,"ra2_gen_GenMHT/D");
    
    reducedValues->Branch("ra2_genBosons",     &m_genBosons,     "ra2_genBosons/I" );
    reducedValues->Branch("ra2_genBoson1Pt",   &m_genBoson1Pt,   "ra2_genBoson1Pt/D" );
    reducedValues->Branch("ra2_genBoson1Eta",  &m_genBoson1Eta,  "ra2_genBoson1Eta/D" );
    reducedValues->Branch("ra2_genBoson1M",    &m_genBoson1M,    "ra2_genBoson1M/D" );
    reducedValues->Branch("ra2_genBoson1MinDR",&m_genBoson1MinDR,"ra2_genBoson1MinDR/D" );

    reducedValues->Branch("ra2_gen_dPhiMHT1",   &m_gen_dPhiMHT1,   "ra2_gen_dPhiMHT1/D");
    reducedValues->Branch("ra2_gen_dPhiMHT2",   &m_gen_dPhiMHT2,   "ra2_gen_dPhiMHT2/D");
    reducedValues->Branch("ra2_gen_dPhiMHT3",   &m_gen_dPhiMHT3,   "ra2_gen_dPhiMHT3/D");
    reducedValues->Branch("ra2_gen_dPhiMHT4",   &m_gen_dPhiMHT4,   "ra2_gen_dPhiMHT4/D");
    
    if (computeMET_) {
      reducedValues->Branch("ra2_gen_MET",      &m_gen_MET,      "ra2_gen_MET/D");
      reducedValues->Branch("ra2_gen_GenMET",   &m_gen_GenMET,   "ra2_gen_GenMET/D");
      reducedValues->Branch("ra2_gen_dPhiMET1", &m_gen_dPhiMET1, "ra2_gen_dPhiMET1/D");
      reducedValues->Branch("ra2_gen_dPhiMET2", &m_gen_dPhiMET2, "ra2_gen_dPhiMET2/D");
      reducedValues->Branch("ra2_gen_dPhiMET3", &m_gen_dPhiMET3, "ra2_gen_dPhiMET3/D");
      reducedValues->Branch("ra2_gen_dPhiMET4", &m_gen_dPhiMET4, "ra2_gen_dPhiMET4/D");
    }

    reducedValues->Branch("ra2_gen_nJetsPt30Eta50",   &m_gen_nJetsPt30Eta50, "ra2_gen_nJetsPt30Eta50/I" );
    reducedValues->Branch("ra2_gen_nJetsPt50Eta25",   &m_gen_nJetsPt50Eta25, "ra2_gen_nJetsPt50Eta25/I" );
    reducedValues->Branch("ra2_gen_nGenJetsPt30Eta50",&m_gen_nGenJetsPt30Eta50, "ra2_gen_nGenJetsPt30Eta50/I" );
    reducedValues->Branch("ra2_gen_nGenJetsPt50Eta25",&m_gen_nGenJetsPt50Eta25, "ra2_gen_nGenJetsPt50Eta25/I" );

    reducedValues->Branch("ra2_gen_genMatchRecoID",      &m_genMatchRecoID,      "ra2_gen_genMatchRecoID/O");
    reducedValues->Branch("ra2_gen_genMatchRecoIDPixV",  &m_genMatchRecoIDPixV,  "ra2_gen_genMatchRecoIDPixV/O");
    reducedValues->Branch("ra2_gen_genMatchRecoIDCSEV",  &m_genMatchRecoIDCSEV,  "ra2_gen_genMatchRecoIDCSEV/O");
    reducedValues->Branch("ra2_gen_genMatchRecoIDIso",   &m_genMatchRecoIDIso,   "ra2_gen_genMatchRecoIDIso/O");
    reducedValues->Branch("ra2_gen_reco1MatchRecoID",    &m_reco1MatchRecoID,    "ra2_gen_reco1MatchRecoID/O");
    reducedValues->Branch("ra2_gen_reco1MatchRecoIDPixV",&m_reco1MatchRecoIDPixV,"ra2_gen_reco1MatchRecoIDPixV/O");
    reducedValues->Branch("ra2_gen_reco1MatchRecoIDCSEV",&m_reco1MatchRecoIDCSEV,"ra2_gen_reco1MatchRecoIDCSEV/O");
    reducedValues->Branch("ra2_gen_reco1MatchRecoIDIso", &m_reco1MatchRecoIDIso, "ra2_gen_reco1MatchRecoIDIso/O");
    reducedValues->Branch("ra2_gen_gen1MatchRecoID",     &m_gen1MatchRecoID,     "ra2_gen_gen1MatchRecoID/O");
    reducedValues->Branch("ra2_gen_gen1MatchRecoIDPixV", &m_gen1MatchRecoIDPixV, "ra2_gen_gen1MatchRecoIDPixV/O");
    reducedValues->Branch("ra2_gen_gen1MatchRecoIDCSEV", &m_gen1MatchRecoIDCSEV, "ra2_gen_gen1MatchRecoIDCSEV/O");
    reducedValues->Branch("ra2_gen_gen1MatchRecoIDIso",  &m_gen1MatchRecoIDIso,  "ra2_gen_gen1MatchRecoIDIso/O");
  }

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

  reducedValues->Branch("ra2_nPhotonsID",        &m_nPhotonsID,        "ra2_nPhotonsID/I");
  reducedValues->Branch("ra2_nPhotonsTightIso",  &m_nPhotonsTightIso,  "ra2_nPhotonsTightIso/I");

  reducedValues->Branch("ra2_Photon1PDGID",&m_Photon1PDGID,"ra2_Photon1PDGID/I");
  reducedValues->Branch("ra2_Photon1pfCH", &m_Photon1pfCH, "ra2_Photon1pfCH/D");
  reducedValues->Branch("ra2_Photon1pfNU", &m_Photon1pfNU, "ra2_Photon1pfNU/D");
  reducedValues->Branch("ra2_Photon1pfGA", &m_Photon1pfGA, "ra2_Photon1pfGA/D");
  reducedValues->Branch("ra2_Photon1Pt",   &m_Photon1Pt,   "ra2_Photon1Pt/D" );
  reducedValues->Branch("ra2_Photon1Eta",  &m_Photon1Eta,  "ra2_Photon1Eta/D");
  reducedValues->Branch("ra2_Photon1Phi",  &m_Photon1Phi,  "ra2_Photon1Phi/D");
  reducedValues->Branch("ra2_Photon1SigmaIetaIeta", &m_Photon1SigmaIetaIeta, "ra2_Photon1SigmaIetaIeta/D");
  reducedValues->Branch("ra2_Photon1HadTowOverEm",  &m_Photon1HadTowOverEm,  "ra2_Photon1HadTowOverEm/D");
  reducedValues->Branch("ra2_Photon1MinDR", &m_Photon1MinDR, "ra2_Photon1MinDR/D");
  reducedValues->Branch("ra2_Photon1DRJet1",&m_Photon1DRJet1,"ra2_Photon1DRJet1/D");
  reducedValues->Branch("ra2_Photon1EConvVeto",   &m_Photon1EConvVeto,   "ra2_Photon1EConvVeto/O" );
  reducedValues->Branch("ra2_Photon1PixelVeto",   &m_Photon1PixelVeto,   "ra2_Photon1PixelVeto/O" );
  reducedValues->Branch("ra2_Photon1IsTightPFIso",&m_Photon1IsTightPFIso,"ra2_Photon1IsTightPFIso/O" );

  reducedValues->Branch("ra2_Photon2PDGID",&m_Photon2PDGID,"ra2_Photon2PDGID/I");
  reducedValues->Branch("ra2_Photon2pfCH", &m_Photon2pfCH, "ra2_Photon2pfCH/D");
  reducedValues->Branch("ra2_Photon2pfNU", &m_Photon2pfNU, "ra2_Photon2pfNU/D");
  reducedValues->Branch("ra2_Photon2pfGA", &m_Photon2pfGA, "ra2_Photon2pfGA/D");
  reducedValues->Branch("ra2_Photon2Pt",   &m_Photon2Pt,   "ra2_Photon2Pt/D" );
  reducedValues->Branch("ra2_Photon2Eta",  &m_Photon2Eta,  "ra2_Photon2Eta/D");
  reducedValues->Branch("ra2_Photon2Phi",  &m_Photon2Phi,  "ra2_Photon2Phi/D");
  reducedValues->Branch("ra2_Photon2SigmaIetaIeta", &m_Photon2SigmaIetaIeta, "ra2_Photon2SigmaIetaIeta/D");
  reducedValues->Branch("ra2_Photon2HadTowOverEm",  &m_Photon2HadTowOverEm,  "ra2_Photon2HadTowOverEm/D");
  reducedValues->Branch("ra2_Photon2MinDR", &m_Photon2MinDR, "ra2_Photon2MinDR/D");
  reducedValues->Branch("ra2_Photon2DRJet1",&m_Photon2DRJet1,"ra2_Photon2DRJet1/D");
  reducedValues->Branch("ra2_Photon2EConvVeto",   &m_Photon2EConvVeto,   "ra2_Photon2EConvVeto/O" );
  reducedValues->Branch("ra2_Photon2PixelVeto",   &m_Photon2PixelVeto,   "ra2_Photon2PixelVeto/O" );
  reducedValues->Branch("ra2_Photon2IsTightPFIso",&m_Photon2IsTightPFIso,"ra2_Photon2IsTightPFIso/O" );

  reducedValues->Branch("ra2_Photon70PFMET100",    &m_Photon70PFMET100,    "ra2_Photon70PFMET100/O" );
  reducedValues->Branch("ra2_Photon70PFHT400",     &m_Photon70PFHT400,     "ra2_Photon70PFHT400/O" );
  reducedValues->Branch("ra2_Photon70PFNoPUHT400", &m_Photon70PFNoPUHT400, "ra2_Photon70PFNoPUHT400/O" );
  reducedValues->Branch("ra2_Photon135",           &m_Photon135,           "ra2_Photon135/O" );
  reducedValues->Branch("ra2_Photon150",           &m_Photon150,           "ra2_Photon150/O" );

  reducedValues->Branch("ra2_nJetsCSVM", &m_nJetsCSVM, "ra2_nJetsCSVM/I");
  reducedValues->Branch("ra2_nJetsCSVT", &m_nJetsCSVT, "ra2_nJetsCSVT/I");
  reducedValues->Branch("ra2_nJetsPt30Eta50", &m_nJetsPt30Eta50, "ra2_nJetsPt30Eta50/I" );
  reducedValues->Branch("ra2_nJetsPt30Eta24", &m_nJetsPt30Eta24, "ra2_nJetsPt30Eta24/I");
  reducedValues->Branch("ra2_nJetsPt50Eta24", &m_nJetsPt50Eta24, "ra2_nJetsPt50Eta24/I");
  reducedValues->Branch("ra2_nJetsPt70Eta24", &m_nJetsPt70Eta24, "ra2_nJetsPt70Eta24/I");
  reducedValues->Branch("ra2_nJetsPt50Eta25", &m_nJetsPt50Eta25, "ra2_nJetsPt50Eta25/I" );
  reducedValues->Branch("ra2_nJetsPt50Eta25MInv", &m_nJetsPt50Eta25MInv, "ra2_nJetsPt50Eta25MInv/I" );


  if (storeExtraVetos_) {
    reducedValues->Branch("ra2_passRA2ElVeto",    &m_passRA2ElVeto   , "ra2_passRA2ElVeto/O"    );
    reducedValues->Branch("ra2_passRA2MuVeto",    &m_passRA2MuVeto   , "ra2_passRA2MuVeto/O"    );
    reducedValues->Branch("ra2_passDirIsoElVeto", &m_passDirIsoElVeto, "ra2_passDirIsoElVeto/O"    );
    reducedValues->Branch("ra2_passDirIsoMuVeto", &m_passDirIsoMuVeto, "ra2_passDirIsoMuVeto/O"    );
    reducedValues->Branch("ra2_passIsoTrkVeto",   &m_passIsoTrkVeto  , "ra2_passIsoTrkVeto/O");
  }
  
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


//  LocalWords:  reco
