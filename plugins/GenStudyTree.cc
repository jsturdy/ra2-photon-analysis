// -*- C++ -*-
//
// Package:    Photons
// Class:      Photons
// 
/**\class Photons GenStudyTree.cc ZInvisibleBkgds/Photons/plugins/GenStudyTree.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: GenStudyTree.cc,v 1.2 2012/09/04 10:26:08 sturdy Exp $
//
//
// system include files
#include <cmath>

#include "ZInvisibleBkgds/Photons/interface/GenStudyTree.h"

#include "DataFormats/Math/interface/deltaPhi.h"
//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
GenStudyTree::GenStudyTree(const edm::ParameterSet& pset) :
  debug_        ( pset.getParameter< bool >( "debug" ) ),
  scale_        ( pset.getParameter<double>("ScaleFactor") ),
  debugString_  ( pset.getParameter< std::string >( "debugString" ) ),
  genSrc_     ( pset.getParameter< edm::InputTag >( "genSrc" ) ),
  genJetSrc_  ( pset.getParameter< edm::InputTag >( "genJetSrc" ) ),
  genMETSrc_  ( pset.getParameter< edm::InputTag >( "genMETSrc" ) ),

  doPUReweight_ ( pset.getParameter< bool >( "doPUReweight" ) ),
  puWeightSrc_( pset.getParameter< edm::InputTag >( "puWeights" ) ),
  eventWeightSrc_( pset.getParameter< edm::InputTag >( "eventWeights" ) ),

  studyAcc_     ( pset.getParameter< bool >( "studyAcceptance" ) ),
  studyRecoIso_ ( pset.getParameter< bool >( "studyRecoIso" ) ),

  bosonMinPt_   ( pset.getParameter< double >( "bosonMinPt" ) ),
  bosonEBMaxEta_( pset.getParameter< double >( "bosonEBMaxEta" ) ),
  bosonEEMinEta_( pset.getParameter< double >( "bosonEEMinEta" ) ),
  bosonEEMaxEta_( pset.getParameter< double >( "bosonEEMaxEta" ) )
{
  vertexSrc_      = pset.getParameter<edm::InputTag>("VertexSrc");

  if (studyRecoIso_) 
    recoPhotonSrc_ = pset.getParameter< edm::InputTag >( "recoPhotonSrc");
  recoJetSrc_    = pset.getParameter< edm::InputTag >( "recoJetSrc" );
  htJetSrc_      = pset.getParameter< edm::InputTag >( "htJetSrc" );
  bJetSrc_       = pset.getParameter< edm::InputTag >( "bJetSrc" );
  htSrc_         = pset.getParameter<edm::InputTag>("htSource");
  mhtSrc_        = pset.getParameter<edm::InputTag>("mhtSource");
  metSrc_        = pset.getParameter<edm::InputTag>("metSource");
  htNoBosonSrc_  = pset.getParameter<edm::InputTag>("htNoBosonSource");
  mhtNoBosonSrc_ = pset.getParameter<edm::InputTag>("mhtNoBosonSource");
  
}


GenStudyTree::~GenStudyTree()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void GenStudyTree::produce(edm::Event& ev, const edm::EventSetup& es)
{
  using namespace edm;
  //read in the gen particles
  edm::Handle<reco::GenParticleCollection> gens;
  ev.getByLabel(genSrc_,gens);

  //get the jets
  edm::Handle<reco::GenJetCollection> genJets;
  ev.getByLabel(genJetSrc_,genJets);

  //get the METs
  edm::Handle<edm::View<reco::GenMET> > genMET;
  ev.getByLabel(genMETSrc_,genMET);

  ///Reco quantities

  edm::Handle<reco::VertexCollection > vertices;
  ev.getByLabel(vertexSrc_, vertices);
  if (debug_)
    std::cout<<"vertex collection has size "<<vertices->size()<<std::endl;

  //get the reco products
  edm::Handle<edm::View<pat::Photon> > recoPhotons;
  edm::Handle<edm::View<pat::Jet> > recoJets;
  edm::Handle<edm::View<pat::Jet> > htJets;
  edm::Handle<edm::View<pat::Jet> > bJets;
  edm::Handle<double > ht;
  edm::Handle<edm::View<reco::MET> > mht;
  edm::Handle<double > htNoBoson;
  edm::Handle<edm::View<reco::MET> > mhtNoBoson;
  edm::Handle<edm::View<reco::MET> > met;

  if (studyRecoIso_) {
    ev.getByLabel(recoPhotonSrc_,recoPhotons);
    //if (debug_)
    //  std::cout<<"RecoIso Photon collection has size "<<recoPhotons->size()<<std::endl;
  }
  //get the jets
  ev.getByLabel(recoJetSrc_,recoJets);
  //if (debug_)
  //  std::cout<<"all-Jet collection has size "<<recoJets->size()<<std::endl;
  ev.getByLabel(htJetSrc_,htJets);
  //if (debug_)
  //  std::cout<<"ht-Jet collection has size "<<htJets->size()<<std::endl;
  ev.getByLabel(bJetSrc_,bJets);
  //if (debug_)
  //  std::cout<<"b-Jet collection has size "<<bJets->size()<<std::endl;
  ev.getByLabel(htSrc_, ht);
  //if (debug_)
  //  std::cout<<"HT value "<<*ht<<std::endl;
  ev.getByLabel(mhtSrc_, mht);
  //if (debug_)
  //    std::cout<<"MHT value "<<(*mht)[0].pt()<<std::endl;
  ev.getByLabel(htNoBosonSrc_, htNoBoson);
  //if (debug_)
  //  std::cout<<"HT no boson value "<<*htNoBoson<<std::endl;
  ev.getByLabel(mhtNoBosonSrc_, mhtNoBoson);
  //if (debug_)
  //  std::cout<<"MHT no boson value "<<(*mhtNoBoson)[0].pt()<<std::endl;
  ev.getByLabel(metSrc_, met);
  //if (debug_)
  //    std::cout<<"MET value "<<(*met)[0].pt()<<std::endl;
  
  double pu_event_wt = 1.;
  edm::Handle<double> puWeight;
  if (doPUReweight_) {
    ev.getByLabel(puWeightSrc_,puWeight);
    pu_event_wt = *puWeight;
  }

  double event_wt = 1.;
  edm::Handle<double> eventWeight;
  ev.getByLabel(eventWeightSrc_,eventWeight);
  event_wt = *eventWeight;

  //m_EventWt = scale_;
  m_EventWt = event_wt;
  m_PUWt    = pu_event_wt;
  m_Vertices = vertices->size();
  
  //do the analysis on the gen bosons
  m_genBoson1Pt  = -10.;  //m_genBoson2Pt  = -10.;
  m_genBoson1Eta = -10.;  //m_genBoson2Eta = -10.;
  m_genBoson1M   = -10.;  //m_genBoson2M   = -10.;
  m_genBoson1MinDR = 10.;

  m_genBosons = gens->size();
  if (m_genBosons < 1) {
    std::cout<<debugString_<<"::Unable to find any gen bosons"<<std::endl;
    return;
  }

  //compute event varibles from gen jets
  double htGenJets(0.);//, mhtGenJets(0.);
  m_genHTMInv = 0.;
  //double htNoBosonGenJets(0.), mhtNoBosonGenJets(0.);
  //double mhtGenJetsNoPhot(0.), meffGenJetsNoPhot(0.);
  //double dPhi1(9.0), dPhi2(9.0), dPhi3(9.0);
  //double dPhi1NoPhot(9.0), dPhi2NoPhot(9.0), dPhi3NoPhot(9.0);
  std::vector<const reco::GenJet*> htJetsGen;
  std::vector<const reco::GenJet*> mhtJetsGen;
  std::vector<const reco::GenJet*> bJetsGen;
  reco::MET::LorentzVector mhtGen(0,0,0,0);
  reco::MET::LorentzVector mhtGenNoBoson(0,0,0,0);
  reco::GenJetCollection::const_iterator jet = genJets->begin();

  m_nJetsGenPt30Eta50 = 0; m_nJetsGenPt30Eta24 = 0;
  m_bJetsGenPt30Eta24 = 0; m_nJetsGenPt50Eta25 = 0;
  m_nJetsGenPt50Eta25MInv = 0;

  int jetNum = 0;
  for (; jet != genJets->end(); ++jet) {
    if (jet->pt() > 30) {
      //++nJetsPt30;
      if (fabs(jet->eta()) < 5.0) {
	++m_nJetsGenPt30Eta50 ;
	mhtGen -= jet->p4();
	mhtGenNoBoson -= jet->p4();
	mhtJetsGen.push_back(&((*genJets)[jetNum]));
	if (jet->pt() > 50. && fabs(jet->eta()) < 2.5) {
	  ++m_nJetsGenPt50Eta25;
	  htGenJets += jet->pt();
	  htJetsGen.push_back(&((*genJets)[jetNum]));
	  if (gens->size() > 0 )
	    if ((jet->p4()+(*gens)[0].p4()).mass() > 90.0) {
	      m_genHTMInv  += jet->pt();
	      ++m_nJetsGenPt50Eta25MInv;
	    }
	}
	if (fabs(jet->eta()) < 2.4) {
	  ++m_nJetsGenPt30Eta24;
	  //find the b-jets
	  std::vector<const reco::GenParticle*> jetParts = jet->getGenConstituents();
	  bool hasB = false;
	  std::vector<const reco::GenParticle*>::const_iterator part = jetParts.begin();
	  //if (debug_) 
	  //  std::cout<<"jet has "<<jetParts.size()<<" particle constituents"<<std::endl;
	  for (; part != jetParts.end(); ++part) {
	    //const reco::GenParticle* constituent(*part);
	    //if (fabs(constituent->pdgId()) == 5)

	    const reco::GenParticle* p = *part;
	    while ( true ) {
		if ( !p ) break;
		
		if ( p->status() == 3 && fabs(p->pdgId()) == 5 ) {
		    // This is B jet: A jet constituent came from b-quark in hard process
		  hasB = true;
		  break;
		}
		else {
		  // Continue following up mothers
		  p = dynamic_cast<const reco::GenParticle*>(p->mother());
		}
	    }
	    
	    
	    //if (debug_) 
	    //  std::cout<<"jet has particle constituent "<<(*part)->pdgId()<<std::endl;
	    //if (fabs((*part)->pdgId()) == 5)
	    //  hasB = true;
	  }
	  if (hasB) {
	    bJetsGen.push_back(&((*genJets)[jetNum]));
	    ++m_bJetsGenPt30Eta24;
	  }
	}
      }
    }
    ++jetNum;
  }

  reco::MET MHT = reco::MET(mhtGen, reco::MET::Point());
  m_genMHT = MHT.pt();
  m_genMET = (*genMET)[0].pt();
  m_genMETNoBoson = ((*genMET)[0].p4()+(*gens)[0].p4()).pt();
  m_genHT  = htGenJets;
  //remove leading boson from the mht
  //if (removePhot_ && (gens->size() > 0))
  //  mhtNoPhot -= (*gens)[0].p4();

  //reco::MET MHTNoPhot = reco::MET(mhtNoPhot, reco::MET::Point());
  //mhtGenJetsNoPhot = MHTNoPhot.pt();
  ////meffGenJetsNoPhot = mhtGenJetsNoPhot + htGenJets;
  m_genDPhiMHT1 = 10.0;   m_genDPhiMET1 = 10.0;   m_genDPhiMETNoBoson1 = 10.0;
  m_genDPhiMHT2 = 10.0;	  m_genDPhiMET2 = 10.0;	  m_genDPhiMETNoBoson2 = 10.0;
  m_genDPhiMHT3 = 10.0;	  m_genDPhiMET3 = 10.0;	  m_genDPhiMETNoBoson3 = 10.0;
  m_genDPhiMHT4 = 10.0;	  m_genDPhiMET4 = 10.0;	  m_genDPhiMETNoBoson4 = 10.0;
  m_genJet1Pt  = -10.;  m_genJet3Pt  = -10.;
  m_genJet1Eta = -10.;  m_genJet3Eta = -10.;
  m_genJet2Pt  = -10.;  m_genJet4Pt  = -10.;
  m_genJet2Eta = -10.;  m_genJet4Eta = -10.;

  
  if (mhtJetsGen.size() > 0) {
    m_genDPhiMHT1 = fabs(reco::deltaPhi(mhtJetsGen.at(0)->phi(), MHT.phi()));
    m_genDPhiMET1 = fabs(reco::deltaPhi(mhtJetsGen.at(0)->phi(), (*genMET)[0].phi()));
    m_genDPhiMETNoBoson1 = fabs(reco::deltaPhi(mhtJetsGen.at(0)->phi(), ((*genMET)[0].p4()+(*gens)[0].p4()).phi()));
    m_genJet1Pt = mhtJetsGen.at(0)->pt();
    m_genJet1Eta = mhtJetsGen.at(0)->eta();
    //m_genDPhiMHT1NoPhot = fabs(reco::deltaPhi(mhtJetsGen.at(0)->phi(), MHTNoPhot.phi()));
    if (mhtJetsGen.size() > 1) {
      m_genDPhiMHT2 = fabs(reco::deltaPhi(mhtJetsGen.at(1)->phi(), MHT.phi()));
      m_genDPhiMET2 = fabs(reco::deltaPhi(mhtJetsGen.at(1)->phi(), (*genMET)[0].phi()));
      m_genDPhiMETNoBoson2 = fabs(reco::deltaPhi(mhtJetsGen.at(0)->phi(), ((*genMET)[0].p4()+(*gens)[0].p4()).phi()));
      m_genJet2Pt = mhtJetsGen.at(1)->pt();
      m_genJet2Eta = mhtJetsGen.at(1)->eta();
      //m_genDPhiMHT2NoPhot = fabs(reco::deltaPhi(mhtJetsGen.at(1)->phi(), MHTNoPhot.phi()));
      if (mhtJetsGen.size() > 2) {
	m_genDPhiMHT3 = fabs(reco::deltaPhi(mhtJetsGen.at(2)->phi(), MHT.phi()));
	m_genDPhiMET3 = fabs(reco::deltaPhi(mhtJetsGen.at(2)->phi(), (*genMET)[0].phi()));
	m_genDPhiMETNoBoson3 = fabs(reco::deltaPhi(mhtJetsGen.at(0)->phi(), ((*genMET)[0].p4()+(*gens)[0].p4()).phi()));
	m_genJet3Pt = mhtJetsGen.at(2)->pt();
	m_genJet3Eta = mhtJetsGen.at(2)->eta();
	//m_genDPhiMHT3NoPhot = fabs(reco::deltaPhi(mhtJetsGen.at(2)->phi(), MHTNoPhot.phi()));
	if (mhtJetsGen.size() > 3) {
	  m_genDPhiMHT4 = fabs(reco::deltaPhi(mhtJetsGen.at(3)->phi(), MHT.phi()));
	  m_genDPhiMET4 = fabs(reco::deltaPhi(mhtJetsGen.at(3)->phi(), (*genMET)[0].phi()));
	  m_genDPhiMETNoBoson4 = fabs(reco::deltaPhi(mhtJetsGen.at(0)->phi(), ((*genMET)[0].p4()+(*gens)[0].p4()).phi()));
	  m_genJet4Pt = mhtJetsGen.at(3)->pt();
	  m_genJet4Eta = mhtJetsGen.at(3)->eta();
	  //m_genDPhiMHT4NoPhot = fabs(reco::deltaPhi(mhtJetsGen.at(3)->phi(), MHTNoPhot.phi()));
	}
      }
    }
  }

  if (debug_) {
    reco::GenParticleCollection::const_iterator genp = gens->begin();
    for (; genp!= gens->end(); ++genp) {
      std::cout<<"pt("<<genp->pt()<<"), eta("<<genp->eta()<<"), phi("<<genp->phi()<<")"<<std::endl;
      std::cout<<"status("<<genp->status()<<"), pdgId("<<genp->pdgId()<<"), numberOfDaughters("<<genp->numberOfDaughters()<<")"<<std::endl;
      reco::GenParticle::const_iterator daus = genp->begin();
      std::cout<<"list of daughters:"<<std::endl;
      for (; daus!= genp->end(); ++daus) {
	std::cout<<"pt("<<daus->pt()<<"), eta("<<daus->eta()<<"), phi("<<daus->phi()<<")"<<std::endl;
	std::cout<<"status("<<daus->status()<<"), pdgId("<<daus->pdgId()<<"), numberOfDaughters("<<daus->numberOfDaughters()<<")"<<std::endl;
      }
    }
    if (studyRecoIso_) 
      std::cout<<"RecoIso Photon collection has size "<<recoPhotons->size()<<std::endl;
    std::cout<<"gen-Jet collection has size "<<genJets->size()<<std::endl;
    reco::GenJetCollection::const_iterator genj = genJets->begin();
    for (; genj!= genJets->end(); ++genj) 
      std::cout<<"pt("<<genj->pt()<<"), eta("<<genj->eta()<<"), phi("<<genj->phi()<<")"<<std::endl;
    std::cout<<"all-Jet collection has size "<<recoJets->size()<<std::endl;
    edm::View<pat::Jet>::const_iterator recoj = recoJets->begin();
    for (; recoj!= recoJets->end(); ++recoj) 
      std::cout<<"pt("<<recoj->pt()<<"), eta("<<recoj->eta()<<"), phi("<<recoj->phi()<<")"<<std::endl;
    std::cout<<"ht-Jet collection has size "<<htJets->size()<<std::endl;
    //edm::View<pat::Jet>::const_iterator htj = htJets->begin();
    //for (; htj!= htJets->end(); ++htj) 
    //  std::cout<<"pt("<<htj->pt()<<"), eta("<<htj->eta()<<"), phi("<<htj->phi()<<")"<<std::endl;
    std::cout<<"b-Jet collection has size "<<bJets->size()<<std::endl;
    std::cout<<"HT value "<<*ht<<std::endl;
    std::cout<<"MHT value "<<(*mht)[0].pt()<<std::endl;
    std::cout<<"MET value "<<(*met)[0].pt()<<std::endl;
  }
  //bool passAcc = false;
  //bool passRecoIso = false;
  if (m_genBosons > 0) {
    m_genBoson1Pt  = (*gens)[0].pt();
    m_genBoson1Eta = (*gens)[0].eta();
    double m_genBoson1Phi = (*gens)[0].phi();
    m_genBoson1M   = (*gens)[0].mass();
    
    if (debug_) std::cout<<"found gen boson"<<std::endl;
    jet = genJets->begin();
    for (; jet != genJets->end(); ++jet) {
      if (debug_) std::cout<<"looping gen jets"<<std::endl;
      if (jet->pt() > 30) {
	if (debug_) std::cout<<"pt > 30"<<std::endl;
	//++nJetsPt30;
	if (fabs(jet->eta()) < 5.0) {
	  if (debug_) std::cout<<"|eta| < 5.0"<<std::endl;
	  double dR = reco::deltaR((*gens)[0].eta(),(*gens)[0].phi(),jet->eta(), jet->phi());
	  if (debug_) std::cout<<"dR::"<<dR<<std::endl;
	  if (m_genBoson1MinDR > dR) 
	    m_genBoson1MinDR = dR;
	}
      }
    }
    //if (m_genBosons > 1) {
    //  m_genBoson2Pt  = (*gens)[1].pt();
    //  m_genBoson2Eta = (*gens)[1].eta();
    //  m_genBoson2M   = (*gens)[1].mass();
    //}
    
    
    if (m_genBoson1Pt > bosonMinPt_ && 
	(fabs(m_genBoson1Eta) < bosonEBMaxEta_ || 
	(fabs(m_genBoson1Eta) > bosonEEMinEta_ && 
	 fabs(m_genBoson1Eta) < bosonEEMaxEta_ )
	)
       )
      m_genPassAcc = true;
    //if (!studyAcc_)
    //  m_genPassAcc = true;
    m_daughter1Pt = -10;    m_daughter1Eta = -10;    m_daughter1M = -10;   m_daughter1ID = 0;
    m_daughter2Pt = -10;    m_daughter2Eta = -10;    m_daughter2M = -10;   m_daughter2ID = 0;
    m_combDaughterPt = -10; m_combDaughterEta = -10; m_combDaughterM = -10;

    reco::GenParticle::const_iterator daus = (*gens)[0].begin();
    const reco::Candidate *dau1,*dau2;
    int daucounter = 0;
    int founddaus = 0;
    //for (; daus!= (*gens)[0].end(); ++daus) {
    for (unsigned int dau = 0; dau < (*gens)[0].numberOfDaughters(); ++dau) {
      if ((*gens)[0].daughter(dau)->status()==3) {
	if (founddaus==0) 
	  dau1 = (*gens)[0].daughter(dau);
	if (founddaus==1) 
	  dau2 = (*gens)[0].daughter(dau);
	++founddaus;
      }
      ++daucounter;
    }
    if (founddaus == 2) {
      m_daughter1Pt = dau1->pt();    m_daughter1Eta = dau1->eta();    m_daughter1M = dau1->mass();   m_daughter1ID = dau1->pdgId();
      m_daughter2Pt = dau2->pt();    m_daughter2Eta = dau2->eta();    m_daughter2M = dau2->mass();   m_daughter2ID = dau2->pdgId();
      m_combDaughterPt = (dau1->p4()+dau2->p4()).pt(); m_combDaughterEta = (dau1->p4()+dau2->p4()).eta(); m_combDaughterM = (dau1->p4()+dau2->p4()).mass();
      
    }
    //Match gen boson to reco/isolated photon
    reco::Candidate::LorentzVector metNoBosonV = (*met)[0].p4();
    reco::MET metNoBoson = (*met)[0];
    m_genPassRecoIso = false;
    m_gen1PassRecoIso = false;
    if (!studyRecoIso_) {
      m_genPassRecoIso = true;
      m_gen1PassRecoIso = true;
    }
    else {

      m_boson1Pt  = -10.;  //m_boson2Pt  = -10.;
      m_boson1Eta = -10.;  //m_boson2Eta = -10.;
      m_boson1M   = -10.;  //m_boson2M   = -10.;
      m_boson1MinDR = 10.;
      m_nBosons   = recoPhotons->size();
      if (debug_)
	std::cout<<"found "<<m_nBosons<<" reconstructed photons"<<std::endl;

      if (recoPhotons->size()>0) {
	metNoBosonV += (*recoPhotons)[0].p4();
	metNoBoson.setP4(metNoBoson.p4() + (*recoPhotons)[0].p4());

	edm::View<pat::Jet>::const_iterator jet = recoJets->begin();
	for (; jet != recoJets->end(); ++jet) {
	  if (recoPhotons->size() > 0){
	    double dR = reco::deltaR((*recoPhotons)[0].eta(),(*recoPhotons)[0].phi(),jet->eta(), jet->phi());
	    if (m_boson1MinDR > dR)
	      m_boson1MinDR = dR;
	  }
	}
	jet = htJets->begin();
	for (; jet != htJets->end(); ++jet) {
	  if ((jet->p4()+(*recoPhotons)[0].p4()).mass() > 90.0) {
	    m_HTMInv += jet->pt();
	    ++m_nJetsPt50Eta25MInv;
	  }
	}
	
	m_boson1Pt  = (*recoPhotons)[0].pt();
	m_boson1Eta = (*recoPhotons)[0].eta();
	m_boson1M   = (*recoPhotons)[0].mass();
	if (debug_) {
	  std::cout<<"pt1::"<<m_boson1Pt
		   <<"eta1::"<<m_boson1Eta
		   <<std::endl;
	}
	
	//if (recoPhotons->size()>1) {
	//  m_boson2Pt  = (*recoPhotons)[1].pt();
	//  m_boson2Eta = (*recoPhotons)[1].eta();
	//  m_boson2M   = (*recoPhotons)[1].mass();
	//  if (debug_) {
	//    std::cout<<"pt2::"<<m_boson2Pt
	//	     <<"eta2::"<<m_boson2Eta
	//	     <<std::endl;
	//  }
	//}
      }

      //int bestDRPhot = -1;
      reco::GenParticleCollection::const_iterator genp = gens->begin();
      int gphot = 0;
      for (; genp != gens->end(); ++genp) {
	double bestDRMin = 999.0;
	int phot = 0;
	//here or outside the genp loop?
	edm::View<pat::Photon>::const_iterator recop = recoPhotons->begin();
	for (; recop != recoPhotons->end(); ++recop){
	  double dR = reco::deltaR(genp->eta(),genp->phi(),recop->eta(), recop->phi());
	  //double dR = reco::deltaR(m_genBoson1Eta,m_genBoson1Phi,recop->eta(), recop->phi());
	  if (dR < bestDRMin) {
	    //bestDRPhot = phot;
	    bestDRMin = dR;
	  }
	  ++phot;
	}
	if (bestDRMin < 0.2) {
	  //recoMatched = &((*recoPhotons)[bestDRPhot]);
	  m_genPassRecoIso = true;
	}
	if (m_genPassRecoIso && gphot==0)
	  m_gen1PassRecoIso = true;
	++gphot;
      }//end loop over gen particles
    }
    
    //    //remove the leading boson in the case of photons
    //    double htGenJetsNoPhot(htGenJets), mhtGenJetsNoPhot(0.), meffGenJetsNoPhot(0.);
    //    reco::MET::LorentzVector mhtNoPhot = mht;
    //    mhtNoPhot -= (*gens)[0].p4();
    //    reco::MET MHTNoPhot = reco::MET(mhtNoPhot, reco::MET::Point());
    //    mhtGenJetsNoPhot = MHTNoPhot.pt();
    //    meffGenJetsNoPhot = htGenJetsNoPhot + mhtGenJetsNoPhot;
    m_nJetsPt30Eta50 = recoJets->size();
    m_nJetsPt50Eta25 = htJets  ->size();

    m_bJetsPt30Eta24 = bJets   ->size();
    m_nJetsPt30Eta24 = 0;
    m_nJetsPt50Eta25MInv = 0;

    m_HT  = *ht;
    m_MHT = (*mht)[0].pt();
    m_MET = (*met)[0].pt();
    m_METNoBoson = metNoBosonV.pt();
    //m_METNoBoson = metNoBoson.pt();
    
    m_HTMInv  = 0;
    
    edm::View<pat::Jet>::const_iterator jet = recoJets->begin();
    for (; jet != recoJets->end(); ++jet) {
      if (fabs(jet->eta()) < 2.4)
	++m_nJetsPt30Eta24;
    }
    const pat::Jet  *r1, *r2, *r3, *r4;
    m_dPhiMHT1 = 10.0;    m_dPhiMET1 = 10.0;    m_dPhiMETNoBoson1 = 10.0;
    m_dPhiMHT2 = 10.0;    m_dPhiMET2 = 10.0;    m_dPhiMETNoBoson2 = 10.0;
    m_dPhiMHT3 = 10.0;    m_dPhiMET3 = 10.0;    m_dPhiMETNoBoson3 = 10.0;
    m_dPhiMHT4 = 10.0;    m_dPhiMET4 = 10.0;    m_dPhiMETNoBoson4 = 10.0;
    m_Jet1Pt  = -10.;    m_Jet3Pt  = -10.;
    m_Jet1Eta = -10.;    m_Jet3Eta = -10.;
    m_Jet2Pt  = -10.;    m_Jet4Pt  = -10.;
    m_Jet2Eta = -10.;    m_Jet4Eta = -10.;

    if(m_nJetsPt30Eta50 >= 1) {
      r1 = &((*recoJets)[0]);
      m_dPhiMHT1 = fabs(reco::deltaPhi(r1->phi(),(*mht)[0].phi()));
      m_dPhiMET1 = fabs(reco::deltaPhi(r1->phi(),(*met)[0].phi()));
      m_dPhiMETNoBoson1 = fabs(reco::deltaPhi(r1->phi(),metNoBoson.phi()));
      m_Jet1Pt = r1->pt();
      m_Jet1Eta = r1->eta();
      if(m_nJetsPt30Eta50 >= 2) {
	r2 = &((*recoJets)[1]);
	m_dPhiMHT2 = fabs(reco::deltaPhi(r2->phi(),(*mht)[0].phi())); 
	m_dPhiMET2 = fabs(reco::deltaPhi(r2->phi(),(*met)[0].phi())); 
	m_dPhiMETNoBoson2 = fabs(reco::deltaPhi(r1->phi(),metNoBoson.phi()));
	m_Jet2Pt = r2->pt();
	m_Jet2Eta = r2->eta();
	
	if(m_nJetsPt30Eta50 >= 3) {
	  r3 = &((*recoJets)[2]);
	  m_dPhiMHT3 = fabs(reco::deltaPhi(r3->phi(),(*mht)[0].phi()));
	  m_dPhiMET3 = fabs(reco::deltaPhi(r3->phi(),(*met)[0].phi()));
	  m_dPhiMETNoBoson3 = fabs(reco::deltaPhi(r1->phi(),metNoBoson.phi()));
	  m_Jet3Pt = r3->pt();
	  m_Jet3Eta = r3->eta();

	  if(m_nJetsPt30Eta50 >= 4) {
	    r4 = &((*recoJets)[3]);
	    m_dPhiMHT4 = fabs(reco::deltaPhi(r4->phi(),(*mht)[0].phi()));
	    m_dPhiMET4 = fabs(reco::deltaPhi(r4->phi(),(*met)[0].phi()));
	    m_dPhiMETNoBoson4 = fabs(reco::deltaPhi(r1->phi(),metNoBoson.phi()));
	    m_Jet4Pt = r4->pt();
	    m_Jet4Eta = r4->eta();
	  }
	}
      }
    }
    
  }
  //if (reducedValues)
  if (debug_ && (m_genBoson1MinDR>10 || m_genBoson1MinDR < 0)) std::cout<<"m_genBoson1MinDR::"<<m_genBoson1MinDR<<std::endl;
  if (debug_ && (m_boson1MinDR>10 || m_boson1MinDR < 0)) std::cout<<"m_boson1MinDR::"<<m_boson1MinDR<<std::endl;
  reducedValues->Fill();
}
// ------------ method called once each job just before starting event loop  ------------
void GenStudyTree::beginJob()
{
  edm::Service< TFileService > fs;
  
  reducedValues = fs->make<TTree>( "RA2Values", "Variables for reduced studies" );

  reducedValues->Branch("ra2_Vertices", &m_Vertices, "m_Vertices/I");
  reducedValues->Branch("ra2_PUWt",     &m_PUWt,     "m_PUWt/D");
  reducedValues->Branch("ra2_EventWt",  &m_EventWt,  "m_EventWt/D");

  ///Generator level quantities
  reducedValues->Branch("ra2_genHT",    &m_genHT,    "m_genHT/D" );
  reducedValues->Branch("ra2_genMHT",   &m_genMHT,   "m_genMHT/D");
  reducedValues->Branch("ra2_genMET",   &m_genMET,   "m_genMET/D");
  reducedValues->Branch("ra2_genMETNoBoson",   &m_genMETNoBoson,   "m_genMETNoBoson/D");
  reducedValues->Branch("ra2_genHTMInv",    &m_genHTMInv,    "m_genHTMInv/D" );

  reducedValues->Branch("ra2_genBosons",   &m_genBosons,   "m_genBosons/I" );
  reducedValues->Branch("ra2_genBoson1Pt", &m_genBoson1Pt, "m_genBoson1Pt/D" );
  reducedValues->Branch("ra2_genBoson1Eta",&m_genBoson1Eta,"m_genBoson1Eta/D" );
  reducedValues->Branch("ra2_genBoson1M",  &m_genBoson1M,  "m_genBoson1M/D" );
  reducedValues->Branch("ra2_genBoson1MinDR",  &m_genBoson1MinDR,  "m_genBoson1MinDR/D" );
  //reducedValues->Branch("ra2_genBoson2Pt", &m_genBoson2Pt, "m_genBoson2Pt/D" );
  //reducedValues->Branch("ra2_genBoson2Eta",&m_genBoson2Eta,"m_genBoson2Eta/D" );
  //reducedValues->Branch("ra2_genBoson2M",  &m_genBoson2M,  "m_genBoson2M/D" );

  reducedValues->Branch("ra2_daughter1Pt", &m_daughter1Pt, "m_daughter1Pt/D" );
  reducedValues->Branch("ra2_daughter1Eta",&m_daughter1Eta,"m_daughter1Eta/D" );
  reducedValues->Branch("ra2_daughter1M",  &m_daughter1M,  "m_daughter1M/D" );
  reducedValues->Branch("ra2_daughter1ID",  &m_daughter1ID,  "m_daughter1ID/D" );
  reducedValues->Branch("ra2_daughter2Pt", &m_daughter2Pt, "m_daughter2Pt/D" );
  reducedValues->Branch("ra2_daughter2Eta",&m_daughter2Eta,"m_daughter2Eta/D" );
  reducedValues->Branch("ra2_daughter2M",  &m_daughter2M,  "m_daughter2M/D" );
  reducedValues->Branch("ra2_daughter2ID",  &m_daughter2ID,  "m_daughter2ID/D" );
  reducedValues->Branch("ra2_combDaughterPt", &m_combDaughterPt, "m_combDaughterPt/D" );
  reducedValues->Branch("ra2_combDaughterEta",&m_combDaughterEta,"m_combDaughterEta/D" );
  reducedValues->Branch("ra2_combDaughterM",  &m_combDaughterM,  "m_combDaughterM/D" );

  reducedValues->Branch("ra2_genDPhiMHT1", &m_genDPhiMHT1, "m_genDPhiMHT1/D");
  reducedValues->Branch("ra2_genDPhiMHT2", &m_genDPhiMHT2, "m_genDPhiMHT2/D");
  reducedValues->Branch("ra2_genDPhiMHT3", &m_genDPhiMHT3, "m_genDPhiMHT3/D");
  reducedValues->Branch("ra2_genDPhiMHT4", &m_genDPhiMHT4, "m_genDPhiMHT4/D");

  reducedValues->Branch("ra2_genDPhiMET1", &m_genDPhiMET1, "m_genDPhiMET1/D");
  reducedValues->Branch("ra2_genDPhiMET2", &m_genDPhiMET2, "m_genDPhiMET2/D");
  reducedValues->Branch("ra2_genDPhiMET3", &m_genDPhiMET3, "m_genDPhiMET3/D");
  reducedValues->Branch("ra2_genDPhiMET4", &m_genDPhiMET4, "m_genDPhiMET4/D");

  reducedValues->Branch("ra2_genDPhiMETNoBoson1", &m_genDPhiMETNoBoson1, "m_genDPhiMETNoBoson1/D");
  reducedValues->Branch("ra2_genDPhiMETNoBoson2", &m_genDPhiMETNoBoson2, "m_genDPhiMETNoBoson2/D");
  reducedValues->Branch("ra2_genDPhiMETNoBoson3", &m_genDPhiMETNoBoson3, "m_genDPhiMETNoBoson3/D");
  reducedValues->Branch("ra2_genDPhiMETNoBoson4", &m_genDPhiMETNoBoson4, "m_genDPhiMETNoBoson4/D");

  reducedValues->Branch("ra2_genJet1Pt",  &m_genJet1Pt,  "m_genJet1Pt/D");
  reducedValues->Branch("ra2_genJet1Eta", &m_genJet1Eta, "m_genJet1Eta/D");
  reducedValues->Branch("ra2_genJet2Pt",  &m_genJet2Pt,  "m_genJet2Pt/D");
  reducedValues->Branch("ra2_genJet2Eta", &m_genJet2Eta, "m_genJet2Eta/D");
  reducedValues->Branch("ra2_genJet3Pt",  &m_genJet3Pt,  "m_genJet3Pt/D");
  reducedValues->Branch("ra2_genJet3Eta", &m_genJet3Eta, "m_genJet3Eta/D");
  reducedValues->Branch("ra2_genJet4Pt",  &m_genJet4Pt,  "m_genJet4Pt/D");
  reducedValues->Branch("ra2_genJet4Eta", &m_genJet4Eta, "m_genJet4Eta/D");

  reducedValues->Branch("ra2_nJetsGenPt30Eta50", &m_nJetsGenPt30Eta50, "nJetsGenPt30Eta50/I" );
  reducedValues->Branch("ra2_nJetsGenPt30Eta24", &m_nJetsGenPt30Eta24, "nJetsGenPt30Eta24/I" );
  reducedValues->Branch("ra2_bJetsGenPt30Eta24", &m_bJetsGenPt30Eta24, "bJetsGenPt30Eta24/I");
  reducedValues->Branch("ra2_nJetsGenPt50Eta25", &m_nJetsGenPt50Eta25, "nJetsGenPt50Eta25/I" );
  reducedValues->Branch("ra2_nJetsGenPt50Eta25MInv", &m_nJetsGenPt50Eta25MInv, "nJetsGenPt50Eta25MInv/I" );

  ///Reco quantities
  reducedValues->Branch("ra2_HT",    &m_HT,    "m_HT/D" );
  reducedValues->Branch("ra2_MHT",   &m_MHT,   "m_MHT/D");
  reducedValues->Branch("ra2_MET",   &m_MET,   "m_MET/D");
  reducedValues->Branch("ra2_HTMInv",    &m_HTMInv,    "m_HTMInv/D" );
  
  reducedValues->Branch("ra2_genPassRecoIso",  &m_genPassRecoIso,  "m_genPassRecoIso/O");
  reducedValues->Branch("ra2_gen1PassRecoIso",  &m_gen1PassRecoIso,  "m_gen1PassRecoIso/O");
  if (studyRecoIso_) {
    reducedValues->Branch("ra2_nBosons",  &m_nBosons,  "m_nBosons/I" );
    reducedValues->Branch("ra2_boson1Pt", &m_boson1Pt, "m_boson1Pt/D" );
    reducedValues->Branch("ra2_boson1Eta",&m_boson1Eta,"m_boson1Eta/D" );
    reducedValues->Branch("ra2_boson1M",  &m_boson1M,  "m_boson1M/D" );
    reducedValues->Branch("ra2_boson1MinDR",  &m_boson1MinDR,  "m_boson1MinDR/D" );
    //reducedValues->Branch("ra2_boson2Pt", &m_boson2Pt, "m_boson2Pt/D" );
    //reducedValues->Branch("ra2_boson2Eta",&m_boson2Eta,"m_boson2Eta/D" );
    //reducedValues->Branch("ra2_boson2M",  &m_boson2M,  "m_boson2M/D" );

    reducedValues->Branch("ra2_METNoBoson",   &m_METNoBoson,   "m_METNoBoson/D");
    reducedValues->Branch("ra2_dPhiMETNoBoson1", &m_dPhiMETNoBoson1, "m_dPhiMETNoBoson1/D");
    reducedValues->Branch("ra2_dPhiMETNoBoson2", &m_dPhiMETNoBoson2, "m_dPhiMETNoBoson2/D");
    reducedValues->Branch("ra2_dPhiMETNoBoson3", &m_dPhiMETNoBoson3, "m_dPhiMETNoBoson3/D");
    reducedValues->Branch("ra2_dPhiMETNoBoson4", &m_dPhiMETNoBoson4, "m_dPhiMETNoBoson4/D");
  }
  reducedValues->Branch("ra2_dPhiMHT1", &m_dPhiMHT1, "m_dPhiMHT1/D");
  reducedValues->Branch("ra2_dPhiMHT2", &m_dPhiMHT2, "m_dPhiMHT2/D");
  reducedValues->Branch("ra2_dPhiMHT3", &m_dPhiMHT3, "m_dPhiMHT3/D");
  reducedValues->Branch("ra2_dPhiMHT4", &m_dPhiMHT4, "m_dPhiMHT4/D");
  
  reducedValues->Branch("ra2_dPhiMET1", &m_dPhiMET1, "m_dPhiMET1/D");
  reducedValues->Branch("ra2_dPhiMET2", &m_dPhiMET2, "m_dPhiMET2/D");
  reducedValues->Branch("ra2_dPhiMET3", &m_dPhiMET3, "m_dPhiMET3/D");
  reducedValues->Branch("ra2_dPhiMET4", &m_dPhiMET4, "m_dPhiMET4/D");
  
  reducedValues->Branch("ra2_Jet1Pt",  &m_Jet1Pt,  "m_Jet1Pt/D");
  reducedValues->Branch("ra2_Jet1Eta", &m_Jet1Eta, "m_Jet1Eta/D");
  reducedValues->Branch("ra2_Jet2Pt",  &m_Jet2Pt,  "m_Jet2Pt/D");
  reducedValues->Branch("ra2_Jet2Eta", &m_Jet2Eta, "m_Jet2Eta/D");
  reducedValues->Branch("ra2_Jet3Pt",  &m_Jet3Pt,  "m_Jet3Pt/D");
  reducedValues->Branch("ra2_Jet3Eta", &m_Jet3Eta, "m_Jet3Eta/D");
  reducedValues->Branch("ra2_Jet4Pt",  &m_Jet4Pt,  "m_Jet4Pt/D");
  reducedValues->Branch("ra2_Jet4Eta", &m_Jet4Eta, "m_Jet4Eta/D");
  
  reducedValues->Branch("ra2_nJetsPt30Eta50", &m_nJetsPt30Eta50, "nJetsPt30Eta50/I" );
  reducedValues->Branch("ra2_nJetsPt30Eta24", &m_nJetsPt30Eta24, "nJetsPt30Eta24/I" );
  reducedValues->Branch("ra2_bJetsPt30Eta24", &m_bJetsPt30Eta24, "bJetsPt30Eta24/I");
  reducedValues->Branch("ra2_nJetsPt50Eta25", &m_nJetsPt50Eta25, "nJetsPt50Eta25/I" );
  reducedValues->Branch("ra2_nJetsPt50Eta25MInv", &m_nJetsPt50Eta25MInv, "nJetsPt50Eta25MInv/I" );
  

  reducedValues->SetAutoSave(1);
}

// ------------ method called once each job just after ending the event loop  ------------
void GenStudyTree::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenStudyTree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenStudyTree);
