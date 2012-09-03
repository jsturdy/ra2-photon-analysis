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
// $Id: GenStudyTree.cc,v 1.1 2012/08/30 09:45:14 sturdy Exp $
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

  doPUReweight_ ( pset.getParameter< bool >( "doPUReweight" ) ),
  puWeightSrc_( pset.getParameter< edm::InputTag >( "puWeight" ) ),

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
  
  double pu_event_wt = 1.;
  edm::Handle<double> puWeight;
  if (doPUReweight_) {
    ev.getByLabel(puWeightSrc_,puWeight);
    pu_event_wt = *puWeight;
  }

  m_EventWt = scale_;
  m_PUWt    = pu_event_wt;
  m_Vertices = vertices->size();
  

  //compute event varibles from gen jets
  double htGenJets(0.);//, mhtGenJets(0.);
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
  m_genHT  = htGenJets;
  //remove leading boson from the mht
  //if (removePhot_ && (gens->size() > 0))
  //  mhtNoPhot -= (*gens)[0].p4();

  //reco::MET MHTNoPhot = reco::MET(mhtNoPhot, reco::MET::Point());
  //mhtGenJetsNoPhot = MHTNoPhot.pt();
  ////meffGenJetsNoPhot = mhtGenJetsNoPhot + htGenJets;
  m_genDPhi1 = -1.0;
  m_genDPhi2 = -1.0;
  m_genDPhi3 = -1.0;
  m_genJet1Pt  = -10.;
  m_genJet1Eta = -10.;
  m_genJet2Pt  = -10.;
  m_genJet2Eta = -10.;
  m_genJet3Pt  = -10.;
  m_genJet3Eta = -10.;
  
  if (mhtJetsGen.size() > 0) {
    m_genDPhi1 = fabs(reco::deltaPhi(mhtJetsGen.at(0)->phi(), MHT.phi()));
    m_genJet1Pt = mhtJetsGen.at(0)->pt();
    m_genJet1Eta = mhtJetsGen.at(0)->eta();
    //m_genDPhi1NoPhot = fabs(reco::deltaPhi(mhtJetsGen.at(0)->phi(), MHTNoPhot.phi()));
    if (mhtJetsGen.size() > 1) {
      m_genDPhi2 = fabs(reco::deltaPhi(mhtJetsGen.at(1)->phi(), MHT.phi()));
      m_genJet2Pt = mhtJetsGen.at(1)->pt();
      m_genJet2Eta = mhtJetsGen.at(1)->eta();
      //m_genDPhi2NoPhot = fabs(reco::deltaPhi(mhtJetsGen.at(1)->phi(), MHTNoPhot.phi()));
      if (mhtJetsGen.size() > 2) {
	m_genDPhi3 = fabs(reco::deltaPhi(mhtJetsGen.at(2)->phi(), MHT.phi()));
	m_genJet3Pt = mhtJetsGen.at(2)->pt();
	m_genJet3Eta = mhtJetsGen.at(2)->eta();
	//m_genDPhi3NoPhot = fabs(reco::deltaPhi(mhtJetsGen.at(2)->phi(), MHTNoPhot.phi()));
      }
    }
  }
  //do the analysis on the gen bosons
  m_genBoson1Pt  = -10.;  m_genBoson2Pt  = -10.;
  m_genBoson1Eta = -10.;  m_genBoson2Eta = -10.;
  m_genBoson1M   = -10.;  m_genBoson2M   = -10.;

  m_genBosons = gens->size();
  if (m_genBosons < 1) {
    std::cout<<debugString_<<"::Unable to find any gen bosons"<<std::endl;
    return;
  }

  if (debug_) {
    if (studyRecoIso_) 
      std::cout<<"RecoIso Photon collection has size "<<recoPhotons->size()<<std::endl;
    std::cout<<"gen-Jet collection has size "<<genJets->size()<<std::endl;
    //reco::GenJetCollection::const_iterator genj = genJets->begin();
    //for (; genj!= genJets->end(); ++genj) 
    //  std::cout<<"pt("<<genj->pt()<<"), eta("<<genj->eta()<<"), phi("<<genj->phi()<<")"<<std::endl;
    std::cout<<"all-Jet collection has size "<<recoJets->size()<<std::endl;
    //edm::View<pat::Jet>::const_iterator recoj = recoJets->begin();
    //for (; recoj!= recoJets->end(); ++recoj) 
    //  std::cout<<"pt("<<recoj->pt()<<"), eta("<<recoj->eta()<<"), phi("<<recoj->phi()<<")"<<std::endl;
    std::cout<<"ht-Jet collection has size "<<htJets->size()<<std::endl;
    //edm::View<pat::Jet>::const_iterator htj = htJets->begin();
    //for (; htj!= htJets->end(); ++htj) 
    //  std::cout<<"pt("<<htj->pt()<<"), eta("<<htj->eta()<<"), phi("<<htj->phi()<<")"<<std::endl;
    std::cout<<"b-Jet collection has size "<<bJets->size()<<std::endl;
    std::cout<<"HT value "<<*ht<<std::endl;
    std::cout<<"MHT value "<<(*mht)[0].pt()<<std::endl;
    std::cout<<"HT no boson value "<<*htNoBoson<<std::endl;
    std::cout<<"MHT no boson value "<<(*mhtNoBoson)[0].pt()<<std::endl;
  }
  //bool passAcc = false;
  //bool passRecoIso = false;
  if (m_genBosons > 0) {
    m_genBoson1Pt  = (*gens)[0].pt();
    m_genBoson1Eta = (*gens)[0].eta();
    double m_genBoson1Phi = (*gens)[0].phi();
    m_genBoson1M   = (*gens)[0].mass();

    if (m_genBosons > 1) {
      m_genBoson2Pt  = (*gens)[1].pt();
      m_genBoson2Eta = (*gens)[1].eta();
      m_genBoson2M   = (*gens)[1].mass();
    }
    
    
    if (m_genBoson1Pt > bosonMinPt_ && 
	(fabs(m_genBoson1Eta) < bosonEBMaxEta_ || 
	(fabs(m_genBoson1Eta) > bosonEEMinEta_ && 
	 fabs(m_genBoson1Eta) < bosonEEMaxEta_ )
	)
       )
      m_genPassAcc = true;
    //if (!studyAcc_)
    //  m_genPassAcc = true;
    
    //Match gen boson to reco/isolated photon
    if (!studyRecoIso_)
      m_genPassRecoIso = true;
    else {
      edm::View<pat::Photon>::const_iterator recop = recoPhotons->begin();

      m_boson1Pt  = -10.;  m_boson2Pt  = -10.;
      m_boson1Eta = -10.;  m_boson2Eta = -10.;
      m_boson1M   = -10.;  m_boson2M   = -10.;
      m_nBosons   = recoPhotons->size();
      
      if (recoPhotons->size()>0) {
	m_boson1Pt  = (*recoPhotons)[0].pt();
	m_boson1Eta = (*recoPhotons)[0].eta();
	m_boson1M   = (*recoPhotons)[0].mass();
	
	if (recoPhotons->size()>1) {
	  m_boson2Pt  = (*recoPhotons)[1].pt();
	  m_boson2Eta = (*recoPhotons)[1].eta();
	  m_boson2M   = (*recoPhotons)[1].mass();
	}
      }

      //int bestDRPhot = -1;
      double bestDRMin = 999.0;
      int phot = 0;
      for (; recop != recoPhotons->end(); ++recop){
	double dR = reco::deltaR(m_genBoson1Eta,m_genBoson1Phi,recop->eta(), recop->phi());
	if (dR < bestDRMin) {
	  //bestDRPhot = phot;
	  bestDRMin = dR;
	}
	++phot;
      }
      if (bestDRMin < 0.1) {
	//recoMatched = &((*recoPhotons)[bestDRPhot]);
	m_genPassRecoIso = true;
      }
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
    m_HT  = *ht;
    m_MHT = (*mht)[0].pt();
    m_HTNoBoson  = *htNoBoson;
    m_MHTNoBoson = (*mhtNoBoson)[0].pt();
    
    edm::View<pat::Jet>::const_iterator jet = recoJets->begin();
    for (; jet != recoJets->end(); ++jet) {
      if (fabs(jet->eta()) < 2.4)
	++m_nJetsPt30Eta24;
    }
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
      r1 = &((*recoJets)[0]);
      m_dPhi1 = fabs(reco::deltaPhi(r1->phi(),(*mht)[0].phi()));
      m_Jet1Pt = r1->pt();
      m_Jet1Eta = r1->eta();
      if(m_nJetsPt30Eta50 >= 2) {
	r2 = &((*recoJets)[1]);
	m_dPhi2 = fabs(reco::deltaPhi(r2->phi(),(*mht)[0].phi())); 
	m_Jet2Pt = r2->pt();
	m_Jet2Eta = r2->eta();
	
	if(m_nJetsPt30Eta50 >= 3) {
	  r3 = &((*recoJets)[2]);
	  m_dPhi3 = fabs(reco::deltaPhi(r3->phi(),(*mht)[0].phi()));
	  m_Jet3Pt = r3->pt();
	  m_Jet3Eta = r3->eta();
	}
      }
    }
    
  }
  //if (reducedValues)
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
  reducedValues->Branch("ra2_genHTNoBoson",    &m_genHTNoBoson,    "m_genHTNoBoson/D" );
  reducedValues->Branch("ra2_genMHTNoBoson",   &m_genMHTNoBoson,   "m_genMHTNoBoson/D");

  reducedValues->Branch("ra2_genBosons",   &m_genBosons,   "m_genBosons/I" );
  reducedValues->Branch("ra2_genBoson1Pt", &m_genBoson1Pt, "m_genBoson1Pt/D" );
  reducedValues->Branch("ra2_genBoson1Eta",&m_genBoson1Eta,"m_genBoson1Eta/D" );
  reducedValues->Branch("ra2_genBoson1M",  &m_genBoson1M,  "m_genBoson1M/D" );
  reducedValues->Branch("ra2_genBoson2Pt", &m_genBoson2Pt, "m_genBoson2Pt/D" );
  reducedValues->Branch("ra2_genBoson2Eta",&m_genBoson2Eta,"m_genBoson2Eta/D" );
  reducedValues->Branch("ra2_genBoson2M",  &m_genBoson2M,  "m_genBoson2M/D" );

  reducedValues->Branch("ra2_genDPhi1", &m_genDPhi1, "m_genDPhi1/D");
  reducedValues->Branch("ra2_genDPhi2", &m_genDPhi2, "m_genDPhi2/D");
  reducedValues->Branch("ra2_genDPhi3", &m_genDPhi3, "m_genDPhi3/D");
  //reducedValues->Branch("ra2_genDPhiNoBoson1", &m_genDPhiNoBoson1, "m_genDPhiNoBoson1/D");
  //reducedValues->Branch("ra2_genDPhiNoBoson2", &m_genDPhiNoBoson2, "m_genDPhiNoBoson2/D");
  //reducedValues->Branch("ra2_genDPhiNoBoson3", &m_genDPhiNoBoson3, "m_genDPhiNoBoson3/D");

  reducedValues->Branch("ra2_genJet1Pt",  &m_genJet1Pt,  "m_genJet1Pt/D");
  reducedValues->Branch("ra2_genJet1Eta", &m_genJet1Eta, "m_genJet1Eta/D");
  reducedValues->Branch("ra2_genJet2Pt",  &m_genJet2Pt,  "m_genJet2Pt/D");
  reducedValues->Branch("ra2_genJet2Eta", &m_genJet2Eta, "m_genJet2Eta/D");
  reducedValues->Branch("ra2_genJet3Pt",  &m_genJet3Pt,  "m_genJet3Pt/D");
  reducedValues->Branch("ra2_genJet3Eta", &m_genJet3Eta, "m_genJet3Eta/D");

  reducedValues->Branch("ra2_nJetsGenPt30Eta50", &m_nJetsGenPt30Eta50, "nJetsGenPt30Eta50/I" );
  reducedValues->Branch("ra2_nJetsGenPt30Eta24", &m_nJetsGenPt30Eta24, "nJetsGenPt30Eta24/I" );
  reducedValues->Branch("ra2_bJetsGenPt30Eta24", &m_bJetsGenPt30Eta24, "bJetsGenPt30Eta24/I");
  reducedValues->Branch("ra2_nJetsGenPt50Eta25", &m_nJetsGenPt50Eta25, "nJetsGenPt50Eta25/I" );

  ///Reco quantities
  if (studyRecoIso_) {
    reducedValues->Branch("ra2_HT",    &m_HT,    "m_HT/D" );
    reducedValues->Branch("ra2_MHT",   &m_MHT,   "m_MHT/D");
    reducedValues->Branch("ra2_HTNoBoson",    &m_HTNoBoson,    "m_HTNoBoson/D" );
    reducedValues->Branch("ra2_MHTNoBoson",   &m_MHTNoBoson,   "m_MHTNoBoson/D");
    
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
    //reducedValues->Branch("ra2_dPhiNoBoson1", &m_dPhiNoBoson1, "m_dPhiNoBoson1/D");
    //reducedValues->Branch("ra2_dPhiNoBoson2", &m_dPhiNoBoson2, "m_dPhiNoBoson2/D");
    //reducedValues->Branch("ra2_dPhiNoBoson3", &m_dPhiNoBoson3, "m_dPhiNoBoson3/D");
    
    reducedValues->Branch("ra2_Jet1Pt",  &m_Jet1Pt,  "m_Jet1Pt/D");
    reducedValues->Branch("ra2_Jet1Eta", &m_Jet1Eta, "m_Jet1Eta/D");
    reducedValues->Branch("ra2_Jet2Pt",  &m_Jet2Pt,  "m_Jet2Pt/D");
    reducedValues->Branch("ra2_Jet2Eta", &m_Jet2Eta, "m_Jet2Eta/D");
    reducedValues->Branch("ra2_Jet3Pt",  &m_Jet3Pt,  "m_Jet3Pt/D");
    reducedValues->Branch("ra2_Jet3Eta", &m_Jet3Eta, "m_Jet3Eta/D");
    
    reducedValues->Branch("ra2_nJetsPt30Eta50", &m_nJetsPt30Eta50, "nJetsPt30Eta50/I" );
    reducedValues->Branch("ra2_nJetsPt30Eta24", &m_nJetsPt30Eta24, "nJetsPt30Eta24/I" );
    reducedValues->Branch("ra2_bJetsPt30Eta24", &m_bJetsPt30Eta24, "bJetsPt30Eta24/I");
    reducedValues->Branch("ra2_nJetsPt50Eta25", &m_nJetsPt50Eta25, "nJetsPt50Eta25/I" );
  }

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
