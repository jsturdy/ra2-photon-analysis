// -*- C++ -*-
//
// Package:    Photons
// Class:      Photons
// 
/**\class Photons GenStudy.cc ZInvisibleBkgds/Photons/plugins/GenStudy.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: GenStudy.cc,v 1.2 2012/07/20 11:34:12 sturdy Exp $
//
//
// system include files
#include <cmath>

#include "ZInvisibleBkgds/Photons/interface/GenStudy.h"

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
GenStudy::GenStudy(const edm::ParameterSet& pset) :
  debug_        ( pset.getParameter< bool >( "debug" ) ),
  debugString_  ( pset.getParameter< std::string >( "debugString" ) ),
  genLabel_     ( pset.getParameter< edm::InputTag >( "genLabel" ) ),
  pdgId_        ( pset.getParameter< int >( "pdgId" ) ),
  genStatus_    ( pset.getParameter< int >( "genStatus" ) ),
  mom_pdgId_    ( pset.getParameter< int >( "mompdgId" ) ),
  mom_genStatus_( pset.getParameter< int >( "momgenStatus" ) ),
  genJetLabel_  ( pset.getParameter< edm::InputTag >( "genJetLabel" ) ),

  doPUReweight_ ( pset.getParameter< bool >( "doPUReweight" ) ),
  puWeightLabel_( pset.getParameter< edm::InputTag >( "puWeight" ) ),

  minHT_ ( pset.getParameter< double >( "minHT" ) ),
  minMHT_( pset.getParameter< double >( "minMHT" ) ),
  
  studyAcc_     ( pset.getParameter< bool >( "studyAcceptance" ) ),
  studyRecoIso_ ( pset.getParameter< bool >( "studyRecoIso" ) ),
  removePhot_   ( pset.getParameter< bool >( "removePhoton" ) ),
  bosonMinPt_   ( pset.getParameter< double >( "bosonMinPt" ) ),
  bosonEBMaxEta_( pset.getParameter< double >( "bosonEBMaxEta" ) ),
  bosonEEMinEta_( pset.getParameter< double >( "bosonEEMinEta" ) ),
  bosonEEMaxEta_( pset.getParameter< double >( "bosonEEMaxEta" ) )
{
  bosonPtBins_  = pset.getParameter< std::vector<double> >( "bosonPtBins" );
  bosonEtaBins_ = pset.getParameter< std::vector<double> >( "bosonEtaBins" );
  nJetBins_     = pset.getParameter< std::vector<double> >( "nJetBins" );
  htBins_       = pset.getParameter< std::vector<double> >( "htBins" );
  mhtBins_      = pset.getParameter< std::vector<double> >( "mhtBins" );

  if (studyRecoIso_) {
    recoJetLabel_    = pset.getParameter< edm::InputTag >( "recoJetLabel" );
    recoPhotonLabel_ = pset.getParameter< edm::InputTag >( "recoPhotonLabel");
  }
}


GenStudy::~GenStudy()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void GenStudy::produce(edm::Event& ev, const edm::EventSetup& es)
{
  using namespace edm;
  //read in the gen particles
  edm::Handle<reco::GenParticleCollection> gens;
  ev.getByLabel(genLabel_,gens);

  //get the jets
  edm::Handle<reco::GenJetCollection> genJets;
  ev.getByLabel(genJetLabel_,genJets);

  ///Reco quantities
  //get the photons
  edm::Handle<edm::View<pat::Photon> > recoPhotons;
  edm::Handle<edm::View<pat::Jet> > recoJets;
  if (studyRecoIso_) {
    ev.getByLabel(recoPhotonLabel_,recoPhotons);
    //get the jets
    ev.getByLabel(recoJetLabel_,recoJets);
  }
  
  double eventWeightPU = 1.;
  edm::Handle<double> puWeight;
  if (doPUReweight_) {
    ev.getByLabel(puWeightLabel_,puWeight);
    eventWeightPU = *puWeight;
  }

  //compute event varibles from gen jets
  int nJetsPt30(0), nJetsPt50Eta25(0);
  double htGenJets(0.), mhtGenJets(0.), meffGenJets(0.);
  //double mhtGenJetsNoPhot(0.), meffGenJetsNoPhot(0.);
  double dPhi1(9.0), dPhi2(9.0), dPhi3(9.0);
  //double dPhi1NoPhot(9.0), dPhi2NoPhot(9.0), dPhi3NoPhot(9.0);
  std::vector<const reco::GenJet*> htJetsGen;
  std::vector<const reco::GenJet*> mhtJetsGen;
  reco::MET::LorentzVector mht(0,0,0,0);
  //reco::MET::LorentzVector mhtNoPhot(0,0,0,0);
  reco::GenJetCollection::const_iterator jet = genJets->begin();
  int jetNum = 0;
  for (; jet != genJets->end(); ++jet) {
    if (jet->pt() > 30) {
      ++nJetsPt30;
      mht -= jet->p4();
      //mhtNoPhot -= jet->p4();
      mhtJetsGen.push_back(&((*genJets)[jetNum]));
      if (jet->pt() > 50. && fabs(jet->eta()) < 2.5) {
	++nJetsPt50Eta25;
	htGenJets += jet->pt();
	htJetsGen.push_back(&((*genJets)[jetNum]));
      }
    }
    ++jetNum;
  }

  reco::MET MHT = reco::MET(mht, reco::MET::Point());
  mhtGenJets = MHT.pt();
  meffGenJets = mhtGenJets + htGenJets;

  //remove leading boson from the mht
  //if (removePhot_ && (gens->size() > 0))
  //  mhtNoPhot -= (*gens)[0].p4();

  //reco::MET MHTNoPhot = reco::MET(mhtNoPhot, reco::MET::Point());
  //mhtGenJetsNoPhot = MHTNoPhot.pt();
  ////meffGenJetsNoPhot = mhtGenJetsNoPhot + htGenJets;

  if (mhtJetsGen.size() > 0) {
    dPhi1 = fabs(reco::deltaPhi(mhtJetsGen.at(0)->phi(), MHT.phi()));
    //dPhi1NoPhot = fabs(reco::deltaPhi(mhtJetsGen.at(0)->phi(), MHTNoPhot.phi()));
    if (mhtJetsGen.size() > 1) {
      dPhi2 = fabs(reco::deltaPhi(mhtJetsGen.at(1)->phi(), MHT.phi()));
      //dPhi2NoPhot = fabs(reco::deltaPhi(mhtJetsGen.at(1)->phi(), MHTNoPhot.phi()));
      if (mhtJetsGen.size() > 2) {
	dPhi3 = fabs(reco::deltaPhi(mhtJetsGen.at(2)->phi(), MHT.phi()));
	//dPhi3NoPhot = fabs(reco::deltaPhi(mhtJetsGen.at(2)->phi(), MHTNoPhot.phi()));
      }
    }
  }
  //do the analysis on the gen bosons
  int nBosons(0);
  double bosonPt[2] = {-10.,-10.}, bosonEta[2] = {-10.,-10.}, bosonPhi[2] = {-10.,-10.};
  /*
    std::vector<const reco::GenParticle*> particles;
    reco::GenParticleCollection::const_iterator part = gens->begin();
    for (; part != gens->end(); ++part) {
    if (part->pdgId() == pdgId_)
    if (part->status() == genStatus_) {
    // definition of direct photon 
    bool directPhoton    = ( std::fabs(part->pdgId()) == 22 && part->status() == 3 &&
    (std::fabs(mother->pdgId())<10 || mother->pdgId()==21 )  ) ;
    
    // definition of secondary photon or decay photons
    bool secondaryPhoton = ( std::fabs(part->pdgId()) == 22 && part->status() == 1 && std::fabs(mother->pdgId())>100 );
    
    // definition of fragmentation photon
    bool fragmenPhoton   = ( std::fabs(part->pdgId()) == 22 && part->status() == 1 && (std::fabs(mother->pdgId())<10 || mother->pdgId()==21 )  );
    
    // definition of electron mistags
    // to be done
    }
    particles.push_back(*part);
    }
  */
  nBosons = gens->size();
  if (nBosons < 1) {
    std::cout<<debugString_<<"::Unable to find any gen bosons"<<std::endl;
    return;
  }
  bool passAcc = false;
  bool passRecoIso = false;
  if (nBosons > 0) {
    bosonPt[0]  = (*gens)[0].pt();
    bosonEta[0] = (*gens)[0].eta();
    bosonPhi[0] = (*gens)[0].phi();
    
    if (bosonPt[0] > bosonMinPt_ && (fabs(bosonEta[0]) < bosonEBMaxEta_ || 
				     (fabs(bosonEta[0]) > bosonEEMinEta_ && 
				      fabs(bosonEta[0]) < bosonEEMaxEta_ )
				     )
	)
      passAcc = true;
    if (!studyAcc_)
      passAcc = true;
    
    //Match gen boson to reco/isolated photon
    if (!studyRecoIso_)
      passRecoIso = true;
    else {
      edm::View<pat::Photon>::const_iterator recop = recoPhotons->begin();
      //const pat::Photon *recoMatched;
      //int bestDRPhot = -1;
      double bestDRMin = 999.0;
      int phot = 0;
      for (; recop != recoPhotons->end(); ++recop){
	double dR = reco::deltaR(bosonEta[0],bosonPhi[0],recop->eta(), recop->phi());
	if (dR < bestDRMin) {
	  //bestDRPhot = phot;
	  bestDRMin = dR;
	}
	++phot;
      }
      if (bestDRMin < 0.1) {
	//recoMatched = &((*recoPhotons)[bestDRPhot]);
	passRecoIso = true;
      }
    }
    
    //    //remove the leading boson in the case of photons
    //    double htGenJetsNoPhot(htGenJets), mhtGenJetsNoPhot(0.), meffGenJetsNoPhot(0.);
    //    reco::MET::LorentzVector mhtNoPhot = mht;
    //    mhtNoPhot -= (*gens)[0].p4();
    //    reco::MET MHTNoPhot = reco::MET(mhtNoPhot, reco::MET::Point());
    //    mhtGenJetsNoPhot = MHTNoPhot.pt();
    //    meffGenJetsNoPhot = htGenJetsNoPhot + mhtGenJetsNoPhot;

    if (nBosons > 1) {
      bosonPt[1]  = (*gens)[1].pt();
      bosonEta[1] = (*gens)[1].eta();
      bosonPhi[1] = (*gens)[1].phi();
    }
    
    histos2D_["preRA2_boson1Pt" ]    ->Fill( bosonPt[0],  nJetsPt30, eventWeightPU);
    histos2D_["preRA2_boson1Eta"]    ->Fill( bosonEta[0], nJetsPt30, eventWeightPU);
    histos2D_["preRA2_boson1Phi"]    ->Fill( bosonPhi[0], nJetsPt30, eventWeightPU);
    histos2D_["preRA2_boson1PtBins" ]->Fill( bosonPt[0],  nJetsPt30, eventWeightPU);
    histos2D_["preRA2_boson1EtaBins"]->Fill( bosonEta[0], nJetsPt30, eventWeightPU);

    histos3D_["preRA2_boson1PtvsHTBins" ] ->Fill( htGenJets,  bosonPt[0], nJetsPt30, eventWeightPU);
    histos3D_["preRA2_boson1PtvsMHTBins" ]->Fill( mhtGenJets, bosonPt[0], nJetsPt30, eventWeightPU);
    histos3D_["preRA2_boson1PtvsHT" ]     ->Fill( htGenJets,  bosonPt[0], nJetsPt30, eventWeightPU);
    histos3D_["preRA2_boson1PtvsMHT" ]    ->Fill( mhtGenJets, bosonPt[0], nJetsPt30, eventWeightPU);

    /*    
    if (nBosons > 1) {
      histos2D_["preRA2_boson2Pt" ]    ->Fill( bosonPt[1],  nJetsPt30, eventWeightPU);
      histos2D_["preRA2_boson2Eta"]    ->Fill( bosonEta[1], nJetsPt30, eventWeightPU);
      histos2D_["preRA2_boson2Phi"]    ->Fill( bosonPhi[1], nJetsPt30, eventWeightPU);
      histos2D_["preRA2_boson2PtBins" ]->Fill( bosonPt[1],  nJetsPt30, eventWeightPU);
      histos2D_["preRA2_boson2EtaBins"]->Fill( bosonEta[1], nJetsPt30, eventWeightPU);
      
      histos3D_["preRA2_boson2PtvsHTBins" ] ->Fill( htGenJets,  bosonPt[1], nJetsPt30, eventWeightPU);
      histos3D_["preRA2_boson2PtvsMHTBins" ]->Fill( mhtGenJets, bosonPt[1], nJetsPt30, eventWeightPU);
      histos3D_["preRA2_boson2PtvsHT" ]     ->Fill( htGenJets,  bosonPt[1], nJetsPt30, eventWeightPU);
      histos3D_["preRA2_boson2PtvsMHT" ]    ->Fill( mhtGenJets, bosonPt[1], nJetsPt30, eventWeightPU);
    }
    */
    histos2D_["preRA2_nGenJets_Pt30"     ]->Fill( nJetsPt30,      nJetsPt30, eventWeightPU);
    histos2D_["preRA2_nGenJets_Pt50Eta25"]->Fill( nJetsPt50Eta25, nJetsPt30, eventWeightPU);
    
    histos2D_["preRA2_dPhiJet1MHT" ]->Fill( dPhi1,   nJetsPt30, eventWeightPU);
    histos2D_["preRA2_dPhiJet2MHT" ]->Fill( dPhi2,   nJetsPt30, eventWeightPU);
    histos2D_["preRA2_dPhiJet3MHT" ]->Fill( dPhi3,   nJetsPt30, eventWeightPU);

    histos2D_["preRA2_genHT" ]     ->Fill( htGenJets,   nJetsPt30, eventWeightPU);
    histos2D_["preRA2_genMHT"]     ->Fill( mhtGenJets,  nJetsPt30, eventWeightPU);
    histos2D_["preRA2_genMeff"]    ->Fill( meffGenJets, nJetsPt30, eventWeightPU);
    histos2D_["preRA2_genHTBins" ] ->Fill( htGenJets,   nJetsPt30, eventWeightPU);
    histos2D_["preRA2_genMHTBins"] ->Fill( mhtGenJets,  nJetsPt30, eventWeightPU);
    histos3D_["preRA2_MHTvsHT"]    ->Fill( htGenJets,   mhtGenJets, nJetsPt30, eventWeightPU);
    histos3D_["preRA2_MHTvsHTBins"]->Fill( htGenJets,   mhtGenJets, nJetsPt30, eventWeightPU);
    histos2D_["preRA2_nGenBoson"]  ->Fill( nBosons,     nJetsPt30, eventWeightPU);

    /*
    //No photon in mht
    histos3D_["preRA2_boson1PtvsMHTNoPhotBins" ]->Fill( mhtGenJetsNoPhot, bosonPt[0], nJetsPt30, eventWeightPU);
    histos3D_["preRA2_boson1PtvsMHTNoPhot" ]    ->Fill( mhtGenJetsNoPhot, bosonPt[0], nJetsPt30, eventWeightPU);
    if (nBosons > 1) {
      histos3D_["preRA2_boson2PtvsMHTNoPhotBins" ]->Fill( mhtGenJetsNoPhot, bosonPt[1], nJetsPt30, eventWeightPU);
      histos3D_["preRA2_boson2PtvsMHTNoPhot" ]    ->Fill( mhtGenJetsNoPhot, bosonPt[1], nJetsPt30, eventWeightPU);
    }
    histos2D_["preRA2_dPhiJet1MHTNoPhot" ]->Fill( dPhi1NoPhot,   nJetsPt30, eventWeightPU);
    histos2D_["preRA2_dPhiJet2MHTNoPhot" ]->Fill( dPhi2NoPhot,   nJetsPt30, eventWeightPU);
    histos2D_["preRA2_dPhiJet3MHTNoPhot" ]->Fill( dPhi3NoPhot,   nJetsPt30, eventWeightPU);

    histos2D_["preRA2_genMHTNoPhot"]     ->Fill( mhtGenJetsNoPhot,  nJetsPt30, eventWeightPU);
    histos2D_["preRA2_genMeffNoPhot"]    ->Fill( meffGenJetsNoPhot, nJetsPt30, eventWeightPU);
    histos2D_["preRA2_genMHTNoPhotBins"] ->Fill( mhtGenJetsNoPhot,  nJetsPt30, eventWeightPU);
    histos3D_["preRA2_MHTNoPhotvsHT"]    ->Fill( htGenJets,   mhtGenJetsNoPhot, nJetsPt30, eventWeightPU);
    histos3D_["preRA2_MHTNoPhotvsHTBins"]->Fill( htGenJets,   mhtGenJetsNoPhot, nJetsPt30, eventWeightPU);
    */
    //acceptance cuts
    if ( passAcc) {
      histos2D_["accepted_boson1Pt" ]    ->Fill( bosonPt[0],  nJetsPt30, eventWeightPU);
      histos2D_["accepted_boson1Eta"]    ->Fill( bosonEta[0], nJetsPt30, eventWeightPU);
      histos2D_["accepted_boson1Phi"]    ->Fill( bosonPhi[0], nJetsPt30, eventWeightPU);
      histos2D_["accepted_boson1PtBins" ]->Fill( bosonPt[0],  nJetsPt30, eventWeightPU);
      histos2D_["accepted_boson1EtaBins"]->Fill( bosonEta[0], nJetsPt30, eventWeightPU);

      histos3D_["accepted_boson1PtvsHTBins" ] ->Fill( htGenJets,  bosonPt[0], nJetsPt30, eventWeightPU);
      histos3D_["accepted_boson1PtvsMHTBins" ]->Fill( mhtGenJets, bosonPt[0], nJetsPt30, eventWeightPU);
      histos3D_["accepted_boson1PtvsHT" ]     ->Fill( htGenJets,  bosonPt[0], nJetsPt30, eventWeightPU);
      histos3D_["accepted_boson1PtvsMHT" ]    ->Fill( mhtGenJets, bosonPt[0], nJetsPt30, eventWeightPU);
      /*
      if (nBosons > 1) {
	histos2D_["accepted_boson2Pt" ]    ->Fill( bosonPt[1],  nJetsPt30, eventWeightPU);
	histos2D_["accepted_boson2Eta"]    ->Fill( bosonEta[1], nJetsPt30, eventWeightPU);
	histos2D_["accepted_boson2Phi"]    ->Fill( bosonPhi[1], nJetsPt30, eventWeightPU);
	histos2D_["accepted_boson2PtBins" ]->Fill( bosonPt[1],  nJetsPt30, eventWeightPU);
	histos2D_["accepted_boson2EtaBins"]->Fill( bosonEta[1], nJetsPt30, eventWeightPU);
      
	histos3D_["accepted_boson2PtvsHTBins" ] ->Fill( htGenJets,  bosonPt[1], nJetsPt30, eventWeightPU);
	histos3D_["accepted_boson2PtvsMHTBins" ]->Fill( mhtGenJets, bosonPt[1], nJetsPt30, eventWeightPU);
	histos3D_["accepted_boson2PtvsHT" ]     ->Fill( htGenJets,  bosonPt[1], nJetsPt30, eventWeightPU);
	histos3D_["accepted_boson2PtvsMHT" ]    ->Fill( mhtGenJets, bosonPt[1], nJetsPt30, eventWeightPU);
      }
      */
      histos2D_["accepted_nGenJets_Pt30"     ]->Fill( nJetsPt30,      nJetsPt30, eventWeightPU);
      histos2D_["accepted_nGenJets_Pt50Eta25"]->Fill( nJetsPt50Eta25, nJetsPt30, eventWeightPU);
    
      histos2D_["accepted_dPhiJet1MHT" ]->Fill( dPhi1,   nJetsPt30, eventWeightPU);
      histos2D_["accepted_dPhiJet2MHT" ]->Fill( dPhi2,   nJetsPt30, eventWeightPU);
      histos2D_["accepted_dPhiJet3MHT" ]->Fill( dPhi3,   nJetsPt30, eventWeightPU);

      histos2D_["accepted_genHT" ]     ->Fill( htGenJets,   nJetsPt30, eventWeightPU);
      histos2D_["accepted_genMHT"]     ->Fill( mhtGenJets,  nJetsPt30, eventWeightPU);
      histos2D_["accepted_genMeff"]    ->Fill( meffGenJets, nJetsPt30, eventWeightPU);
      histos2D_["accepted_genHTBins" ] ->Fill( htGenJets,   nJetsPt30, eventWeightPU);
      histos2D_["accepted_genMHTBins"] ->Fill( mhtGenJets,  nJetsPt30, eventWeightPU);
      histos3D_["accepted_MHTvsHT"]    ->Fill( htGenJets,   mhtGenJets, nJetsPt30, eventWeightPU);
      histos3D_["accepted_MHTvsHTBins"]->Fill( htGenJets,   mhtGenJets, nJetsPt30, eventWeightPU);
      histos2D_["accepted_nGenBoson"]  ->Fill( nBosons,     nJetsPt30, eventWeightPU);

      /*  
      //No photon in mht
      histos3D_["accepted_boson1PtvsMHTNoPhotBins" ]->Fill( mhtGenJetsNoPhot, bosonPt[0], nJetsPt30, eventWeightPU);
      histos3D_["accepted_boson1PtvsMHTNoPhot" ]    ->Fill( mhtGenJetsNoPhot, bosonPt[0], nJetsPt30, eventWeightPU);
      if (nBosons > 1) {
	histos3D_["accepted_boson2PtvsMHTNoPhotBins" ]->Fill( mhtGenJetsNoPhot, bosonPt[1], nJetsPt30, eventWeightPU);
	histos3D_["accepted_boson2PtvsMHTNoPhot" ]    ->Fill( mhtGenJetsNoPhot, bosonPt[1], nJetsPt30, eventWeightPU);
      }
      histos2D_["accepted_dPhiJet1MHTNoPhot" ]->Fill( dPhi1NoPhot,   nJetsPt30, eventWeightPU);
      histos2D_["accepted_dPhiJet2MHTNoPhot" ]->Fill( dPhi2NoPhot,   nJetsPt30, eventWeightPU);
      histos2D_["accepted_dPhiJet3MHTNoPhot" ]->Fill( dPhi3NoPhot,   nJetsPt30, eventWeightPU);
      
      histos2D_["accepted_genMHTNoPhot"]     ->Fill( mhtGenJetsNoPhot,  nJetsPt30, eventWeightPU);
      histos2D_["accepted_genMeffNoPhot"]    ->Fill( meffGenJetsNoPhot, nJetsPt30, eventWeightPU);
      histos2D_["accepted_genMHTNoPhotBins"] ->Fill( mhtGenJetsNoPhot,  nJetsPt30, eventWeightPU);
      histos3D_["accepted_MHTNoPhotvsHT"]    ->Fill( htGenJets,   mhtGenJetsNoPhot, nJetsPt30, eventWeightPU);
      histos3D_["accepted_MHTNoPhotvsHTBins"]->Fill( htGenJets,   mhtGenJetsNoPhot, nJetsPt30, eventWeightPU);
      */

      //reco/iso cuts
      if ( passRecoIso) {
	histos2D_["recoiso_boson1Pt" ]    ->Fill( bosonPt[0],  nJetsPt30, eventWeightPU);
	histos2D_["recoiso_boson1Eta"]    ->Fill( bosonEta[0], nJetsPt30, eventWeightPU);
	histos2D_["recoiso_boson1Phi"]    ->Fill( bosonPhi[0], nJetsPt30, eventWeightPU);
	histos2D_["recoiso_boson1PtBins" ]->Fill( bosonPt[0],  nJetsPt30, eventWeightPU);
	histos2D_["recoiso_boson1EtaBins"]->Fill( bosonEta[0], nJetsPt30, eventWeightPU);
	  
	histos3D_["recoiso_boson1PtvsHTBins" ] ->Fill( htGenJets,  bosonPt[0], nJetsPt30, eventWeightPU);
	histos3D_["recoiso_boson1PtvsMHTBins" ]->Fill( mhtGenJets, bosonPt[0], nJetsPt30, eventWeightPU);
	histos3D_["recoiso_boson1PtvsHT" ]     ->Fill( htGenJets,  bosonPt[0], nJetsPt30, eventWeightPU);
	histos3D_["recoiso_boson1PtvsMHT" ]    ->Fill( mhtGenJets, bosonPt[0], nJetsPt30, eventWeightPU);
	/*
	  if (nBosons > 1) {
	  histos2D_["recoiso_boson2Pt" ]    ->Fill( bosonPt[1],  nJetsPt30, eventWeightPU);
	  histos2D_["recoiso_boson2Eta"]    ->Fill( bosonEta[1], nJetsPt30, eventWeightPU);
	  histos2D_["recoiso_boson2Phi"]    ->Fill( bosonPhi[1], nJetsPt30, eventWeightPU);
	  histos2D_["recoiso_boson2PtBins" ]->Fill( bosonPt[1],  nJetsPt30, eventWeightPU);
	  histos2D_["recoiso_boson2EtaBins"]->Fill( bosonEta[1], nJetsPt30, eventWeightPU);
	    
	  histos3D_["recoiso_boson2PtvsHTBins" ] ->Fill( htGenJets,  bosonPt[1], nJetsPt30, eventWeightPU);
	  histos3D_["recoiso_boson2PtvsMHTBins" ]->Fill( mhtGenJets, bosonPt[1], nJetsPt30, eventWeightPU);
	  histos3D_["recoiso_boson2PtvsHT" ]     ->Fill( htGenJets,  bosonPt[1], nJetsPt30, eventWeightPU);
	  histos3D_["recoiso_boson2PtvsMHT" ]    ->Fill( mhtGenJets, bosonPt[1], nJetsPt30, eventWeightPU);
	  }
	*/
	histos2D_["recoiso_nGenJets_Pt30"     ]->Fill( nJetsPt30,      nJetsPt30, eventWeightPU);
	histos2D_["recoiso_nGenJets_Pt50Eta25"]->Fill( nJetsPt50Eta25, nJetsPt30, eventWeightPU);
	  
	histos2D_["recoiso_dPhiJet1MHT" ]->Fill( dPhi1,   nJetsPt30, eventWeightPU);
	histos2D_["recoiso_dPhiJet2MHT" ]->Fill( dPhi2,   nJetsPt30, eventWeightPU);
	histos2D_["recoiso_dPhiJet3MHT" ]->Fill( dPhi3,   nJetsPt30, eventWeightPU);
	  
	histos2D_["recoiso_genHT" ]     ->Fill( htGenJets,   nJetsPt30, eventWeightPU);
	histos2D_["recoiso_genMHT"]     ->Fill( mhtGenJets,  nJetsPt30, eventWeightPU);
	histos2D_["recoiso_genMeff"]    ->Fill( meffGenJets, nJetsPt30, eventWeightPU);
	histos2D_["recoiso_genHTBins" ] ->Fill( htGenJets,   nJetsPt30, eventWeightPU);
	histos2D_["recoiso_genMHTBins"] ->Fill( mhtGenJets,  nJetsPt30, eventWeightPU);
	histos3D_["recoiso_MHTvsHT"]    ->Fill( htGenJets,   mhtGenJets, nJetsPt30, eventWeightPU);
	histos3D_["recoiso_MHTvsHTBins"]->Fill( htGenJets,   mhtGenJets, nJetsPt30, eventWeightPU);
	histos2D_["recoiso_nGenBoson"]  ->Fill( nBosons,     nJetsPt30, eventWeightPU);

	/*	  
	//No photon in mht
	histos3D_["recoiso_boson1PtvsMHTNoPhotBins" ]->Fill( mhtGenJetsNoPhot, bosonPt[0], nJetsPt30, eventWeightPU);
	histos3D_["recoiso_boson1PtvsMHTNoPhot" ]    ->Fill( mhtGenJetsNoPhot, bosonPt[0], nJetsPt30, eventWeightPU);
	if (nBosons > 1) {
	histos3D_["recoiso_boson2PtvsMHTNoPhotBins" ]->Fill( mhtGenJetsNoPhot, bosonPt[1], nJetsPt30, eventWeightPU);
	histos3D_["recoiso_boson2PtvsMHTNoPhot" ]    ->Fill( mhtGenJetsNoPhot, bosonPt[1], nJetsPt30, eventWeightPU);
	}
	histos2D_["recoiso_dPhiJet1MHTNoPhot" ]->Fill( dPhi1NoPhot,   nJetsPt30, eventWeightPU);
	histos2D_["recoiso_dPhiJet2MHTNoPhot" ]->Fill( dPhi2NoPhot,   nJetsPt30, eventWeightPU);
	histos2D_["recoiso_dPhiJet3MHTNoPhot" ]->Fill( dPhi3NoPhot,   nJetsPt30, eventWeightPU);
	  
	histos2D_["recoiso_genMHTNoPhot"]     ->Fill( mhtGenJetsNoPhot,  nJetsPt30, eventWeightPU);
	histos2D_["recoiso_genMeffNoPhot"]    ->Fill( meffGenJetsNoPhot, nJetsPt30, eventWeightPU);
	histos2D_["recoiso_genMHTNoPhotBins"] ->Fill( mhtGenJetsNoPhot,  nJetsPt30, eventWeightPU);
	histos3D_["recoiso_MHTNoPhotvsHT"]    ->Fill( htGenJets,   mhtGenJetsNoPhot, nJetsPt30, eventWeightPU);
	histos3D_["recoiso_MHTNoPhotvsHTBins"]->Fill( htGenJets,   mhtGenJetsNoPhot, nJetsPt30, eventWeightPU);
	*/
	//apply RA2 cuts
	if (htGenJets > minHT_ && mhtGenJets > minMHT_ && nJetsPt30 > 1) {
	  histos2D_["preDPhi_boson1Pt" ]    ->Fill( bosonPt[0],  nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_boson1Eta"]    ->Fill( bosonEta[0], nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_boson1Phi"]    ->Fill( bosonPhi[0], nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_boson1PtBins" ]->Fill( bosonPt[0],  nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_boson1EtaBins"]->Fill( bosonEta[0], nJetsPt30, eventWeightPU);
      
	  histos3D_["preDPhi_boson1PtvsHTBins" ] ->Fill( htGenJets,  bosonPt[0], nJetsPt30, eventWeightPU);
	  histos3D_["preDPhi_boson1PtvsMHTBins" ]->Fill( mhtGenJets, bosonPt[0], nJetsPt30, eventWeightPU);
	  histos3D_["preDPhi_boson1PtvsHT" ]     ->Fill( htGenJets,  bosonPt[0], nJetsPt30, eventWeightPU);
	  histos3D_["preDPhi_boson1PtvsMHT" ]    ->Fill( mhtGenJets, bosonPt[0], nJetsPt30, eventWeightPU);
	  /*
	    if (nBosons > 1) {
	    histos2D_["preDPhi_boson2Pt" ]    ->Fill( bosonPt[1],  nJetsPt30, eventWeightPU);
	    histos2D_["preDPhi_boson2Eta"]    ->Fill( bosonEta[1], nJetsPt30, eventWeightPU);
	    histos2D_["preDPhi_boson2Phi"]    ->Fill( bosonPhi[1], nJetsPt30, eventWeightPU);
	    histos2D_["preDPhi_boson2PtBins" ]->Fill( bosonPt[1],  nJetsPt30, eventWeightPU);
	    histos2D_["preDPhi_boson2EtaBins"]->Fill( bosonEta[1], nJetsPt30, eventWeightPU);
	
	    histos3D_["preDPhi_boson2PtvsHTBins" ] ->Fill( htGenJets,  bosonPt[1], nJetsPt30, eventWeightPU);
	    histos3D_["preDPhi_boson2PtvsMHTBins" ]->Fill( mhtGenJets, bosonPt[1], nJetsPt30, eventWeightPU);
	    histos3D_["preDPhi_boson2PtvsHT" ]     ->Fill( htGenJets,  bosonPt[1], nJetsPt30, eventWeightPU);
	    histos3D_["preDPhi_boson2PtvsMHT" ]    ->Fill( mhtGenJets, bosonPt[1], nJetsPt30, eventWeightPU);
	    }
	  */
	  histos2D_["preDPhi_nGenJets_Pt30"     ]->Fill( nJetsPt30,      nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_nGenJets_Pt50Eta25"]->Fill( nJetsPt50Eta25, nJetsPt30, eventWeightPU);
      
	  histos2D_["preDPhi_dPhiJet1MHT" ]->Fill( dPhi1,   nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_dPhiJet2MHT" ]->Fill( dPhi2,   nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_dPhiJet3MHT" ]->Fill( dPhi3,   nJetsPt30, eventWeightPU);

	  histos2D_["preDPhi_genHT" ]     ->Fill( htGenJets,   nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_genMHT"]     ->Fill( mhtGenJets,  nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_genMeff"]    ->Fill( meffGenJets, nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_genHTBins" ] ->Fill( htGenJets,   nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_genMHTBins"] ->Fill( mhtGenJets,  nJetsPt30, eventWeightPU);
	  histos3D_["preDPhi_MHTvsHT"]    ->Fill( htGenJets,   mhtGenJets, nJetsPt30, eventWeightPU);
	  histos3D_["preDPhi_MHTvsHTBins"]->Fill( htGenJets,   mhtGenJets, nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_nGenBoson"]  ->Fill( nBosons,     nJetsPt30, eventWeightPU);

	  /*
	  //No photon in mht
	  histos3D_["preDPhi_boson1PtvsMHTNoPhotBins" ]->Fill( mhtGenJetsNoPhot, bosonPt[0], nJetsPt30, eventWeightPU);
	  histos3D_["preDPhi_boson1PtvsMHTNoPhot" ]    ->Fill( mhtGenJetsNoPhot, bosonPt[0], nJetsPt30, eventWeightPU);
	  if (nBosons > 1) {
	  histos3D_["preDPhi_boson2PtvsMHTNoPhotBins" ]->Fill( mhtGenJetsNoPhot, bosonPt[1], nJetsPt30, eventWeightPU);
	  histos3D_["preDPhi_boson2PtvsMHTNoPhot" ]    ->Fill( mhtGenJetsNoPhot, bosonPt[1], nJetsPt30, eventWeightPU);
	  }
	  histos2D_["preDPhi_dPhiJet1MHTNoPhot" ]->Fill( dPhi1NoPhot,   nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_dPhiJet2MHTNoPhot" ]->Fill( dPhi2NoPhot,   nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_dPhiJet3MHTNoPhot" ]->Fill( dPhi3NoPhot,   nJetsPt30, eventWeightPU);

	  histos2D_["preDPhi_genMHTNoPhot"]     ->Fill( mhtGenJetsNoPhot,  nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_genMeffNoPhot"]    ->Fill( meffGenJetsNoPhot, nJetsPt30, eventWeightPU);
	  histos2D_["preDPhi_genMHTNoPhotBins"] ->Fill( mhtGenJetsNoPhot,  nJetsPt30, eventWeightPU);
	  histos3D_["preDPhi_MHTNoPhotvsHT"]    ->Fill( htGenJets,   mhtGenJetsNoPhot, nJetsPt30, eventWeightPU);
	  histos3D_["preDPhi_MHTNoPhotvsHTBins"]->Fill( htGenJets,   mhtGenJetsNoPhot, nJetsPt30, eventWeightPU);
	  */
	  //apply DPhi cuts
	  if (dPhi1 > 0.5 && dPhi2 > 0.5 && dPhi3 > 0.3) {
	    histos2D_["boson1Pt" ]    ->Fill( bosonPt[0],  nJetsPt30, eventWeightPU);
	    histos2D_["boson1Eta"]    ->Fill( bosonEta[0], nJetsPt30, eventWeightPU);
	    histos2D_["boson1Phi"]    ->Fill( bosonPhi[0], nJetsPt30, eventWeightPU);
	    histos2D_["boson1PtBins" ]->Fill( bosonPt[0],  nJetsPt30, eventWeightPU);
	    histos2D_["boson1EtaBins"]->Fill( bosonEta[0], nJetsPt30, eventWeightPU);
	    
	    histos3D_["boson1PtvsHTBins" ] ->Fill( htGenJets,  bosonPt[0], nJetsPt30, eventWeightPU);
	    histos3D_["boson1PtvsMHTBins" ]->Fill( mhtGenJets, bosonPt[0], nJetsPt30, eventWeightPU);
	    histos3D_["boson1PtvsHT" ]     ->Fill( htGenJets,  bosonPt[0], nJetsPt30, eventWeightPU);
	    histos3D_["boson1PtvsMHT" ]    ->Fill( mhtGenJets, bosonPt[0], nJetsPt30, eventWeightPU);
	    /*
	      if (nBosons > 1) {
	      histos2D_["boson2Pt" ]    ->Fill( bosonPt[1],  nJetsPt30, eventWeightPU);
	      histos2D_["boson2Eta"]    ->Fill( bosonEta[1], nJetsPt30, eventWeightPU);
	      histos2D_["boson2Phi"]    ->Fill( bosonPhi[1], nJetsPt30, eventWeightPU);
	      histos2D_["boson2PtBins" ]->Fill( bosonPt[1],  nJetsPt30, eventWeightPU);
	      histos2D_["boson2EtaBins"]->Fill( bosonEta[1], nJetsPt30, eventWeightPU);
	      
	      histos3D_["boson2PtvsHTBins" ] ->Fill( htGenJets,  bosonPt[1], nJetsPt30, eventWeightPU);
	      histos3D_["boson2PtvsMHTBins" ]->Fill( mhtGenJets, bosonPt[1], nJetsPt30, eventWeightPU);
	      histos3D_["boson2PtvsHT" ]     ->Fill( htGenJets,  bosonPt[1], nJetsPt30, eventWeightPU);
	      histos3D_["boson2PtvsMHT" ]    ->Fill( mhtGenJets, bosonPt[1], nJetsPt30, eventWeightPU);
	      }
	    */
	    histos2D_["nGenJets_Pt30"     ]->Fill( nJetsPt30,      nJetsPt30, eventWeightPU);
	    histos2D_["nGenJets_Pt50Eta25"]->Fill( nJetsPt50Eta25, nJetsPt30, eventWeightPU);
	    
	    histos2D_["dPhiJet1MHT" ]->Fill( dPhi1,   nJetsPt30, eventWeightPU);
	    histos2D_["dPhiJet2MHT" ]->Fill( dPhi2,   nJetsPt30, eventWeightPU);
	    histos2D_["dPhiJet3MHT" ]->Fill( dPhi3,   nJetsPt30, eventWeightPU);
	    
	    histos2D_["genHT" ]     ->Fill( htGenJets,   nJetsPt30, eventWeightPU);
	    histos2D_["genMHT"]     ->Fill( mhtGenJets,  nJetsPt30, eventWeightPU);
	    histos2D_["genMeff"]    ->Fill( meffGenJets, nJetsPt30, eventWeightPU);
	    histos2D_["genHTBins" ] ->Fill( htGenJets,   nJetsPt30, eventWeightPU);
	    histos2D_["genMHTBins"] ->Fill( mhtGenJets,  nJetsPt30, eventWeightPU);
	    histos3D_["MHTvsHT"]    ->Fill( htGenJets,   mhtGenJets, nJetsPt30, eventWeightPU);
	    histos3D_["MHTvsHTBins"]->Fill( htGenJets,   mhtGenJets, nJetsPt30, eventWeightPU);
	    histos2D_["nGenBoson"]  ->Fill( nBosons,     nJetsPt30, eventWeightPU);
	    
	    /*	      
	    //No photon in mht
	    histos3D_["boson1PtvsMHTNoPhotBins" ]->Fill( mhtGenJetsNoPhot, bosonPt[0], nJetsPt30, eventWeightPU);
	    histos3D_["boson1PtvsMHTNoPhot" ]    ->Fill( mhtGenJetsNoPhot, bosonPt[0], nJetsPt30, eventWeightPU);
	    if (nBosons > 1) {
	    histos3D_["boson2PtvsMHTNoPhotBins" ]->Fill( mhtGenJetsNoPhot, bosonPt[1], nJetsPt30, eventWeightPU);
	    histos3D_["boson2PtvsMHTNoPhot" ]    ->Fill( mhtGenJetsNoPhot, bosonPt[1], nJetsPt30, eventWeightPU);
	    }
	    histos2D_["dPhiJet1MHTNoPhot" ]->Fill( dPhi1NoPhot,   nJetsPt30, eventWeightPU);
	    histos2D_["dPhiJet2MHTNoPhot" ]->Fill( dPhi2NoPhot,   nJetsPt30, eventWeightPU);
	    histos2D_["dPhiJet3MHTNoPhot" ]->Fill( dPhi3NoPhot,   nJetsPt30, eventWeightPU);
	    
	    histos2D_["genMHTNoPhot"]     ->Fill( mhtGenJetsNoPhot,  nJetsPt30, eventWeightPU);
	    histos2D_["genMeffNoPhot"]    ->Fill( meffGenJetsNoPhot, nJetsPt30, eventWeightPU);
	    histos2D_["genMHTNoPhotBins"] ->Fill( mhtGenJetsNoPhot,  nJetsPt30, eventWeightPU);
	    histos3D_["MHTNoPhotvsHT"]    ->Fill( htGenJets,   mhtGenJetsNoPhot, nJetsPt30, eventWeightPU);
	    histos3D_["MHTNoPhotvsHTBins"]->Fill( htGenJets,   mhtGenJetsNoPhot, nJetsPt30, eventWeightPU);
	    */
	  }
	}    
      }
    }
  }
}
// ------------ method called once each job just before starting event loop  ------------
void GenStudy::beginJob()
{
  edm::Service< TFileService > fs;
  
  //char title[128];
  //sprintf(title,"# of events passing HLT_%s_v*",unbiasedTrigger_.c_str());

  /////pre RA2 Cuts
  histos2D_[ "preRA2_nGenBoson" ] = fs->make< TH2D >( "preRA2_numGenBoson", "#boson^{GEN}", 5, -0.5,4.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_nGenBoson" ]->SetXTitle( "#boson^{GEN}" );
  histos2D_[ "preRA2_nGenBoson" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_nGenBoson" ]->Sumw2();

  histos2D_[ "preRA2_nGenJets_Pt30" ] = fs->make< TH2D >( "preRA2_nGenJets_Pt30", "#Jets^{GEN} (p_{T} > 30 GeV, |#eta| < 5.0)", 25, -0.5,24.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_nGenJets_Pt30" ]->SetXTitle( "#Jets^{GEN}" );
  histos2D_[ "preRA2_nGenJets_Pt30" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_nGenJets_Pt30" ]->Sumw2();

  histos2D_[ "preRA2_nGenJets_Pt50Eta25" ] = fs->make< TH2D >( "preRA2_nGenJets_Pt50Eta25", "#Jets^{GEN} (p_{T} > 50 GeV, |#eta| < 2.5)", 25, -0.5,24.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_nGenJets_Pt50Eta25" ]->SetXTitle( "#Jets^{GEN}" );
  histos2D_[ "preRA2_nGenJets_Pt50Eta25" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_nGenJets_Pt50Eta25" ]->Sumw2();

  histos2D_[ "preRA2_dPhiJet1MHT" ] = fs->make< TH2D >( "preRA2_dPhiJet1MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_dPhiJet1MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "preRA2_dPhiJet1MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_dPhiJet1MHT" ]->Sumw2();

  histos2D_[ "preRA2_dPhiJet2MHT" ] = fs->make< TH2D >( "preRA2_dPhiJet2MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_dPhiJet2MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "preRA2_dPhiJet2MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_dPhiJet2MHT" ]->Sumw2();

  histos2D_[ "preRA2_dPhiJet3MHT" ] = fs->make< TH2D >( "preRA2_dPhiJet3MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_dPhiJet3MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "preRA2_dPhiJet3MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_dPhiJet3MHT" ]->Sumw2();

  histos2D_[ "preRA2_genHTBins" ] = fs->make< TH2D >( "preRA2_genHTBins", "Event H^{GEN}_{T}", htBins_.size()-1, htBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_genHTBins" ]->SetXTitle( "Event H^{GEN}_{T}" );
  histos2D_[ "preRA2_genHTBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_genHTBins" ]->Sumw2();

  histos2D_[ "preRA2_genMHTBins" ] = fs->make< TH2D >( "preRA2_genMHTBins", "Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_genMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos2D_[ "preRA2_genMHTBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_genMHTBins" ]->Sumw2();

  histos2D_[ "preRA2_genHT" ] = fs->make< TH2D >( "preRA2_genHT", "Event H^{GEN}_{T}", 100, 0, 5000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_genHT" ]->SetXTitle( "Event H^{GEN}_{T}" );
  histos2D_[ "preRA2_genHT" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_genHT" ]->Sumw2();

  histos2D_[ "preRA2_genMHT" ] = fs->make< TH2D >( "preRA2_genMHT", "Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_genMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos2D_[ "preRA2_genMHT" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_genMHT" ]->Sumw2();

  histos2D_[ "preRA2_genMeff" ] = fs->make< TH2D >( "preRA2_genMeff", "Event M^{GEN}_{eff}", 150, 0.0, 7500.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_genMeff" ]->SetXTitle( "Event M^{GEN}_{eff}" );
  histos2D_[ "preRA2_genMeff" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_genMeff" ]->Sumw2();

  histos2D_[ "preRA2_boson1Pt" ] = fs->make< TH2D >( "preRA2_boson1Pt", "boson^{GEN} p_{T}", 60, 0.0, 600.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_boson1Pt" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "preRA2_boson1Pt" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_boson1Pt" ]->Sumw2();

  histos2D_[ "preRA2_boson1Eta" ] = fs->make< TH2D >( "preRA2_boson1Eta", "boson^{GEN} #eta", 50, -5.0, 5.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_boson1Eta" ]->SetXTitle( "boson #eta" );
  histos2D_[ "preRA2_boson1Eta" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_boson1Eta" ]->Sumw2();

  histos2D_[ "preRA2_boson1PtBins" ] = fs->make< TH2D >( "preRA2_boson1PtBins", "boson^{GEN} p_{T}", bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_boson1PtBins" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "preRA2_boson1PtBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_boson1PtBins" ]->Sumw2();

  histos2D_[ "preRA2_boson1EtaBins" ] = fs->make< TH2D >( "preRA2_boson1EtaBins", "boson^{GEN} #eta", bosonEtaBins_.size()-1, bosonEtaBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_boson1EtaBins" ]->SetXTitle( "boson #eta" );
  histos2D_[ "preRA2_boson1EtaBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_boson1EtaBins" ]->Sumw2();

  histos2D_[ "preRA2_boson1Phi" ] = fs->make< TH2D >( "preRA2_boson1Phi", "boson^{GEN} #phi", 50, -M_PI, M_PI , nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_boson1Phi" ]->SetXTitle( "boson #phi" );
  histos2D_[ "preRA2_boson1Phi" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_boson1Phi" ]->Sumw2();

  /*
  histos2D_[ "preRA2_boson2Pt" ] = fs->make< TH2D >( "preRA2_boson2Pt", "boson^{GEN} p_{T}", 60, 0.0, 600.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_boson2Pt" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "preRA2_boson2Pt" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_boson2Pt" ]->Sumw2();

  histos2D_[ "preRA2_boson2Eta" ] = fs->make< TH2D >( "preRA2_boson2Eta", "boson^{GEN} #eta", 50, -5.0, 5.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_boson2Eta" ]->SetXTitle( "boson #eta" );
  histos2D_[ "preRA2_boson2Eta" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_boson2Eta" ]->Sumw2();

  histos2D_[ "preRA2_boson2PtBins" ] = fs->make< TH2D >( "preRA2_boson2PtBins", "boson^{GEN} p_{T}", bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_boson2PtBins" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "preRA2_boson2PtBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_boson2PtBins" ]->Sumw2();

  histos2D_[ "preRA2_boson2EtaBins" ] = fs->make< TH2D >( "preRA2_boson2EtaBins", "boson^{GEN} #eta", bosonEtaBins_.size()-1, bosonEtaBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_boson2EtaBins" ]->SetXTitle( "boson #eta" );
  histos2D_[ "preRA2_boson2EtaBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_boson2EtaBins" ]->Sumw2();

  histos2D_[ "preRA2_boson2Phi" ] = fs->make< TH2D >( "preRA2_boson2Phi", "boson^{GEN} #phi", 50, -M_PI, M_PI , nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_boson2Phi" ]->SetXTitle( "boson #phi" );
  histos2D_[ "preRA2_boson2Phi" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_boson2Phi" ]->Sumw2();
  */
  //eta vs. pt
  histos3D_[ "preRA2_boson1PtvsHT" ] = fs->make< TH3D >( "preRA2_boson1PtvsHT", "boson^{GEN} p_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preRA2_boson1PtvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preRA2_boson1PtvsHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson1PtvsHT" ]->Sumw2();

  histos3D_[ "preRA2_boson1PtvsMHT" ] = fs->make< TH3D >( "preRA2_boson1PtvsMHT", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preRA2_boson1PtvsMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "preRA2_boson1PtvsMHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson1PtvsMHT" ]->Sumw2();

  /*
  histos3D_[ "preRA2_boson2PtvsHT" ] = fs->make< TH3D >( "preRA2_boson2PtvsHT", "boson^{GEN} p_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preRA2_boson2PtvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preRA2_boson2PtvsHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson2PtvsHT" ]->Sumw2();

  histos3D_[ "preRA2_boson2PtvsMHT" ] = fs->make< TH3D >( "preRA2_boson2PtvsMHT", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preRA2_boson2PtvsMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "preRA2_boson2PtvsMHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson2PtvsMHT" ]->Sumw2();
  */
  histos3D_[ "preRA2_MHTvsHT" ] = fs->make< TH3D >( "preRA2_MHTvsHT", "Event #slashH_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 200, 0.0, 2000.0, 15,-0.5,14.5);
  histos3D_[ "preRA2_MHTvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preRA2_MHTvsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos3D_[ "preRA2_MHTvsHT" ]->Sumw2();

  histos3D_[ "preRA2_boson1PtvsHTBins" ] = fs->make< TH3D >( "preRA2_boson1PtvsHTBins", "boson^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preRA2_boson1PtvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preRA2_boson1PtvsHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson1PtvsHTBins" ]->Sumw2();

  histos3D_[ "preRA2_boson1PtvsMHTBins" ] = fs->make< TH3D >( "preRA2_boson1PtvsMHTBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preRA2_boson1PtvsMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "preRA2_boson1PtvsMHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson1PtvsMHTBins" ]->Sumw2();
  /*
  histos3D_[ "preRA2_boson2PtvsHTBins" ] = fs->make< TH3D >( "preRA2_boson2PtvsHTBins", "boson^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preRA2_boson2PtvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preRA2_boson2PtvsHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson2PtvsHTBins" ]->Sumw2();

  histos3D_[ "preRA2_boson2PtvsMHTBins" ] = fs->make< TH3D >( "preRA2_boson2PtvsMHTBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preRA2_boson2PtvsMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "preRA2_boson2PtvsMHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson2PtvsMHTBins" ]->Sumw2();
  */
  histos3D_[ "preRA2_MHTvsHTBins" ] = fs->make< TH3D >( "preRA2_MHTvsHTBins", "Event #slashH_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preRA2_MHTvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preRA2_MHTvsHTBins" ]->SetYTitle( "Event #slashH_{T}" );
  histos3D_[ "preRA2_MHTvsHTBins" ]->Sumw2();

  /*
  //plots with photon removed mht 
  histos2D_[ "preRA2_dPhiJet1MHTNoPhot" ] = fs->make< TH2D >( "preRA2_dPhiJet1MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_dPhiJet1MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "preRA2_dPhiJet1MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_dPhiJet1MHTNoPhot" ]->Sumw2();

  histos2D_[ "preRA2_dPhiJet2MHTNoPhot" ] = fs->make< TH2D >( "preRA2_dPhiJet2MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_dPhiJet2MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "preRA2_dPhiJet2MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_dPhiJet2MHTNoPhot" ]->Sumw2();

  histos2D_[ "preRA2_dPhiJet3MHTNoPhot" ] = fs->make< TH2D >( "preRA2_dPhiJet3MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_dPhiJet3MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "preRA2_dPhiJet3MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_dPhiJet3MHTNoPhot" ]->Sumw2();

  histos2D_[ "preRA2_genMHTNoPhotBins" ] = fs->make< TH2D >( "preRA2_genMHTNoPhotBins", "Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_genMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos2D_[ "preRA2_genMHTNoPhotBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_genMHTNoPhotBins" ]->Sumw2();

  histos2D_[ "preRA2_genMHTNoPhot" ] = fs->make< TH2D >( "preRA2_genMHTNoPhot", "Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_genMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos2D_[ "preRA2_genMHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_genMHTNoPhot" ]->Sumw2();

  histos2D_[ "preRA2_genMeffNoPhot" ] = fs->make< TH2D >( "preRA2_genMeffNoPhot", "Event M^{GEN,no #gamma}_{eff}", 150, 0.0, 7500.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preRA2_genMeffNoPhot" ]->SetXTitle( "Event M^{GEN,no #gamma}_{eff}" );
  histos2D_[ "preRA2_genMeffNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "preRA2_genMeffNoPhot" ]->Sumw2();

  histos3D_[ "preRA2_boson1PtvsMHTNoPhot" ] = fs->make< TH3D >( "preRA2_boson1PtvsMHTNoPhot", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preRA2_boson1PtvsMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "preRA2_boson1PtvsMHTNoPhot" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson1PtvsMHTNoPhot" ]->Sumw2();

  histos3D_[ "preRA2_boson2PtvsMHTNoPhot" ] = fs->make< TH3D >( "preRA2_boson2PtvsMHTNoPhot", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preRA2_boson2PtvsMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "preRA2_boson2PtvsMHTNoPhot" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson2PtvsMHTNoPhot" ]->Sumw2();

  histos3D_[ "preRA2_MHTNoPhotvsHT" ] = fs->make< TH3D >( "preRA2_MHTNoPhotvsHT", "Event #slashH^{no #gamma}_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 200, 0.0, 2000.0, 15,-0.5,14.5);
  histos3D_[ "preRA2_MHTNoPhotvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preRA2_MHTNoPhotvsHT" ]->SetYTitle( "Event #slashH^{no #gamma}_{T}" );
  histos3D_[ "preRA2_MHTNoPhotvsHT" ]->Sumw2();

  histos3D_[ "preRA2_boson1PtvsMHTNoPhotBins" ] = fs->make< TH3D >( "preRA2_boson1PtvsMHTNoPhotBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preRA2_boson1PtvsMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "preRA2_boson1PtvsMHTNoPhotBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson1PtvsMHTNoPhotBins" ]->Sumw2();

  histos3D_[ "preRA2_boson2PtvsMHTNoPhotBins" ] = fs->make< TH3D >( "preRA2_boson2PtvsMHTNoPhotBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preRA2_boson2PtvsMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "preRA2_boson2PtvsMHTNoPhotBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preRA2_boson2PtvsMHTNoPhotBins" ]->Sumw2();

  histos3D_[ "preRA2_MHTNoPhotvsHTBins" ] = fs->make< TH3D >( "preRA2_MHTNoPhotvsHTBins", "Event #slashH^{no #gamma}_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preRA2_MHTNoPhotvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preRA2_MHTNoPhotvsHTBins" ]->SetYTitle( "Event #slashH^{no #gamma}_{T}" );
  histos3D_[ "preRA2_MHTNoPhotvsHTBins" ]->Sumw2();
  */
  /////pass acceptance Cuts
  histos2D_[ "accepted_nGenBoson" ] = fs->make< TH2D >( "accepted_numGenBoson", "#boson^{GEN}", 5, -0.5,4.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_nGenBoson" ]->SetXTitle( "#boson^{GEN}" );
  histos2D_[ "accepted_nGenBoson" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_nGenBoson" ]->Sumw2();

  histos2D_[ "accepted_nGenJets_Pt30" ] = fs->make< TH2D >( "accepted_nGenJets_Pt30", "#Jets^{GEN} (p_{T} > 30 GeV, |#eta| < 5.0)", 25, -0.5,24.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_nGenJets_Pt30" ]->SetXTitle( "#Jets^{GEN}" );
  histos2D_[ "accepted_nGenJets_Pt30" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_nGenJets_Pt30" ]->Sumw2();

  histos2D_[ "accepted_nGenJets_Pt50Eta25" ] = fs->make< TH2D >( "accepted_nGenJets_Pt50Eta25", "#Jets^{GEN} (p_{T} > 50 GeV, |#eta| < 2.5)", 25, -0.5,24.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_nGenJets_Pt50Eta25" ]->SetXTitle( "#Jets^{GEN}" );
  histos2D_[ "accepted_nGenJets_Pt50Eta25" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_nGenJets_Pt50Eta25" ]->Sumw2();

  histos2D_[ "accepted_dPhiJet1MHT" ] = fs->make< TH2D >( "accepted_dPhiJet1MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_dPhiJet1MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "accepted_dPhiJet1MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_dPhiJet1MHT" ]->Sumw2();

  histos2D_[ "accepted_dPhiJet2MHT" ] = fs->make< TH2D >( "accepted_dPhiJet2MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_dPhiJet2MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "accepted_dPhiJet2MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_dPhiJet2MHT" ]->Sumw2();

  histos2D_[ "accepted_dPhiJet3MHT" ] = fs->make< TH2D >( "accepted_dPhiJet3MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_dPhiJet3MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "accepted_dPhiJet3MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_dPhiJet3MHT" ]->Sumw2();

  histos2D_[ "accepted_genHTBins" ] = fs->make< TH2D >( "accepted_genHTBins", "Event H^{GEN}_{T}", htBins_.size()-1, htBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_genHTBins" ]->SetXTitle( "Event H^{GEN}_{T}" );
  histos2D_[ "accepted_genHTBins" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_genHTBins" ]->Sumw2();

  histos2D_[ "accepted_genMHTBins" ] = fs->make< TH2D >( "accepted_genMHTBins", "Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_genMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos2D_[ "accepted_genMHTBins" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_genMHTBins" ]->Sumw2();

  histos2D_[ "accepted_genHT" ] = fs->make< TH2D >( "accepted_genHT", "Event H^{GEN}_{T}", 100, 0, 5000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_genHT" ]->SetXTitle( "Event H^{GEN}_{T}" );
  histos2D_[ "accepted_genHT" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_genHT" ]->Sumw2();

  histos2D_[ "accepted_genMHT" ] = fs->make< TH2D >( "accepted_genMHT", "Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_genMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos2D_[ "accepted_genMHT" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_genMHT" ]->Sumw2();

  histos2D_[ "accepted_genMeff" ] = fs->make< TH2D >( "accepted_genMeff", "Event M^{GEN}_{eff}", 150, 0.0, 7500.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_genMeff" ]->SetXTitle( "Event M^{GEN}_{eff}" );
  histos2D_[ "accepted_genMeff" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_genMeff" ]->Sumw2();

  histos2D_[ "accepted_boson1Pt" ] = fs->make< TH2D >( "accepted_boson1Pt", "boson^{GEN} p_{T}", 60, 0.0, 600.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_boson1Pt" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "accepted_boson1Pt" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_boson1Pt" ]->Sumw2();

  histos2D_[ "accepted_boson1Eta" ] = fs->make< TH2D >( "accepted_boson1Eta", "boson^{GEN} #eta", 50, -5.0, 5.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_boson1Eta" ]->SetXTitle( "boson #eta" );
  histos2D_[ "accepted_boson1Eta" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_boson1Eta" ]->Sumw2();

  histos2D_[ "accepted_boson1PtBins" ] = fs->make< TH2D >( "accepted_boson1PtBins", "boson^{GEN} p_{T}", bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_boson1PtBins" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "accepted_boson1PtBins" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_boson1PtBins" ]->Sumw2();

  histos2D_[ "accepted_boson1EtaBins" ] = fs->make< TH2D >( "accepted_boson1EtaBins", "boson^{GEN} #eta", bosonEtaBins_.size()-1, bosonEtaBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_boson1EtaBins" ]->SetXTitle( "boson #eta" );
  histos2D_[ "accepted_boson1EtaBins" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_boson1EtaBins" ]->Sumw2();

  histos2D_[ "accepted_boson1Phi" ] = fs->make< TH2D >( "accepted_boson1Phi", "boson^{GEN} #phi", 50, -M_PI, M_PI , nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_boson1Phi" ]->SetXTitle( "boson #phi" );
  histos2D_[ "accepted_boson1Phi" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_boson1Phi" ]->Sumw2();

  /*
  histos2D_[ "accepted_boson2Pt" ] = fs->make< TH2D >( "accepted_boson2Pt", "boson^{GEN} p_{T}", 60, 0.0, 600.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_boson2Pt" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "accepted_boson2Pt" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_boson2Pt" ]->Sumw2();

  histos2D_[ "accepted_boson2Eta" ] = fs->make< TH2D >( "accepted_boson2Eta", "boson^{GEN} #eta", 50, -5.0, 5.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_boson2Eta" ]->SetXTitle( "boson #eta" );
  histos2D_[ "accepted_boson2Eta" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_boson2Eta" ]->Sumw2();

  histos2D_[ "accepted_boson2PtBins" ] = fs->make< TH2D >( "accepted_boson2PtBins", "boson^{GEN} p_{T}", bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_boson2PtBins" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "accepted_boson2PtBins" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_boson2PtBins" ]->Sumw2();

  histos2D_[ "accepted_boson2EtaBins" ] = fs->make< TH2D >( "accepted_boson2EtaBins", "boson^{GEN} #eta", bosonEtaBins_.size()-1, bosonEtaBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_boson2EtaBins" ]->SetXTitle( "boson #eta" );
  histos2D_[ "accepted_boson2EtaBins" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_boson2EtaBins" ]->Sumw2();

  histos2D_[ "accepted_boson2Phi" ] = fs->make< TH2D >( "accepted_boson2Phi", "boson^{GEN} #phi", 50, -M_PI, M_PI , nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_boson2Phi" ]->SetXTitle( "boson #phi" );
  histos2D_[ "accepted_boson2Phi" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_boson2Phi" ]->Sumw2();
  */
  //eta vs. pt
  histos3D_[ "accepted_boson1PtvsHT" ] = fs->make< TH3D >( "accepted_boson1PtvsHT", "boson^{GEN} p_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "accepted_boson1PtvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "accepted_boson1PtvsHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson1PtvsHT" ]->Sumw2();

  histos3D_[ "accepted_boson1PtvsMHT" ] = fs->make< TH3D >( "accepted_boson1PtvsMHT", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "accepted_boson1PtvsMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "accepted_boson1PtvsMHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson1PtvsMHT" ]->Sumw2();

  /*
  histos3D_[ "accepted_boson2PtvsHT" ] = fs->make< TH3D >( "accepted_boson2PtvsHT", "boson^{GEN} p_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "accepted_boson2PtvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "accepted_boson2PtvsHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson2PtvsHT" ]->Sumw2();

  histos3D_[ "accepted_boson2PtvsMHT" ] = fs->make< TH3D >( "accepted_boson2PtvsMHT", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "accepted_boson2PtvsMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "accepted_boson2PtvsMHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson2PtvsMHT" ]->Sumw2();
  */
  histos3D_[ "accepted_MHTvsHT" ] = fs->make< TH3D >( "accepted_MHTvsHT", "Event #slashH_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 200, 0.0, 2000.0, 15,-0.5,14.5);
  histos3D_[ "accepted_MHTvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "accepted_MHTvsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos3D_[ "accepted_MHTvsHT" ]->Sumw2();
  
  histos3D_[ "accepted_boson1PtvsHTBins" ] = fs->make< TH3D >( "accepted_boson1PtvsHTBins", "boson^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "accepted_boson1PtvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "accepted_boson1PtvsHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson1PtvsHTBins" ]->Sumw2();
  
  histos3D_[ "accepted_boson1PtvsMHTBins" ] = fs->make< TH3D >( "accepted_boson1PtvsMHTBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "accepted_boson1PtvsMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "accepted_boson1PtvsMHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson1PtvsMHTBins" ]->Sumw2();

  /*
  histos3D_[ "accepted_boson2PtvsHTBins" ] = fs->make< TH3D >( "accepted_boson2PtvsHTBins", "boson^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "accepted_boson2PtvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "accepted_boson2PtvsHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson2PtvsHTBins" ]->Sumw2();
  
  histos3D_[ "accepted_boson2PtvsMHTBins" ] = fs->make< TH3D >( "accepted_boson2PtvsMHTBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "accepted_boson2PtvsMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "accepted_boson2PtvsMHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson2PtvsMHTBins" ]->Sumw2();
  */

  histos3D_[ "accepted_MHTvsHTBins" ] = fs->make< TH3D >( "accepted_MHTvsHTBins", "Event #slashH_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "accepted_MHTvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "accepted_MHTvsHTBins" ]->SetYTitle( "Event #slashH_{T}" );
  histos3D_[ "accepted_MHTvsHTBins" ]->Sumw2();
  /*
  //plots with photon removed mht 
  histos2D_[ "accepted_dPhiJet1MHTNoPhot" ] = fs->make< TH2D >( "accepted_dPhiJet1MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_dPhiJet1MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "accepted_dPhiJet1MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_dPhiJet1MHTNoPhot" ]->Sumw2();

  histos2D_[ "accepted_dPhiJet2MHTNoPhot" ] = fs->make< TH2D >( "accepted_dPhiJet2MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_dPhiJet2MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "accepted_dPhiJet2MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_dPhiJet2MHTNoPhot" ]->Sumw2();

  histos2D_[ "accepted_dPhiJet3MHTNoPhot" ] = fs->make< TH2D >( "accepted_dPhiJet3MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_dPhiJet3MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "accepted_dPhiJet3MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_dPhiJet3MHTNoPhot" ]->Sumw2();

  histos2D_[ "accepted_genMHTNoPhotBins" ] = fs->make< TH2D >( "accepted_genMHTNoPhotBins", "Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_genMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos2D_[ "accepted_genMHTNoPhotBins" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_genMHTNoPhotBins" ]->Sumw2();

  histos2D_[ "accepted_genMHTNoPhot" ] = fs->make< TH2D >( "accepted_genMHTNoPhot", "Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_genMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos2D_[ "accepted_genMHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_genMHTNoPhot" ]->Sumw2();

  histos2D_[ "accepted_genMeffNoPhot" ] = fs->make< TH2D >( "accepted_genMeffNoPhot", "Event M^{GEN,no #gamma}_{eff}", 150, 0.0, 7500.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "accepted_genMeffNoPhot" ]->SetXTitle( "Event M^{GEN,no #gamma}_{eff}" );
  histos2D_[ "accepted_genMeffNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "accepted_genMeffNoPhot" ]->Sumw2();

  histos3D_[ "accepted_boson1PtvsMHTNoPhot" ] = fs->make< TH3D >( "accepted_boson1PtvsMHTNoPhot", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "accepted_boson1PtvsMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "accepted_boson1PtvsMHTNoPhot" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson1PtvsMHTNoPhot" ]->Sumw2();

  histos3D_[ "accepted_boson2PtvsMHTNoPhot" ] = fs->make< TH3D >( "accepted_boson2PtvsMHTNoPhot", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "accepted_boson2PtvsMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "accepted_boson2PtvsMHTNoPhot" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson2PtvsMHTNoPhot" ]->Sumw2();

  histos3D_[ "accepted_MHTNoPhotvsHT" ] = fs->make< TH3D >( "accepted_MHTNoPhotvsHT", "Event #slashH^{no #gamma}_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 200, 0.0, 2000.0, 15,-0.5,14.5);
  histos3D_[ "accepted_MHTNoPhotvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "accepted_MHTNoPhotvsHT" ]->SetYTitle( "Event #slashH^{no #gamma}_{T}" );
  histos3D_[ "accepted_MHTNoPhotvsHT" ]->Sumw2();

  histos3D_[ "accepted_boson1PtvsMHTNoPhotBins" ] = fs->make< TH3D >( "accepted_boson1PtvsMHTNoPhotBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "accepted_boson1PtvsMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "accepted_boson1PtvsMHTNoPhotBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson1PtvsMHTNoPhotBins" ]->Sumw2();

  histos3D_[ "accepted_boson2PtvsMHTNoPhotBins" ] = fs->make< TH3D >( "accepted_boson2PtvsMHTNoPhotBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "accepted_boson2PtvsMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "accepted_boson2PtvsMHTNoPhotBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "accepted_boson2PtvsMHTNoPhotBins" ]->Sumw2();

  histos3D_[ "accepted_MHTNoPhotvsHTBins" ] = fs->make< TH3D >( "accepted_MHTNoPhotvsHTBins", "Event #slashH^{no #gamma}_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "accepted_MHTNoPhotvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "accepted_MHTNoPhotvsHTBins" ]->SetYTitle( "Event #slashH^{no #gamma}_{T}" );
  histos3D_[ "accepted_MHTNoPhotvsHTBins" ]->Sumw2();
  */

  /////pass reco/iso Cuts
  histos2D_[ "recoiso_nGenBoson" ] = fs->make< TH2D >( "recoiso_numGenBoson", "#boson^{GEN}", 5, -0.5,4.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_nGenBoson" ]->SetXTitle( "#boson^{GEN}" );
  histos2D_[ "recoiso_nGenBoson" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_nGenBoson" ]->Sumw2();

  histos2D_[ "recoiso_nGenJets_Pt30" ] = fs->make< TH2D >( "recoiso_nGenJets_Pt30", "#Jets^{GEN} (p_{T} > 30 GeV, |#eta| < 5.0)", 25, -0.5,24.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_nGenJets_Pt30" ]->SetXTitle( "#Jets^{GEN}" );
  histos2D_[ "recoiso_nGenJets_Pt30" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_nGenJets_Pt30" ]->Sumw2();

  histos2D_[ "recoiso_nGenJets_Pt50Eta25" ] = fs->make< TH2D >( "recoiso_nGenJets_Pt50Eta25", "#Jets^{GEN} (p_{T} > 50 GeV, |#eta| < 2.5)", 25, -0.5,24.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_nGenJets_Pt50Eta25" ]->SetXTitle( "#Jets^{GEN}" );
  histos2D_[ "recoiso_nGenJets_Pt50Eta25" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_nGenJets_Pt50Eta25" ]->Sumw2();

  histos2D_[ "recoiso_dPhiJet1MHT" ] = fs->make< TH2D >( "recoiso_dPhiJet1MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_dPhiJet1MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "recoiso_dPhiJet1MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_dPhiJet1MHT" ]->Sumw2();

  histos2D_[ "recoiso_dPhiJet2MHT" ] = fs->make< TH2D >( "recoiso_dPhiJet2MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_dPhiJet2MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "recoiso_dPhiJet2MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_dPhiJet2MHT" ]->Sumw2();

  histos2D_[ "recoiso_dPhiJet3MHT" ] = fs->make< TH2D >( "recoiso_dPhiJet3MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_dPhiJet3MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "recoiso_dPhiJet3MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_dPhiJet3MHT" ]->Sumw2();

  histos2D_[ "recoiso_genHTBins" ] = fs->make< TH2D >( "recoiso_genHTBins", "Event H^{GEN}_{T}", htBins_.size()-1, htBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_genHTBins" ]->SetXTitle( "Event H^{GEN}_{T}" );
  histos2D_[ "recoiso_genHTBins" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_genHTBins" ]->Sumw2();

  histos2D_[ "recoiso_genMHTBins" ] = fs->make< TH2D >( "recoiso_genMHTBins", "Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_genMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos2D_[ "recoiso_genMHTBins" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_genMHTBins" ]->Sumw2();

  histos2D_[ "recoiso_genHT" ] = fs->make< TH2D >( "recoiso_genHT", "Event H^{GEN}_{T}", 100, 0, 5000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_genHT" ]->SetXTitle( "Event H^{GEN}_{T}" );
  histos2D_[ "recoiso_genHT" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_genHT" ]->Sumw2();

  histos2D_[ "recoiso_genMHT" ] = fs->make< TH2D >( "recoiso_genMHT", "Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_genMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos2D_[ "recoiso_genMHT" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_genMHT" ]->Sumw2();

  histos2D_[ "recoiso_genMeff" ] = fs->make< TH2D >( "recoiso_genMeff", "Event M^{GEN}_{eff}", 150, 0.0, 7500.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_genMeff" ]->SetXTitle( "Event M^{GEN}_{eff}" );
  histos2D_[ "recoiso_genMeff" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_genMeff" ]->Sumw2();

  histos2D_[ "recoiso_boson1Pt" ] = fs->make< TH2D >( "recoiso_boson1Pt", "boson^{GEN} p_{T}", 60, 0.0, 600.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_boson1Pt" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "recoiso_boson1Pt" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_boson1Pt" ]->Sumw2();

  histos2D_[ "recoiso_boson1Eta" ] = fs->make< TH2D >( "recoiso_boson1Eta", "boson^{GEN} #eta", 50, -5.0, 5.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_boson1Eta" ]->SetXTitle( "boson #eta" );
  histos2D_[ "recoiso_boson1Eta" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_boson1Eta" ]->Sumw2();

  histos2D_[ "recoiso_boson1PtBins" ] = fs->make< TH2D >( "recoiso_boson1PtBins", "boson^{GEN} p_{T}", bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_boson1PtBins" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "recoiso_boson1PtBins" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_boson1PtBins" ]->Sumw2();

  histos2D_[ "recoiso_boson1EtaBins" ] = fs->make< TH2D >( "recoiso_boson1EtaBins", "boson^{GEN} #eta", bosonEtaBins_.size()-1, bosonEtaBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_boson1EtaBins" ]->SetXTitle( "boson #eta" );
  histos2D_[ "recoiso_boson1EtaBins" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_boson1EtaBins" ]->Sumw2();

  histos2D_[ "recoiso_boson1Phi" ] = fs->make< TH2D >( "recoiso_boson1Phi", "boson^{GEN} #phi", 50, -M_PI, M_PI , nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_boson1Phi" ]->SetXTitle( "boson #phi" );
  histos2D_[ "recoiso_boson1Phi" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_boson1Phi" ]->Sumw2();

  /*
  histos2D_[ "recoiso_boson2Pt" ] = fs->make< TH2D >( "recoiso_boson2Pt", "boson^{GEN} p_{T}", 60, 0.0, 600.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_boson2Pt" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "recoiso_boson2Pt" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_boson2Pt" ]->Sumw2();

  histos2D_[ "recoiso_boson2Eta" ] = fs->make< TH2D >( "recoiso_boson2Eta", "boson^{GEN} #eta", 50, -5.0, 5.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_boson2Eta" ]->SetXTitle( "boson #eta" );
  histos2D_[ "recoiso_boson2Eta" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_boson2Eta" ]->Sumw2();

  histos2D_[ "recoiso_boson2PtBins" ] = fs->make< TH2D >( "recoiso_boson2PtBins", "boson^{GEN} p_{T}", bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_boson2PtBins" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "recoiso_boson2PtBins" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_boson2PtBins" ]->Sumw2();

  histos2D_[ "recoiso_boson2EtaBins" ] = fs->make< TH2D >( "recoiso_boson2EtaBins", "boson^{GEN} #eta", bosonEtaBins_.size()-1, bosonEtaBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_boson2EtaBins" ]->SetXTitle( "boson #eta" );
  histos2D_[ "recoiso_boson2EtaBins" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_boson2EtaBins" ]->Sumw2();

  histos2D_[ "recoiso_boson2Phi" ] = fs->make< TH2D >( "recoiso_boson2Phi", "boson^{GEN} #phi", 50, -M_PI, M_PI , nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_boson2Phi" ]->SetXTitle( "boson #phi" );
  histos2D_[ "recoiso_boson2Phi" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_boson2Phi" ]->Sumw2();
  */
  //eta vs. pt
  histos3D_[ "recoiso_boson1PtvsHT" ] = fs->make< TH3D >( "recoiso_boson1PtvsHT", "boson^{GEN} p_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "recoiso_boson1PtvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "recoiso_boson1PtvsHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson1PtvsHT" ]->Sumw2();

  histos3D_[ "recoiso_boson1PtvsMHT" ] = fs->make< TH3D >( "recoiso_boson1PtvsMHT", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "recoiso_boson1PtvsMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "recoiso_boson1PtvsMHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson1PtvsMHT" ]->Sumw2();

  /*
  histos3D_[ "recoiso_boson2PtvsHT" ] = fs->make< TH3D >( "recoiso_boson2PtvsHT", "boson^{GEN} p_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "recoiso_boson2PtvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "recoiso_boson2PtvsHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson2PtvsHT" ]->Sumw2();

  histos3D_[ "recoiso_boson2PtvsMHT" ] = fs->make< TH3D >( "recoiso_boson2PtvsMHT", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "recoiso_boson2PtvsMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "recoiso_boson2PtvsMHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson2PtvsMHT" ]->Sumw2();
  */
  histos3D_[ "recoiso_MHTvsHT" ] = fs->make< TH3D >( "recoiso_MHTvsHT", "Event #slashH_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 200, 0.0, 2000.0, 15,-0.5,14.5);
  histos3D_[ "recoiso_MHTvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "recoiso_MHTvsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos3D_[ "recoiso_MHTvsHT" ]->Sumw2();

  histos3D_[ "recoiso_boson1PtvsHTBins" ] = fs->make< TH3D >( "recoiso_boson1PtvsHTBins", "boson^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "recoiso_boson1PtvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "recoiso_boson1PtvsHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson1PtvsHTBins" ]->Sumw2();

  histos3D_[ "recoiso_boson1PtvsMHTBins" ] = fs->make< TH3D >( "recoiso_boson1PtvsMHTBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "recoiso_boson1PtvsMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "recoiso_boson1PtvsMHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson1PtvsMHTBins" ]->Sumw2();

  /*
  histos3D_[ "recoiso_boson2PtvsHTBins" ] = fs->make< TH3D >( "recoiso_boson2PtvsHTBins", "boson^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "recoiso_boson2PtvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "recoiso_boson2PtvsHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson2PtvsHTBins" ]->Sumw2();

  histos3D_[ "recoiso_boson2PtvsMHTBins" ] = fs->make< TH3D >( "recoiso_boson2PtvsMHTBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "recoiso_boson2PtvsMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "recoiso_boson2PtvsMHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson2PtvsMHTBins" ]->Sumw2();
  */
  histos3D_[ "recoiso_MHTvsHTBins" ] = fs->make< TH3D >( "recoiso_MHTvsHTBins", "Event #slashH_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "recoiso_MHTvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "recoiso_MHTvsHTBins" ]->SetYTitle( "Event #slashH_{T}" );
  histos3D_[ "recoiso_MHTvsHTBins" ]->Sumw2();

  /*
  //plots with photon removed mht 
  histos2D_[ "recoiso_dPhiJet1MHTNoPhot" ] = fs->make< TH2D >( "recoiso_dPhiJet1MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_dPhiJet1MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "recoiso_dPhiJet1MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_dPhiJet1MHTNoPhot" ]->Sumw2();

  histos2D_[ "recoiso_dPhiJet2MHTNoPhot" ] = fs->make< TH2D >( "recoiso_dPhiJet2MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_dPhiJet2MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "recoiso_dPhiJet2MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_dPhiJet2MHTNoPhot" ]->Sumw2();

  histos2D_[ "recoiso_dPhiJet3MHTNoPhot" ] = fs->make< TH2D >( "recoiso_dPhiJet3MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_dPhiJet3MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "recoiso_dPhiJet3MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_dPhiJet3MHTNoPhot" ]->Sumw2();

  histos2D_[ "recoiso_genMHTNoPhotBins" ] = fs->make< TH2D >( "recoiso_genMHTNoPhotBins", "Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_genMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos2D_[ "recoiso_genMHTNoPhotBins" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_genMHTNoPhotBins" ]->Sumw2();

  histos2D_[ "recoiso_genMHTNoPhot" ] = fs->make< TH2D >( "recoiso_genMHTNoPhot", "Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_genMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos2D_[ "recoiso_genMHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_genMHTNoPhot" ]->Sumw2();

  histos2D_[ "recoiso_genMeffNoPhot" ] = fs->make< TH2D >( "recoiso_genMeffNoPhot", "Event M^{GEN,no #gamma}_{eff}", 150, 0.0, 7500.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "recoiso_genMeffNoPhot" ]->SetXTitle( "Event M^{GEN,no #gamma}_{eff}" );
  histos2D_[ "recoiso_genMeffNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "recoiso_genMeffNoPhot" ]->Sumw2();

  histos3D_[ "recoiso_boson1PtvsMHTNoPhot" ] = fs->make< TH3D >( "recoiso_boson1PtvsMHTNoPhot", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "recoiso_boson1PtvsMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "recoiso_boson1PtvsMHTNoPhot" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson1PtvsMHTNoPhot" ]->Sumw2();

  histos3D_[ "recoiso_boson2PtvsMHTNoPhot" ] = fs->make< TH3D >( "recoiso_boson2PtvsMHTNoPhot", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "recoiso_boson2PtvsMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "recoiso_boson2PtvsMHTNoPhot" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson2PtvsMHTNoPhot" ]->Sumw2();

  histos3D_[ "recoiso_MHTNoPhotvsHT" ] = fs->make< TH3D >( "recoiso_MHTNoPhotvsHT", "Event #slashH^{no #gamma}_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 200, 0.0, 2000.0, 15,-0.5,14.5);
  histos3D_[ "recoiso_MHTNoPhotvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "recoiso_MHTNoPhotvsHT" ]->SetYTitle( "Event #slashH^{no #gamma}_{T}" );
  histos3D_[ "recoiso_MHTNoPhotvsHT" ]->Sumw2();

  histos3D_[ "recoiso_boson1PtvsMHTNoPhotBins" ] = fs->make< TH3D >( "recoiso_boson1PtvsMHTNoPhotBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "recoiso_boson1PtvsMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "recoiso_boson1PtvsMHTNoPhotBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson1PtvsMHTNoPhotBins" ]->Sumw2();

  histos3D_[ "recoiso_boson2PtvsMHTNoPhotBins" ] = fs->make< TH3D >( "recoiso_boson2PtvsMHTNoPhotBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "recoiso_boson2PtvsMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "recoiso_boson2PtvsMHTNoPhotBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "recoiso_boson2PtvsMHTNoPhotBins" ]->Sumw2();

  histos3D_[ "recoiso_MHTNoPhotvsHTBins" ] = fs->make< TH3D >( "recoiso_MHTNoPhotvsHTBins", "Event #slashH^{no #gamma}_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "recoiso_MHTNoPhotvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "recoiso_MHTNoPhotvsHTBins" ]->SetYTitle( "Event #slashH^{no #gamma}_{T}" );
  histos3D_[ "recoiso_MHTNoPhotvsHTBins" ]->Sumw2();
  */

  //////pre DPhi Cuts
  histos2D_[ "preDPhi_nGenBoson" ] = fs->make< TH2D >( "preDPhi_numGenBoson", "#boson^{GEN}", 5, -0.5,4.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_nGenBoson" ]->SetXTitle( "#boson^{GEN}" );
  histos2D_[ "preDPhi_nGenBoson" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_nGenBoson" ]->Sumw2();

  histos2D_[ "preDPhi_nGenJets_Pt30" ] = fs->make< TH2D >( "preDPhi_nGenJets_Pt30", "#Jets^{GEN} (p_{T} > 30 GeV, |#eta| < 5.0)", 25, -0.5,24.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_nGenJets_Pt30" ]->SetXTitle( "#Jets^{GEN}" );
  histos2D_[ "preDPhi_nGenJets_Pt30" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_nGenJets_Pt30" ]->Sumw2();

  histos2D_[ "preDPhi_nGenJets_Pt50Eta25" ] = fs->make< TH2D >( "preDPhi_nGenJets_Pt50Eta25", "#Jets^{GEN} (p_{T} > 50 GeV, |#eta| < 2.5)", 25, -0.5,24.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_nGenJets_Pt50Eta25" ]->SetXTitle( "#Jets^{GEN}" );
  histos2D_[ "preDPhi_nGenJets_Pt50Eta25" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_nGenJets_Pt50Eta25" ]->Sumw2();

  histos2D_[ "preDPhi_dPhiJet1MHT" ] = fs->make< TH2D >( "preDPhi_dPhiJet1MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_dPhiJet1MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "preDPhi_dPhiJet1MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_dPhiJet1MHT" ]->Sumw2();

  histos2D_[ "preDPhi_dPhiJet2MHT" ] = fs->make< TH2D >( "preDPhi_dPhiJet2MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_dPhiJet2MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "preDPhi_dPhiJet2MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_dPhiJet2MHT" ]->Sumw2();

  histos2D_[ "preDPhi_dPhiJet3MHT" ] = fs->make< TH2D >( "preDPhi_dPhiJet3MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_dPhiJet3MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "preDPhi_dPhiJet3MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_dPhiJet3MHT" ]->Sumw2();

  histos2D_[ "preDPhi_genHTBins" ] = fs->make< TH2D >( "preDPhi_genHTBins", "Event H^{GEN}_{T}", htBins_.size()-1, htBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_genHTBins" ]->SetXTitle( "Event H^{GEN}_{T}" );
  histos2D_[ "preDPhi_genHTBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_genHTBins" ]->Sumw2();

  histos2D_[ "preDPhi_genMHTBins" ] = fs->make< TH2D >( "preDPhi_genMHTBins", "Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_genMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos2D_[ "preDPhi_genMHTBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_genMHTBins" ]->Sumw2();

  histos2D_[ "preDPhi_genHT" ] = fs->make< TH2D >( "preDPhi_genHT", "Event H^{GEN}_{T}", 100, 0, 5000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_genHT" ]->SetXTitle( "Event H^{GEN}_{T}" );
  histos2D_[ "preDPhi_genHT" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_genHT" ]->Sumw2();

  histos2D_[ "preDPhi_genMHT" ] = fs->make< TH2D >( "preDPhi_genMHT", "Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_genMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos2D_[ "preDPhi_genMHT" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_genMHT" ]->Sumw2();

  histos2D_[ "preDPhi_genMeff" ] = fs->make< TH2D >( "preDPhi_genMeff", "Event M^{GEN}_{eff}", 150, 0.0, 7500.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_genMeff" ]->SetXTitle( "Event M^{GEN}_{eff}" );
  histos2D_[ "preDPhi_genMeff" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_genMeff" ]->Sumw2();

  histos2D_[ "preDPhi_boson1Pt" ] = fs->make< TH2D >( "preDPhi_boson1Pt", "boson^{GEN} p_{T}", 60, 0.0, 600.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_boson1Pt" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "preDPhi_boson1Pt" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_boson1Pt" ]->Sumw2();

  histos2D_[ "preDPhi_boson1Eta" ] = fs->make< TH2D >( "preDPhi_boson1Eta", "boson^{GEN} #eta", 50, -5.0, 5.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_boson1Eta" ]->SetXTitle( "boson #eta" );
  histos2D_[ "preDPhi_boson1Eta" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_boson1Eta" ]->Sumw2();

  histos2D_[ "preDPhi_boson1PtBins" ] = fs->make< TH2D >( "preDPhi_boson1PtBins", "boson^{GEN} p_{T}", bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_boson1PtBins" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "preDPhi_boson1PtBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_boson1PtBins" ]->Sumw2();

  histos2D_[ "preDPhi_boson1EtaBins" ] = fs->make< TH2D >( "preDPhi_boson1EtaBins", "boson^{GEN} #eta", bosonEtaBins_.size()-1, bosonEtaBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_boson1EtaBins" ]->SetXTitle( "boson #eta" );
  histos2D_[ "preDPhi_boson1EtaBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_boson1EtaBins" ]->Sumw2();

  histos2D_[ "preDPhi_boson1Phi" ] = fs->make< TH2D >( "preDPhi_boson1Phi", "boson^{GEN} #phi", 50, -M_PI, M_PI , nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_boson1Phi" ]->SetXTitle( "boson #phi" );
  histos2D_[ "preDPhi_boson1Phi" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_boson1Phi" ]->Sumw2();

  /*
  histos2D_[ "preDPhi_boson2Pt" ] = fs->make< TH2D >( "preDPhi_boson2Pt", "boson^{GEN} p_{T}", 60, 0.0, 600.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_boson2Pt" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "preDPhi_boson2Pt" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_boson2Pt" ]->Sumw2();

  histos2D_[ "preDPhi_boson2Eta" ] = fs->make< TH2D >( "preDPhi_boson2Eta", "boson^{GEN} #eta", 50, -5.0, 5.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_boson2Eta" ]->SetXTitle( "boson #eta" );
  histos2D_[ "preDPhi_boson2Eta" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_boson2Eta" ]->Sumw2();

  histos2D_[ "preDPhi_boson2PtBins" ] = fs->make< TH2D >( "preDPhi_boson2PtBins", "boson^{GEN} p_{T}", bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_boson2PtBins" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "preDPhi_boson2PtBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_boson2PtBins" ]->Sumw2();

  histos2D_[ "preDPhi_boson2EtaBins" ] = fs->make< TH2D >( "preDPhi_boson2EtaBins", "boson^{GEN} #eta", bosonEtaBins_.size()-1, bosonEtaBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_boson2EtaBins" ]->SetXTitle( "boson #eta" );
  histos2D_[ "preDPhi_boson2EtaBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_boson2EtaBins" ]->Sumw2();

  histos2D_[ "preDPhi_boson2Phi" ] = fs->make< TH2D >( "preDPhi_boson2Phi", "boson^{GEN} #phi", 50, -M_PI, M_PI , nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_boson2Phi" ]->SetXTitle( "boson #phi" );
  histos2D_[ "preDPhi_boson2Phi" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_boson2Phi" ]->Sumw2();
  */
  //eta vs. pt
  histos3D_[ "preDPhi_boson1PtvsHT" ] = fs->make< TH3D >( "preDPhi_boson1PtvsHT", "boson^{GEN} p_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preDPhi_boson1PtvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preDPhi_boson1PtvsHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson1PtvsHT" ]->Sumw2();

  histos3D_[ "preDPhi_boson1PtvsMHT" ] = fs->make< TH3D >( "preDPhi_boson1PtvsMHT", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preDPhi_boson1PtvsMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "preDPhi_boson1PtvsMHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson1PtvsMHT" ]->Sumw2();
  /*
  histos3D_[ "preDPhi_boson2PtvsHT" ] = fs->make< TH3D >( "preDPhi_boson2PtvsHT", "boson^{GEN} p_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preDPhi_boson2PtvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preDPhi_boson2PtvsHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson2PtvsHT" ]->Sumw2();

  histos3D_[ "preDPhi_boson2PtvsMHT" ] = fs->make< TH3D >( "preDPhi_boson2PtvsMHT", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preDPhi_boson2PtvsMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "preDPhi_boson2PtvsMHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson2PtvsMHT" ]->Sumw2();
  */
  histos3D_[ "preDPhi_MHTvsHT" ] = fs->make< TH3D >( "preDPhi_MHTvsHT", "Event #slashH_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 200, 0.0, 2000.0, 15,-0.5,14.5);
  histos3D_[ "preDPhi_MHTvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preDPhi_MHTvsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos3D_[ "preDPhi_MHTvsHT" ]->Sumw2();

  histos3D_[ "preDPhi_boson1PtvsHTBins" ] = fs->make< TH3D >( "preDPhi_boson1PtvsHTBins", "boson^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preDPhi_boson1PtvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preDPhi_boson1PtvsHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson1PtvsHTBins" ]->Sumw2();

  histos3D_[ "preDPhi_boson1PtvsMHTBins" ] = fs->make< TH3D >( "preDPhi_boson1PtvsMHTBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preDPhi_boson1PtvsMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "preDPhi_boson1PtvsMHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson1PtvsMHTBins" ]->Sumw2();

  /*
  histos3D_[ "preDPhi_boson2PtvsHTBins" ] = fs->make< TH3D >( "preDPhi_boson2PtvsHTBins", "boson^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preDPhi_boson2PtvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preDPhi_boson2PtvsHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson2PtvsHTBins" ]->Sumw2();

  histos3D_[ "preDPhi_boson2PtvsMHTBins" ] = fs->make< TH3D >( "preDPhi_boson2PtvsMHTBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preDPhi_boson2PtvsMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "preDPhi_boson2PtvsMHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson2PtvsMHTBins" ]->Sumw2();
  */
  histos3D_[ "preDPhi_MHTvsHTBins" ] = fs->make< TH3D >( "preDPhi_MHTvsHTBins", "Event #slashH_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preDPhi_MHTvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preDPhi_MHTvsHTBins" ]->SetYTitle( "Event #slashH_{T}" );
  histos3D_[ "preDPhi_MHTvsHTBins" ]->Sumw2();

  /*
  //plots with photon removed mht 
  histos2D_[ "preDPhi_dPhiJet1MHTNoPhot" ] = fs->make< TH2D >( "preDPhi_dPhiJet1MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_dPhiJet1MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "preDPhi_dPhiJet1MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_dPhiJet1MHTNoPhot" ]->Sumw2();

  histos2D_[ "preDPhi_dPhiJet2MHTNoPhot" ] = fs->make< TH2D >( "preDPhi_dPhiJet2MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_dPhiJet2MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "preDPhi_dPhiJet2MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_dPhiJet2MHTNoPhot" ]->Sumw2();

  histos2D_[ "preDPhi_dPhiJet3MHTNoPhot" ] = fs->make< TH2D >( "preDPhi_dPhiJet3MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_dPhiJet3MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "preDPhi_dPhiJet3MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_dPhiJet3MHTNoPhot" ]->Sumw2();

  histos2D_[ "preDPhi_genMHTNoPhotBins" ] = fs->make< TH2D >( "preDPhi_genMHTNoPhotBins", "Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_genMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos2D_[ "preDPhi_genMHTNoPhotBins" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_genMHTNoPhotBins" ]->Sumw2();

  histos2D_[ "preDPhi_genMHTNoPhot" ] = fs->make< TH2D >( "preDPhi_genMHTNoPhot", "Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_genMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos2D_[ "preDPhi_genMHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_genMHTNoPhot" ]->Sumw2();

  histos2D_[ "preDPhi_genMeffNoPhot" ] = fs->make< TH2D >( "preDPhi_genMeffNoPhot", "Event M^{GEN,no #gamma}_{eff}", 150, 0.0, 7500.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "preDPhi_genMeffNoPhot" ]->SetXTitle( "Event M^{GEN,no #gamma}_{eff}" );
  histos2D_[ "preDPhi_genMeffNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "preDPhi_genMeffNoPhot" ]->Sumw2();

  histos3D_[ "preDPhi_boson1PtvsMHTNoPhot" ] = fs->make< TH3D >( "preDPhi_boson1PtvsMHTNoPhot", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preDPhi_boson1PtvsMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "preDPhi_boson1PtvsMHTNoPhot" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson1PtvsMHTNoPhot" ]->Sumw2();

  histos3D_[ "preDPhi_boson2PtvsMHTNoPhot" ] = fs->make< TH3D >( "preDPhi_boson2PtvsMHTNoPhot", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "preDPhi_boson2PtvsMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "preDPhi_boson2PtvsMHTNoPhot" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson2PtvsMHTNoPhot" ]->Sumw2();

  histos3D_[ "preDPhi_MHTNoPhotvsHT" ] = fs->make< TH3D >( "preDPhi_MHTNoPhotvsHT", "Event #slashH^{no #gamma}_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 200, 0.0, 2000.0, 15,-0.5,14.5);
  histos3D_[ "preDPhi_MHTNoPhotvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preDPhi_MHTNoPhotvsHT" ]->SetYTitle( "Event #slashH^{no #gamma}_{T}" );
  histos3D_[ "preDPhi_MHTNoPhotvsHT" ]->Sumw2();

  histos3D_[ "preDPhi_boson1PtvsMHTNoPhotBins" ] = fs->make< TH3D >( "preDPhi_boson1PtvsMHTNoPhotBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preDPhi_boson1PtvsMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "preDPhi_boson1PtvsMHTNoPhotBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson1PtvsMHTNoPhotBins" ]->Sumw2();

  histos3D_[ "preDPhi_boson2PtvsMHTNoPhotBins" ] = fs->make< TH3D >( "preDPhi_boson2PtvsMHTNoPhotBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preDPhi_boson2PtvsMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "preDPhi_boson2PtvsMHTNoPhotBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "preDPhi_boson2PtvsMHTNoPhotBins" ]->Sumw2();

  histos3D_[ "preDPhi_MHTNoPhotvsHTBins" ] = fs->make< TH3D >( "preDPhi_MHTNoPhotvsHTBins", "Event #slashH^{no #gamma}_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "preDPhi_MHTNoPhotvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "preDPhi_MHTNoPhotvsHTBins" ]->SetYTitle( "Event #slashH^{no #gamma}_{T}" );
  histos3D_[ "preDPhi_MHTNoPhotvsHTBins" ]->Sumw2();
  */

  /////post RA2 Cuts
  histos2D_[ "nGenBoson" ] = fs->make< TH2D >( "numGenBoson", "#boson^{GEN}", 5, -0.5,4.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "nGenBoson" ]->SetXTitle( "#boson^{GEN}" );
  histos2D_[ "nGenBoson" ]->SetYTitle( "# of events" );
  histos2D_[ "nGenBoson" ]->Sumw2();

  histos2D_[ "nGenJets_Pt30" ] = fs->make< TH2D >( "nGenJets_Pt30", "#Jets^{GEN} (p_{T} > 30 GeV, |#eta| < 5.0)", 25, -0.5,24.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "nGenJets_Pt30" ]->SetXTitle( "#Jets^{GEN}" );
  histos2D_[ "nGenJets_Pt30" ]->SetYTitle( "# of events" );
  histos2D_[ "nGenJets_Pt30" ]->Sumw2();

  histos2D_[ "nGenJets_Pt50Eta25" ] = fs->make< TH2D >( "nGenJets_Pt50Eta25", "#Jets^{GEN} (p_{T} > 50 GeV, |#eta| < 2.5)", 25, -0.5,24.5, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "nGenJets_Pt50Eta25" ]->SetXTitle( "#Jets^{GEN}" );
  histos2D_[ "nGenJets_Pt50Eta25" ]->SetYTitle( "# of events" );
  histos2D_[ "nGenJets_Pt50Eta25" ]->Sumw2();

  histos2D_[ "dPhiJet1MHT" ] = fs->make< TH2D >( "dPhiJet1MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "dPhiJet1MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "dPhiJet1MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "dPhiJet1MHT" ]->Sumw2();

  histos2D_[ "dPhiJet2MHT" ] = fs->make< TH2D >( "dPhiJet2MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "dPhiJet2MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "dPhiJet2MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "dPhiJet2MHT" ]->Sumw2();

  histos2D_[ "dPhiJet3MHT" ] = fs->make< TH2D >( "dPhiJet3MHT", "#Delta#phi(J,#slashH_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "dPhiJet3MHT" ]->SetXTitle( "#Delta#phi(J, #slashH_{T})" );
  histos2D_[ "dPhiJet3MHT" ]->SetYTitle( "# of events" );
  histos2D_[ "dPhiJet3MHT" ]->Sumw2();

  histos2D_[ "genHTBins" ] = fs->make< TH2D >( "genHTBins", "Event H^{GEN}_{T}", htBins_.size()-1, htBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "genHTBins" ]->SetXTitle( "Event H^{GEN}_{T}" );
  histos2D_[ "genHTBins" ]->SetYTitle( "# of events" );
  histos2D_[ "genHTBins" ]->Sumw2();

  histos2D_[ "genMHTBins" ] = fs->make< TH2D >( "genMHTBins", "Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "genMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos2D_[ "genMHTBins" ]->SetYTitle( "# of events" );
  histos2D_[ "genMHTBins" ]->Sumw2();

  histos2D_[ "genHT" ] = fs->make< TH2D >( "genHT", "Event H^{GEN}_{T}", 100, 0, 5000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "genHT" ]->SetXTitle( "Event H^{GEN}_{T}" );
  histos2D_[ "genHT" ]->SetYTitle( "# of events" );
  histos2D_[ "genHT" ]->Sumw2();

  histos2D_[ "genMHT" ] = fs->make< TH2D >( "genMHT", "Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "genMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos2D_[ "genMHT" ]->SetYTitle( "# of events" );
  histos2D_[ "genMHT" ]->Sumw2();

  histos2D_[ "genMeff" ] = fs->make< TH2D >( "genMeff", "Event M^{GEN}_{eff}", 150, 0.0, 7500.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "genMeff" ]->SetXTitle( "Event M^{GEN}_{eff}" );
  histos2D_[ "genMeff" ]->SetYTitle( "# of events" );
  histos2D_[ "genMeff" ]->Sumw2();

  histos2D_[ "boson1Pt" ] = fs->make< TH2D >( "boson1Pt", "boson^{GEN} p_{T}", 60, 0.0, 600.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "boson1Pt" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "boson1Pt" ]->SetYTitle( "# of events" );
  histos2D_[ "boson1Pt" ]->Sumw2();

  histos2D_[ "boson1Eta" ] = fs->make< TH2D >( "boson1Eta", "boson^{GEN} #eta", 50, -5.0, 5.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "boson1Eta" ]->SetXTitle( "boson #eta" );
  histos2D_[ "boson1Eta" ]->SetYTitle( "# of events" );
  histos2D_[ "boson1Eta" ]->Sumw2();

  histos2D_[ "boson1PtBins" ] = fs->make< TH2D >( "boson1PtBins", "boson^{GEN} p_{T}", bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "boson1PtBins" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "boson1PtBins" ]->SetYTitle( "# of events" );
  histos2D_[ "boson1PtBins" ]->Sumw2();

  histos2D_[ "boson1EtaBins" ] = fs->make< TH2D >( "boson1EtaBins", "boson^{GEN} #eta", bosonEtaBins_.size()-1, bosonEtaBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "boson1EtaBins" ]->SetXTitle( "boson #eta" );
  histos2D_[ "boson1EtaBins" ]->SetYTitle( "# of events" );
  histos2D_[ "boson1EtaBins" ]->Sumw2();

  histos2D_[ "boson1Phi" ] = fs->make< TH2D >( "boson1Phi", "boson^{GEN} #phi", 50, -M_PI, M_PI , nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "boson1Phi" ]->SetXTitle( "boson #phi" );
  histos2D_[ "boson1Phi" ]->SetYTitle( "# of events" );
  histos2D_[ "boson1Phi" ]->Sumw2();

  /*
  histos2D_[ "boson2Pt" ] = fs->make< TH2D >( "boson2Pt", "boson^{GEN} p_{T}", 60, 0.0, 600.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "boson2Pt" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "boson2Pt" ]->SetYTitle( "# of events" );
  histos2D_[ "boson2Pt" ]->Sumw2();

  histos2D_[ "boson2Eta" ] = fs->make< TH2D >( "boson2Eta", "boson^{GEN} #eta", 50, -5.0, 5.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "boson2Eta" ]->SetXTitle( "boson #eta" );
  histos2D_[ "boson2Eta" ]->SetYTitle( "# of events" );
  histos2D_[ "boson2Eta" ]->Sumw2();

  histos2D_[ "boson2PtBins" ] = fs->make< TH2D >( "boson2PtBins", "boson^{GEN} p_{T}", bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "boson2PtBins" ]->SetXTitle( "boson p_{T}" );
  histos2D_[ "boson2PtBins" ]->SetYTitle( "# of events" );
  histos2D_[ "boson2PtBins" ]->Sumw2();

  histos2D_[ "boson2EtaBins" ] = fs->make< TH2D >( "boson2EtaBins", "boson^{GEN} #eta", bosonEtaBins_.size()-1, bosonEtaBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "boson2EtaBins" ]->SetXTitle( "boson #eta" );
  histos2D_[ "boson2EtaBins" ]->SetYTitle( "# of events" );
  histos2D_[ "boson2EtaBins" ]->Sumw2();

  histos2D_[ "boson2Phi" ] = fs->make< TH2D >( "boson2Phi", "boson^{GEN} #phi", 50, -M_PI, M_PI , nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "boson2Phi" ]->SetXTitle( "boson #phi" );
  histos2D_[ "boson2Phi" ]->SetYTitle( "# of events" );
  histos2D_[ "boson2Phi" ]->Sumw2();
  */
  //eta vs. pt
  histos3D_[ "boson1PtvsHT" ] = fs->make< TH3D >( "boson1PtvsHT", "boson^{GEN} p_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "boson1PtvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "boson1PtvsHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson1PtvsHT" ]->Sumw2();

  histos3D_[ "boson1PtvsMHT" ] = fs->make< TH3D >( "boson1PtvsMHT", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "boson1PtvsMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "boson1PtvsMHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson1PtvsMHT" ]->Sumw2();
  /*
  histos3D_[ "boson2PtvsHT" ] = fs->make< TH3D >( "boson2PtvsHT", "boson^{GEN} p_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "boson2PtvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "boson2PtvsHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson2PtvsHT" ]->Sumw2();

  histos3D_[ "boson2PtvsMHT" ] = fs->make< TH3D >( "boson2PtvsMHT", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "boson2PtvsMHT" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "boson2PtvsMHT" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson2PtvsMHT" ]->Sumw2();
  */
  histos3D_[ "MHTvsHT" ] = fs->make< TH3D >( "MHTvsHT", "Event #slashH_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 200, 0.0, 2000.0, 15,-0.5,14.5);
  histos3D_[ "MHTvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "MHTvsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos3D_[ "MHTvsHT" ]->Sumw2();

  histos3D_[ "boson1PtvsHTBins" ] = fs->make< TH3D >( "boson1PtvsHTBins", "boson^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "boson1PtvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "boson1PtvsHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson1PtvsHTBins" ]->Sumw2();

  histos3D_[ "boson1PtvsMHTBins" ] = fs->make< TH3D >( "boson1PtvsMHTBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "boson1PtvsMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "boson1PtvsMHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson1PtvsMHTBins" ]->Sumw2();

  /*
  histos3D_[ "boson2PtvsHTBins" ] = fs->make< TH3D >( "boson2PtvsHTBins", "boson^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "boson2PtvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "boson2PtvsHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson2PtvsHTBins" ]->Sumw2();

  histos3D_[ "boson2PtvsMHTBins" ] = fs->make< TH3D >( "boson2PtvsMHTBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "boson2PtvsMHTBins" ]->SetXTitle( "Event #slashH^{GEN}_{T}" );
  histos3D_[ "boson2PtvsMHTBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson2PtvsMHTBins" ]->Sumw2();
  */
  histos3D_[ "MHTvsHTBins" ] = fs->make< TH3D >( "MHTvsHTBins", "Event #slashH_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "MHTvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "MHTvsHTBins" ]->SetYTitle( "Event #slashH_{T}" );
  histos3D_[ "MHTvsHTBins" ]->Sumw2();

  /*
  //plots with photon removed mht 
  histos2D_[ "dPhiJet1MHTNoPhot" ] = fs->make< TH2D >( "dPhiJet1MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "dPhiJet1MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "dPhiJet1MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "dPhiJet1MHTNoPhot" ]->Sumw2();

  histos2D_[ "dPhiJet2MHTNoPhot" ] = fs->make< TH2D >( "dPhiJet2MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "dPhiJet2MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "dPhiJet2MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "dPhiJet2MHTNoPhot" ]->Sumw2();

  histos2D_[ "dPhiJet3MHTNoPhot" ] = fs->make< TH2D >( "dPhiJet3MHTNoPhot", "#Delta#phi(J,#slashH^{no #gamma}_{T})", 60, -1.0, 3.2, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "dPhiJet3MHTNoPhot" ]->SetXTitle( "#Delta#phi(J, #slashH^{no #gamma}_{T})" );
  histos2D_[ "dPhiJet3MHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "dPhiJet3MHTNoPhot" ]->Sumw2();

  histos2D_[ "genMHTNoPhotBins" ] = fs->make< TH2D >( "genMHTNoPhotBins", "Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "genMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos2D_[ "genMHTNoPhotBins" ]->SetYTitle( "# of events" );
  histos2D_[ "genMHTNoPhotBins" ]->Sumw2();

  histos2D_[ "genMHTNoPhot" ] = fs->make< TH2D >( "genMHTNoPhot", "Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "genMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos2D_[ "genMHTNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "genMHTNoPhot" ]->Sumw2();

  histos2D_[ "genMeffNoPhot" ] = fs->make< TH2D >( "genMeffNoPhot", "Event M^{GEN,no #gamma}_{eff}", 150, 0.0, 7500.0, nJetBins_.size()-1, nJetBins_.data());
  histos2D_[ "genMeffNoPhot" ]->SetXTitle( "Event M^{GEN,no #gamma}_{eff}" );
  histos2D_[ "genMeffNoPhot" ]->SetYTitle( "# of events" );
  histos2D_[ "genMeffNoPhot" ]->Sumw2();

  histos3D_[ "boson1PtvsMHTNoPhot" ] = fs->make< TH3D >( "boson1PtvsMHTNoPhot", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "boson1PtvsMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "boson1PtvsMHTNoPhot" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson1PtvsMHTNoPhot" ]->Sumw2();

  histos3D_[ "boson2PtvsMHTNoPhot" ] = fs->make< TH3D >( "boson2PtvsMHTNoPhot", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", 200, 0.0, 2000.0, 60, 0.0, 600.0, 15,-0.5,14.5);
  histos3D_[ "boson2PtvsMHTNoPhot" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "boson2PtvsMHTNoPhot" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson2PtvsMHTNoPhot" ]->Sumw2();

  histos3D_[ "MHTNoPhotvsHT" ] = fs->make< TH3D >( "MHTNoPhotvsHT", "Event #slashH^{no #gamma}_{T} vs. Event H_{T}", 100, 0.0, 5000.0, 200, 0.0, 2000.0, 15,-0.5,14.5);
  histos3D_[ "MHTNoPhotvsHT" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "MHTNoPhotvsHT" ]->SetYTitle( "Event #slashH^{no #gamma}_{T}" );
  histos3D_[ "MHTNoPhotvsHT" ]->Sumw2();

  histos3D_[ "boson1PtvsMHTNoPhotBins" ] = fs->make< TH3D >( "boson1PtvsMHTNoPhotBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "boson1PtvsMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "boson1PtvsMHTNoPhotBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson1PtvsMHTNoPhotBins" ]->Sumw2();

  histos3D_[ "boson2PtvsMHTNoPhotBins" ] = fs->make< TH3D >( "boson2PtvsMHTNoPhotBins", "boson^{GEN} p_{T} vs. Event #slashH^{GEN,no #gamma}_{T}", mhtBins_.size()-1, mhtBins_.data(), bosonPtBins_.size()-1, bosonPtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "boson2PtvsMHTNoPhotBins" ]->SetXTitle( "Event #slashH^{GEN,no #gamma}_{T}" );
  histos3D_[ "boson2PtvsMHTNoPhotBins" ]->SetYTitle( "boson p_{T}" );
  histos3D_[ "boson2PtvsMHTNoPhotBins" ]->Sumw2();

  histos3D_[ "MHTNoPhotvsHTBins" ] = fs->make< TH3D >( "MHTNoPhotvsHTBins", "Event #slashH^{no #gamma}_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data(), nJetBins_.size()-1, nJetBins_.data());
  histos3D_[ "MHTNoPhotvsHTBins" ]->SetXTitle( "Event H_{T}" );
  histos3D_[ "MHTNoPhotvsHTBins" ]->SetYTitle( "Event #slashH^{no #gamma}_{T}" );
  histos3D_[ "MHTNoPhotvsHTBins" ]->Sumw2();
  */
}

// ------------ method called once each job just after ending the event loop  ------------
void GenStudy::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenStudy);
