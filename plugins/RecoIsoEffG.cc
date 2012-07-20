// -*- C++ -*-
//
// Package:    Photons
// Class:      Photons
// 
/**\class Photons RecoIsoEffG.cc ZInvisibleBkgds/Photons/plugins/RecoIsoEffG.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: RecoIsoEffG.cc,v 1.1 2012/07/09 14:29:00 sturdy Exp $
//
//


#include "ZInvisibleBkgds/Photons/interface/RecoIsoEffG.h"

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
RecoIsoEffG::RecoIsoEffG(const edm::ParameterSet& pset) :
  debug_           ( pset.getParameter< bool >( "debug" ) ),
  debugString_     ( pset.getParameter< std::string >( "debugString" ) ),
  genLabel_        ( pset.getParameter< edm::InputTag >( "genLabel" ) ),
  genStatus_       ( pset.getParameter< int >( "genStatus" ) ),

  doPUReweight_ ( pset.getParameter< bool >( "doPUReweight" ) ),
  puWeightLabel_( pset.getParameter< edm::InputTag >( "puWeight" ) ),

  jetHTLabel_ ( pset.getParameter< edm::InputTag >( "jetHTLabel" ) ),
  jetMHTLabel_( pset.getParameter< edm::InputTag >( "jetMHTLabel" ) ),
  htLabel_    ( pset.getParameter< edm::InputTag >( "htLabel" ) ),
  mhtLabel_   ( pset.getParameter< edm::InputTag >( "mhtLabel" ) ),

  photonLabel_     ( pset.getParameter< edm::InputTag >( "photonLabel" ) ),
  isoPhotonLabel_  ( pset.getParameter< edm::InputTag >( "isoPhotonLabel" ) ),

  photonMinPt_   ( pset.getParameter< double >( "photonMinPt" ) ),
  photonEBMaxEta_( pset.getParameter< double >( "photonEBMaxEta" ) ),
  photonEEMinEta_( pset.getParameter< double >( "photonEEMinEta" ) ),
  photonEEMaxEta_( pset.getParameter< double >( "photonEEMaxEta" ) )
{
  photonPtBins_  = pset.getParameter< std::vector<double> >( "photonPtBins" );
  photonEtaBins_ = pset.getParameter< std::vector<double> >( "photonEtaBins" );
  htBins_        = pset.getParameter< std::vector<double> >( "htBins" );
  mhtBins_       = pset.getParameter< std::vector<double> >( "mhtBins" );
  //register your products
  /* Examples
     produces<ExampleData2>();

     //if do put with a label
     produces<ExampleData2>("label");
 
     //if you want to put into the Run
     produces<ExampleData2,InRun>();
  */
  //now do what ever other initialization is needed
  
}


RecoIsoEffG::~RecoIsoEffG()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void RecoIsoEffG::produce(edm::Event& ev, const edm::EventSetup& es)
{
  using namespace edm;
  //read in the gen particles
  edm::Handle<reco::GenParticleCollection> gens;
  ev.getByLabel(genLabel_,gens);

  //get the ht
  double eventWeightPU = 1.;
  edm::Handle<double> puWeight;
  if (doPUReweight_) {
    ev.getByLabel(puWeightLabel_,puWeight);
    eventWeightPU = *puWeight;
  }
  //get the jets
  edm::Handle<edm::View<pat::Jet> > jetsHT;
  ev.getByLabel(jetHTLabel_,jetsHT);

  edm::Handle<edm::View<pat::Jet> > jetsMHT;
  ev.getByLabel(jetMHTLabel_,jetsMHT);
  //get the ht
  edm::Handle<double> ht;
  ev.getByLabel(htLabel_,ht);
  //get the mht
  edm::Handle<edm::View<reco::MET> > mht;
  ev.getByLabel(mhtLabel_,mht);

  //get the photons
  edm::Handle<edm::View<pat::Photon> > photons;
  ev.getByLabel(photonLabel_,photons);

  edm::Handle<edm::View<pat::Photon> > photonsIso;
  ev.getByLabel(isoPhotonLabel_,photonsIso);

  /*
    reconstruction/isolation efficiency in each bin found by matching the 
    gen photon to a reco photon passing isolation cuts
  */

  double g_pt = -9.;
  double r_pt = -9.;
  double i_pt = -9.;

  histos1D_["photonRecoIsoEff_NJets_Pt30"     ]->Fill( jetsMHT->size(), eventWeightPU);
  histos1D_["photonRecoIsoEff_NJets_Pt50Eta25"]->Fill( jetsHT->size(),  eventWeightPU);
  histos1D_["photonRecoIsoEff_HT" ]->Fill( *ht,            eventWeightPU);
  histos1D_["photonRecoIsoEff_MHT"]->Fill( (*mht)[0].pt(), eventWeightPU);

  histos1D_["recoPhotonRecoIsoEff_N" ]->Fill( photons->size() , eventWeightPU);
  if (photons->size() > 0) {
    r_pt = (*photons)[0].pt();
    const reco::Candidate* candPhot = (*photons)[0].genPhoton();

    double rg_pt = -9.;
    if (candPhot)
      rg_pt = candPhot->pt();
    //double rg_pt = (*photons)[0].genParticle()->pt();
    histos1D_["recoPhotonRecoIsoEff_Pt" ]->Fill( (*photons)[0].pt() , eventWeightPU);
    histos1D_["recoPhotonRecoIsoEff_Eta"]->Fill( (*photons)[0].eta(), eventWeightPU);
    histos1D_["recoPhotonRecoIsoEff_Phi"]->Fill( (*photons)[0].phi(), eventWeightPU);
    histos2D_["photonRecoIsoEff_EmbGENVsRECOPt"]->Fill( r_pt, rg_pt, eventWeightPU);
  }
  histos1D_["isoPhotonRecoIsoEff_N" ]->Fill( photonsIso->size() , eventWeightPU);
  if (photonsIso->size() > 0) {
    i_pt = (*photonsIso)[0].pt();
    const reco::Candidate* candPhot = (*photonsIso)[0].genPhoton();
    double ig_pt = 9.0;
    if (candPhot)
      ig_pt = candPhot->pt();
    //double ig_pt = (*photonsIso)[0].genParticle()->pt();
    histos1D_["isoPhotonRecoIsoEff_Pt" ]->Fill( (*photonsIso)[0].pt() , eventWeightPU);
    histos1D_["isoPhotonRecoIsoEff_Eta"]->Fill( (*photonsIso)[0].eta(), eventWeightPU);
    histos1D_["isoPhotonRecoIsoEff_Phi"]->Fill( (*photonsIso)[0].phi(), eventWeightPU);
    histos2D_["photonRecoIsoEff_EmbGENVsISOPt"] ->Fill( i_pt, ig_pt, eventWeightPU);
  }
  if (gens->size() > 0)
    g_pt = (*gens)[0].pt();

  histos2D_["photonRecoIsoEff_GENVsRECOPt"]->Fill( r_pt, g_pt, eventWeightPU);
  histos2D_["photonRecoIsoEff_GENVsISOPt"] ->Fill( i_pt, g_pt, eventWeightPU);
  
  //loop over gen particle collection
  int numMatchedGens = 0;
  reco::GenParticleCollection::const_iterator gen = gens->begin();
  if (debug_)
    std::cout<<debugString_<<"::looping over gen particles size = "<<gens->size()<<std::endl;
  for ( ; gen != gens->end(); ++gen) {
    
    if (debug_) {
      std::cout<<"pdgID  --  status  --  #mom  --  {mom(ID,St,#mom,({pdgID,st}),pT,eta),...}"<<std::endl;
      std::cout<<gen->pdgId()<<"  --  "<<gen->status()<<"  --  "<<gen->numberOfMothers()<<"  --  "<<"{";
      for (unsigned int dau = 0; dau < gen->numberOfMothers(); ++dau) {
	std::cout<<"("<<gen->mother(dau)->pdgId()<<","<<gen->mother(dau)->status()<<","<<gen->mother(dau)->numberOfMothers()<<"(";
	for (unsigned int mom = 0; mom < gen->mother(dau)->numberOfMothers(); ++mom) 
	  std::cout<<"{"<<gen->mother(dau)->mother(mom)->pdgId()<<","<<gen->mother(dau)->mother(mom)->status()<<"},";
	std::cout<<"),";
	std::cout<<gen->pt()<<","<<gen->eta()<<") ";
      }
      std::cout<<"}"<<std::endl;
    }
    
    double gphot_eta = gen->eta();
    double gphot_phi = gen->phi();
    
    bool   foundMatch = false;
    double dRMinPhot=99.0;
    int    dRMinPhotIdx = -1;
    int iphot = 0;

    edm::View<pat::Photon>::const_iterator recop = photons->begin();
    if (debug_)
      std::cout<<debugString_<<"::looping over reco photons size = "<<photons->size()<<std::endl;
    for ( ; recop != photons->end(); ++recop) {
      const std::vector<std::string> floatNames = recop->userFloatNames();
      //std::cout<<floatNames.at(0)<<std::endl;
      if (debug_) {
	recop->listUserFloat();
	std::cout<<"RecoIsoEffG::rho25::"  <<recop->userFloat("rho25")<<std::endl;
	std::cout<<"RecoIsoEffG::rhoToPhotonMap::"  <<recop->userFloat("rhoToPhotonMap")<<std::endl;
	std::cout<<"RecoIsoEffG::chargedIsolation::"<<recop->userIsolation(pat::IsolationKeys(pat::UserBaseIso+0))<<std::endl;
	std::cout<<"RecoIsoEffG::neutralIsolation::"<<recop->userIsolation(pat::IsolationKeys(pat::UserBaseIso+2))<<std::endl;
	std::cout<<"RecoIsoEffG::gammaIsolation::"  <<recop->userIsolation(pat::IsolationKeys(pat::UserBaseIso+3))<<std::endl;
	std::cout<<"RecoIsoEffG::chargedIsolation::"<<recop->userIsolation(pat::User1Iso)<<std::endl;
	std::cout<<"RecoIsoEffG::neutralIsolation::"<<recop->userIsolation(pat::User3Iso)<<std::endl;
	std::cout<<"RecoIsoEffG::gammaIsolation::"  <<recop->userIsolation(pat::User4Iso)<<std::endl;
      }
      double rphot_eta = recop->eta();
      double rphot_phi = recop->phi();
      double dR = reco::deltaR(rphot_eta, rphot_phi, gphot_eta, gphot_phi);
      if( dR<dRMinPhot) {
	dRMinPhot    = dR;
	dRMinPhotIdx = iphot;
      }
      ++iphot;
    }
    if( dRMinPhotIdx>=0 && dRMinPhot<0.1 ) {
      foundMatch = true;
      ++numMatchedGens;
    } else {
      //std::cout << "not Matching to leading Phot " << " dRMin " << dRMin << " NRecoPhotons " << IsoPhotons.size() << std::endl;	    
    }

    ////////////// no need to match to jet and clean as this will have been done by the photon-jet cleaning
    //step.  Also, HT/MHT will be only for the cleaned jets, what about only a single isolated photon?
    // in that case, some plots will be doubly filled.  Do I filter out the events with 2 isolated photons
    // for the reco/iso study?
    // maybe we fill these plots for each photon, as a sort of increased statistics measure?
    // needs to be understood, but overall, not my priority at the moment
    histos1D_["photonRecoIsoEff_Pt" ]->Fill( gen->pt() , eventWeightPU);
    histos1D_["photonRecoIsoEff_Eta"]->Fill( gen->eta(), eventWeightPU);
    histos1D_["photonRecoIsoEff_Phi"]->Fill( gen->phi(), eventWeightPU);

    histos2D_["photonRecoIsoEff_PtVsHT" ]->Fill( *ht, gen->pt(),      eventWeightPU);
    histos2D_["photonRecoIsoEff_MHTVsHT"]->Fill( *ht, (*mht)[0].pt(), eventWeightPU);

    if (foundMatch) {
      histos1D_["matchedPhotonRecoIsoEff_Pt" ]->Fill( gen->pt() , eventWeightPU);
      histos1D_["matchedPhotonRecoIsoEff_Eta"]->Fill( gen->eta(), eventWeightPU);
      histos1D_["matchedPhotonRecoIsoEff_Phi"]->Fill( gen->phi(), eventWeightPU);

      histos1D_["matchedPhotonRecoIsoEff_NJets_Pt30"     ]->Fill( jetsMHT->size(), eventWeightPU);
      histos1D_["matchedPhotonRecoIsoEff_NJets_Pt50Eta25"]->Fill( jetsHT->size(),  eventWeightPU);
      histos1D_["matchedPhotonRecoIsoEff_HT" ]->Fill( *ht,             eventWeightPU);
      histos1D_["matchedPhotonRecoIsoEff_MHT"]->Fill( (*mht)[0].pt() , eventWeightPU);
      
      histos2D_["matchedPhotonRecoIsoEff_PtVsHT" ]->Fill( *ht, gen->pt(),      eventWeightPU);
      histos2D_["matchedPhotonRecoIsoEff_MHTVsHT"]->Fill( *ht, (*mht)[0].pt(), eventWeightPU);
    }
  }//end loop over gen particles
  histos1D_["numRecoBoson"       ]->Fill( photons->size(), eventWeightPU);
  histos1D_["numGenMatchedBoson" ]->Fill( numMatchedGens,  eventWeightPU);

}

// ------------ method called once each job just before starting event loop  ------------
void RecoIsoEffG::beginJob()
{
  edm::Service< TFileService > fs;
  
  //char title[128];
  //sprintf(title,"# of events passing HLT_%s_v*",unbiasedTrigger_.c_str());

  histos1D_[ "numGenMatchedBoson" ] = fs->make< TH1D >( "numGenMatchedBoson", "# #gamma^{GEN} matched to #gamma^{RECOISO}", 5, -0.5,4.5);
  histos1D_[ "numGenMatchedBoson" ]->SetXTitle( "##gamma^{GEN}" );
  histos1D_[ "numGenMatchedBoson" ]->SetYTitle( "# of events" );
  histos1D_[ "numGenMatchedBoson" ]->Sumw2();

  histos1D_[ "numRecoBoson" ] = fs->make< TH1D >( "numRecoBoson", "# #gamma^{RECOISO} ", 5, -0.5,4.5);
  histos1D_[ "numRecoBoson" ]->SetXTitle( "# #gamma^{RECOISO}" );
  histos1D_[ "numRecoBoson" ]->SetYTitle( "# of events" );
  histos1D_[ "numRecoBoson" ]->Sumw2();

  histos1D_[ "photonRecoIsoEff_NJets_Pt30" ] = fs->make< TH1D >( "photonRecoIsoEff_NJets_Pt30", "#Jets (p_{T} > 30 GeV, |#eta| < 5.0)", 25, -0.5,24.5);
  histos1D_[ "photonRecoIsoEff_NJets_Pt30" ]->SetXTitle( "#Jets" );
  histos1D_[ "photonRecoIsoEff_NJets_Pt30" ]->SetYTitle( "# of events" );
  histos1D_[ "photonRecoIsoEff_NJets_Pt30" ]->Sumw2();

  histos1D_[ "photonRecoIsoEff_NJets_Pt50Eta25" ] = fs->make< TH1D >( "photonRecoIsoEff_NJets_Pt50Eta25", "#Jets (p_{T} > 50 GeV, |#eta| < 2.5)", 25, -0.5,24.5);
  histos1D_[ "photonRecoIsoEff_NJets_Pt50Eta25" ]->SetXTitle( "#Jets" );
  histos1D_[ "photonRecoIsoEff_NJets_Pt50Eta25" ]->SetYTitle( "# of events" );
  histos1D_[ "photonRecoIsoEff_NJets_Pt50Eta25" ]->Sumw2();

  histos1D_[ "photonRecoIsoEff_HT" ] = fs->make< TH1D >( "photonRecoIsoEff_HT", "Event H_{T}", htBins_.size()-1, htBins_.data());
  histos1D_[ "photonRecoIsoEff_HT" ]->SetXTitle( "Event H_{T}" );
  histos1D_[ "photonRecoIsoEff_HT" ]->SetYTitle( "# of events" );
  histos1D_[ "photonRecoIsoEff_HT" ]->Sumw2();

  histos1D_[ "photonRecoIsoEff_MHT" ] = fs->make< TH1D >( "photonRecoIsoEff_MHT", "Event #slashH_{T}", mhtBins_.size()-1, mhtBins_.data());
  histos1D_[ "photonRecoIsoEff_MHT" ]->SetXTitle( "Event #slashH_{T}" );
  histos1D_[ "photonRecoIsoEff_MHT" ]->SetYTitle( "# of events" );
  histos1D_[ "photonRecoIsoEff_MHT" ]->Sumw2();

  histos1D_[ "photonRecoIsoEff_Pt" ] = fs->make< TH1D >( "photonRecoIsoEff_Pt", "#gamma^{GEN} p_{T}", photonPtBins_.size()-1, photonPtBins_.data());
  histos1D_[ "photonRecoIsoEff_Pt" ]->SetXTitle( "#gamma p_{T}" );
  histos1D_[ "photonRecoIsoEff_Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "photonRecoIsoEff_Pt" ]->Sumw2();

  histos1D_[ "photonRecoIsoEff_Eta" ] = fs->make< TH1D >( "photonRecoIsoEff_Eta", "#gamma^{GEN} #eta", photonEtaBins_.size()-1, photonEtaBins_.data());
  histos1D_[ "photonRecoIsoEff_Eta" ]->SetXTitle( "#gamma #eta" );
  histos1D_[ "photonRecoIsoEff_Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "photonRecoIsoEff_Eta" ]->Sumw2();

  histos1D_[ "photonRecoIsoEff_Phi" ] = fs->make< TH1D >( "photonRecoIsoEff_Phi", "#gamma^{GEN} #phi", 50, -M_PI, M_PI );
  histos1D_[ "photonRecoIsoEff_Phi" ]->SetXTitle( "#gamma #phi" );
  histos1D_[ "photonRecoIsoEff_Phi" ]->SetYTitle( "# of events" );
  histos1D_[ "photonRecoIsoEff_Phi" ]->Sumw2();

  histos1D_[ "recoPhotonRecoIsoEff_N" ] = fs->make< TH1D >( "recoPhotonRecoIsoEff_N", "##gamma^{RECO}", 6, -0.5, 5.5);
  histos1D_[ "recoPhotonRecoIsoEff_N" ]->SetXTitle( "# #gamma^{RECO}" );
  histos1D_[ "recoPhotonRecoIsoEff_N" ]->SetYTitle( "# of events" );
  histos1D_[ "recoPhotonRecoIsoEff_N" ]->Sumw2();

  histos1D_[ "recoPhotonRecoIsoEff_Pt" ] = fs->make< TH1D >( "recoPhotonRecoIsoEff_Pt", "#gamma^{RECO} p_{T}", photonPtBins_.size()-1, photonPtBins_.data());
  histos1D_[ "recoPhotonRecoIsoEff_Pt" ]->SetXTitle( "#gamma p_{T}" );
  histos1D_[ "recoPhotonRecoIsoEff_Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "recoPhotonRecoIsoEff_Pt" ]->Sumw2();

  histos1D_[ "recoPhotonRecoIsoEff_Eta" ] = fs->make< TH1D >( "recoPhotonRecoIsoEff_Eta", "#gamma^{RECO} #eta", photonEtaBins_.size()-1, photonEtaBins_.data());
  histos1D_[ "recoPhotonRecoIsoEff_Eta" ]->SetXTitle( "#gamma #eta" );
  histos1D_[ "recoPhotonRecoIsoEff_Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "recoPhotonRecoIsoEff_Eta" ]->Sumw2();

  histos1D_[ "recoPhotonRecoIsoEff_Phi" ] = fs->make< TH1D >( "recoPhotonRecoIsoEff_Phi", "#gamma^{RECO} #phi", 50, -M_PI, M_PI );
  histos1D_[ "recoPhotonRecoIsoEff_Phi" ]->SetXTitle( "#gamma #phi" );
  histos1D_[ "recoPhotonRecoIsoEff_Phi" ]->SetYTitle( "# of events" );
  histos1D_[ "recoPhotonRecoIsoEff_Phi" ]->Sumw2();

  histos1D_[ "isoPhotonRecoIsoEff_N" ] = fs->make< TH1D >( "isoPhotonRecoIsoEff_N", "##gamma^{ISO}", 6, -0.5, 5.5);
  histos1D_[ "isoPhotonRecoIsoEff_N" ]->SetXTitle( "# #gamma^{ISO}" );
  histos1D_[ "isoPhotonRecoIsoEff_N" ]->SetYTitle( "# of events" );
  histos1D_[ "isoPhotonRecoIsoEff_N" ]->Sumw2();

  histos1D_[ "isoPhotonRecoIsoEff_Pt" ] = fs->make< TH1D >( "isoPhotonRecoIsoEff_Pt", "#gamma^{ISO} p_{T}", photonPtBins_.size()-1, photonPtBins_.data());
  histos1D_[ "isoPhotonRecoIsoEff_Pt" ]->SetXTitle( "#gamma p_{T}" );
  histos1D_[ "isoPhotonRecoIsoEff_Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "isoPhotonRecoIsoEff_Pt" ]->Sumw2();

  histos1D_[ "isoPhotonRecoIsoEff_Eta" ] = fs->make< TH1D >( "isoPhotonRecoIsoEff_Eta", "#gamma^{ISO} #eta", photonEtaBins_.size()-1, photonEtaBins_.data());
  histos1D_[ "isoPhotonRecoIsoEff_Eta" ]->SetXTitle( "#gamma #eta" );
  histos1D_[ "isoPhotonRecoIsoEff_Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "isoPhotonRecoIsoEff_Eta" ]->Sumw2();

  histos1D_[ "isoPhotonRecoIsoEff_Phi" ] = fs->make< TH1D >( "isoPhotonRecoIsoEff_Phi", "#gamma^{ISO} #phi", 50, -M_PI, M_PI );
  histos1D_[ "isoPhotonRecoIsoEff_Phi" ]->SetXTitle( "#gamma #phi" );
  histos1D_[ "isoPhotonRecoIsoEff_Phi" ]->SetYTitle( "# of events" );
  histos1D_[ "isoPhotonRecoIsoEff_Phi" ]->Sumw2();

  //for gen boson matched to reco photon
  histos1D_[ "matchedPhotonRecoIsoEff_NJets_Pt30" ] = fs->make< TH1D >( "matchedPhotonRecoIsoEff_NJets_Pt30", "#Jets (p_{T} > 30 GeV, |#eta| < 5.0)", 25, -0.5,24.5);
  histos1D_[ "matchedPhotonRecoIsoEff_NJets_Pt30" ]->SetXTitle( "#Jets" );
  histos1D_[ "matchedPhotonRecoIsoEff_NJets_Pt30" ]->SetYTitle( "# of events" );
  histos1D_[ "matchedPhotonRecoIsoEff_NJets_Pt30" ]->Sumw2();

  histos1D_[ "matchedPhotonRecoIsoEff_NJets_Pt50Eta25" ] = fs->make< TH1D >( "matchedPhotonRecoIsoEff_NJets_Pt50Eta25", "#Jets (p_{T} > 50 GeV, |#eta| < 2.5)", 25, -0.5,24.5);
  histos1D_[ "matchedPhotonRecoIsoEff_NJets_Pt50Eta25" ]->SetXTitle( "#Jets" );
  histos1D_[ "matchedPhotonRecoIsoEff_NJets_Pt50Eta25" ]->SetYTitle( "# of events" );
  histos1D_[ "matchedPhotonRecoIsoEff_NJets_Pt50Eta25" ]->Sumw2();

  histos1D_[ "matchedPhotonRecoIsoEff_HT" ] = fs->make< TH1D >( "matchedPhotonRecoIsoEff_HT", "Event H_{T}", htBins_.size()-1, htBins_.data());
  histos1D_[ "matchedPhotonRecoIsoEff_HT" ]->SetXTitle( "Event H_{T}" );
  histos1D_[ "matchedPhotonRecoIsoEff_HT" ]->SetYTitle( "# of events" );
  histos1D_[ "matchedPhotonRecoIsoEff_HT" ]->Sumw2();

  histos1D_[ "matchedPhotonRecoIsoEff_MHT" ] = fs->make< TH1D >( "matchedPhotonRecoIsoEff_MHT", "Event #slashH_{T}", mhtBins_.size()-1, mhtBins_.data());
  histos1D_[ "matchedPhotonRecoIsoEff_MHT" ]->SetXTitle( "Event #slashH_{T}" );
  histos1D_[ "matchedPhotonRecoIsoEff_MHT" ]->SetYTitle( "# of events" );
  histos1D_[ "matchedPhotonRecoIsoEff_MHT" ]->Sumw2();

  histos1D_[ "matchedPhotonRecoIsoEff_Pt" ] = fs->make< TH1D >( "matchedPhotonRecoIsoEff_Pt", "#gamma^{GEN} p_{T}", photonPtBins_.size()-1, photonPtBins_.data());
  histos1D_[ "matchedPhotonRecoIsoEff_Pt" ]->SetXTitle( "#gamma p_{T}" );
  histos1D_[ "matchedPhotonRecoIsoEff_Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "matchedPhotonRecoIsoEff_Pt" ]->Sumw2();

  histos1D_[ "matchedPhotonRecoIsoEff_Eta" ] = fs->make< TH1D >( "matchedPhotonRecoIsoEff_Eta", "#gamma^{GEN} #eta", photonEtaBins_.size()-1, photonEtaBins_.data());
  histos1D_[ "matchedPhotonRecoIsoEff_Eta" ]->SetXTitle( "#gamma #eta" );
  histos1D_[ "matchedPhotonRecoIsoEff_Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "matchedPhotonRecoIsoEff_Eta" ]->Sumw2();

  histos1D_[ "matchedPhotonRecoIsoEff_Phi" ] = fs->make< TH1D >( "matchedPhotonRecoIsoEff_Phi", "#gamma^{GEN} #phi", 50, -M_PI, M_PI );
  histos1D_[ "matchedPhotonRecoIsoEff_Phi" ]->SetXTitle( "#gamma #phi" );
  histos1D_[ "matchedPhotonRecoIsoEff_Phi" ]->SetYTitle( "# of events" );
  histos1D_[ "matchedPhotonRecoIsoEff_Phi" ]->Sumw2();

  //eta vs. pt
  histos2D_[ "photonRecoIsoEff_PtVsHT" ] = fs->make< TH2D >( "photonRecoIsoEff_PtVsHT", "#gamma^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), photonPtBins_.size()-1, photonPtBins_.data());
  histos2D_[ "photonRecoIsoEff_PtVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "photonRecoIsoEff_PtVsHT" ]->SetYTitle( "#gamma p_{T}" );
  histos2D_[ "photonRecoIsoEff_PtVsHT" ]->Sumw2();

  histos2D_[ "photonRecoIsoEff_MHTVsHT" ] = fs->make< TH2D >( "photonRecoIsoEff_MHTVsHT", "Event #slashH_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data());
  histos2D_[ "photonRecoIsoEff_MHTVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "photonRecoIsoEff_MHTVsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos2D_[ "photonRecoIsoEff_MHTVsHT" ]->Sumw2();

  //for gen boson matched to reco photon
  histos2D_[ "matchedPhotonRecoIsoEff_PtVsHT" ] = fs->make< TH2D >( "matchedPhotonRecoIsoEff_PtVsHT", "#gamma^{GEN} p_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), photonPtBins_.size()-1, photonPtBins_.data());
  histos2D_[ "matchedPhotonRecoIsoEff_PtVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "matchedPhotonRecoIsoEff_PtVsHT" ]->SetYTitle( "#gamma p_{T}" );
  histos2D_[ "matchedPhotonRecoIsoEff_PtVsHT" ]->Sumw2();

  histos2D_[ "matchedPhotonRecoIsoEff_MHTVsHT" ] = fs->make< TH2D >( "matchedPhotonRecoIsoEff_MHTVsHT", "Event #slashH_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data());
  histos2D_[ "matchedPhotonRecoIsoEff_MHTVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "matchedPhotonRecoIsoEff_MHTVsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos2D_[ "matchedPhotonRecoIsoEff_MHTVsHT" ]->Sumw2();

  //for gen boson matched to reco photon
  histos2D_[ "photonRecoIsoEff_GENVsRECOPt" ] = fs->make< TH2D >( "photonRecoIsoEff_GENVsRECOPt", "#gamma^{GEN} p_{T} vs. #gamma^{RECO} p_{T}", 101, -10., 1000., 101, -10., 1000.);
  histos2D_[ "photonRecoIsoEff_GENVsRECOPt" ]->SetXTitle( "#gamma^{RECO} p_{T}" );
  histos2D_[ "photonRecoIsoEff_GENVsRECOPt" ]->SetYTitle( "#gamma^{GEN} p_{T}" );
  histos2D_[ "photonRecoIsoEff_GENVsRECOPt" ]->Sumw2();

  histos2D_[ "photonRecoIsoEff_GENVsISOPt" ] = fs->make< TH2D >( "photonRecoIsoEff_GENVsISOPt", "#gamma^{GEN} p_{T} vs. #gamma^{ISO} p_{T}", 101, -10., 1000., 101, -10., 1000.);
  histos2D_[ "photonRecoIsoEff_GENVsISOPt" ]->SetXTitle( "#gamma^{ISO} p_{T}" );
  histos2D_[ "photonRecoIsoEff_GENVsISOPt" ]->SetYTitle( "#gamma^{GEN} p_{T}" );
  histos2D_[ "photonRecoIsoEff_GENVsISOPt" ]->Sumw2();

  histos2D_[ "photonRecoIsoEff_EmbGENVsRECOPt" ] = fs->make< TH2D >( "photonRecoIsoEff_EmbGENVsRECOPt", "#gamma^{EmbGEN} p_{T} vs. #gamma^{RECO} p_{T}", 101, -10., 1000., 101, -10., 1000.);
  histos2D_[ "photonRecoIsoEff_EmbGENVsRECOPt" ]->SetXTitle( "#gamma^{RECO} p_{T}" );
  histos2D_[ "photonRecoIsoEff_EmbGENVsRECOPt" ]->SetYTitle( "#gamma^{EmbGEN} p_{T}" );
  histos2D_[ "photonRecoIsoEff_EmbGENVsRECOPt" ]->Sumw2();

  histos2D_[ "photonRecoIsoEff_EmbGENVsISOPt" ] = fs->make< TH2D >( "photonRecoIsoEff_EmbGENVsISOPt", "#gamma^{EmbGEN} p_{T} vs. #gamma^{ISO} p_{T}", 101, -10., 1000., 101, -10., 1000.);
  histos2D_[ "photonRecoIsoEff_EmbGENVsISOPt" ]->SetXTitle( "#gamma^{ISO} p_{T}" );
  histos2D_[ "photonRecoIsoEff_EmbGENVsISOPt" ]->SetYTitle( "#gamma^{EmbGEN} p_{T}" );
  histos2D_[ "photonRecoIsoEff_EmbGENVsISOPt" ]->Sumw2();


}

// ------------ method called once each job just after ending the event loop  ------------
void RecoIsoEffG::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void RecoIsoEffG::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoIsoEffG);
