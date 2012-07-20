// -*- C++ -*-
//
// Package:    Photons
// Class:      Photons
// 
/**\class Photons AcceptanceG.cc ZInvisibleBkgds/Photons/plugins/AcceptanceG.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: AcceptanceG.cc,v 1.1 2012/05/16 20:25:39 sturdy Exp $
//
//


#include "ZInvisibleBkgds/Photons/interface/AcceptanceG.h"

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
AcceptanceG::AcceptanceG(const edm::ParameterSet& pset) :
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

  cleanJetHTLabel_ ( pset.getParameter< edm::InputTag >( "cleanJetHTLabel" ) ),
  cleanJetMHTLabel_( pset.getParameter< edm::InputTag >( "cleanJetMHTLabel" ) ),
  cleanHTLabel_    ( pset.getParameter< edm::InputTag >( "cleanHTLabel" ) ),
  cleanMHTLabel_   ( pset.getParameter< edm::InputTag >( "cleanMHTLabel" ) ),

  photonLabel_  ( pset.getParameter< edm::InputTag >( "photonLabel" ) ),
  photonMinPt_  ( pset.getParameter< double >( "photonMinPt" ) ),
  photonEBMaxEta_( pset.getParameter< double >( "photonEBMaxEta" ) ),
  photonEEMinEta_( pset.getParameter< double >( "photonEEMinEta" ) ),
  photonEEMaxEta_( pset.getParameter< double >( "photonEEMaxEta" ) )
{
  photonPtBins_  = pset.getParameter< std::vector<double> >( "photonPtBins" );
  photonEtaBins_ = pset.getParameter< std::vector<double> >( "photonEtaBins" );
  htBins_        = pset.getParameter< std::vector<double> >( "htBins" );
  mhtBins_       = pset.getParameter< std::vector<double> >( "mhtBins" );
  
}


AcceptanceG::~AcceptanceG()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void AcceptanceG::produce(edm::Event& ev, const edm::EventSetup& es)
{
  //  using namespace edm;
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
  if (!jetsHT.isValid()) {
    std::cout<<"unable to find pat::Jet HT jetCollection "<<jetHTLabel_<<std::endl;
    return;
  }

  edm::Handle<edm::View<pat::Jet> > jetsMHT;
  ev.getByLabel(jetMHTLabel_,jetsMHT);
  if (!jetsMHT.isValid()) {
    std::cout<<"unable to find MHT jetCollection "<<jetMHTLabel_<<std::endl;
    return;
  }
  //get the ht
  edm::Handle<double> ht;
  ev.getByLabel(htLabel_,ht);
  //get the mht
  edm::Handle<edm::View<reco::MET> > mht;
  ev.getByLabel(mhtLabel_,mht);

  //get the cleaned jets
  edm::Handle<edm::View<pat::Jet> > cleanJetsHT;
  ev.getByLabel(cleanJetHTLabel_,cleanJetsHT);
  if (!cleanJetsHT.isValid()) {
    std::cout<<"unable to find pat::Jet HT jetCollection "<<cleanJetHTLabel_<<std::endl;
    return;
  }

  edm::Handle<edm::View<pat::Jet> > cleanJetsMHT;
  ev.getByLabel(cleanJetMHTLabel_,cleanJetsMHT);
  if (!cleanJetsMHT.isValid()) {
    std::cout<<"unable to find MHT jetCollection "<<cleanJetMHTLabel_<<std::endl;
    return;
  }
  //get the ht
  edm::Handle<double> cleanHT;
  ev.getByLabel(cleanHTLabel_,cleanHT);
  //get the mht
  edm::Handle<edm::View<reco::MET> > cleanMHT;
  ev.getByLabel(cleanMHTLabel_,cleanMHT);

  //get the photons
  edm::Handle<edm::View<pat::Photon> > photons;
  ev.getByLabel(photonLabel_,photons);

  ///fill histograms
  double htFillValue  = *ht;
  double mhtFillValue = (*mht)[0].pt();
  //bool htIsOflow  = false;
  //bool mhtIsOflow = false;
  if (*ht > htBins_.at(htBins_.size()-1))
    //htIsOflow = true;
    htFillValue = htBins_.at(htBins_.size()-1)-0.01;
  if ((*mht)[0].pt() > mhtBins_.at(mhtBins_.size()-1))
    //mhtIsOflow = true;
    mhtFillValue = mhtBins_.at(mhtBins_.size()-1)-0.01;
  double cleanHTFillValue  = *cleanHT;
  double cleanMHTFillValue = (*cleanMHT)[0].pt();
  if (*cleanHT > htBins_.at(htBins_.size()-1))
    cleanHTFillValue = htBins_.at(htBins_.size()-1)-0.01;
  if ((*cleanMHT)[0].pt() > mhtBins_.at(mhtBins_.size()-1))
    cleanMHTFillValue = mhtBins_.at(mhtBins_.size()-1)-0.01;
  
  if (*ht > htBins_.at(htBins_.size()-1) || (*mht)[0].pt() > mhtBins_.at(mhtBins_.size()-1) ||
      *cleanHT > htBins_.at(htBins_.size()-1) || (*cleanMHT)[0].pt() > mhtBins_.at(mhtBins_.size()-1))
    {
      std::cout<<"ht jets::"<<jetsHT->size()<<" -- mht jets::"<<jetsMHT->size();
      std::cout<<"  ||  clean ht jets::"<<cleanJetsHT->size()<<" -- clean mht jets::"<<cleanJetsMHT->size()<<std::endl;
      std::cout<<"ht         ::"<<*ht<<" -- mht         ::"<<(*mht)[0].pt();
      std::cout<<"  ||  cleanHT         ::"<<*cleanHT<<" -- cleanMHT         ::"<<(*cleanMHT)[0].pt()<<std::endl;
      std::cout<<"htFillValue::"<<htFillValue<<" -- mhtFillValue::"<<mhtFillValue;
      std::cout<<"  ||  cleanHTFillValue::"<<cleanHTFillValue<<" -- cleanMHTFillValue::"<<cleanMHTFillValue<<std::endl;
    }
  //loop over gen particle collection
  reco::GenParticleCollection::const_iterator gen = gens->begin();
  if (debug_)
    std::cout<<debugString_<<"::looping over gen particles, searching for Ws"<<std::endl;
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
    
    //if (gen->numberOfMothers()>0)
    //  histos1D_["photon_Mother"] ->Fill(gen->mother(0)->pdgId() , eventWeightPU);
    //else 
    //  histos1D_["photon_Mother"] ->Fill(0. , eventWeightPU);
    
    histos1D_["photon_Eta"         ]->Fill(gen->eta() , eventWeightPU);
    histos1D_["photon_Pt"          ]->Fill(gen->pt()  , eventWeightPU);
    histos2D_["photon_EtaVsPt"     ]->Fill(gen->pt(),  gen->eta(),      eventWeightPU);
    //histos2D_["photon_HTJetsVsEta" ]->Fill(gen->eta(), jetsHT->size(),  eventWeightPU);
    //histos2D_["photon_HTJetsVsPt"  ]->Fill(gen->pt(),  jetsHT->size(),  eventWeightPU);
    //histos2D_["photon_MHTJetsVsEta"]->Fill(gen->eta(), jetsMHT->size(), eventWeightPU);
    //histos2D_["photon_MHTJetsVsPt" ]->Fill(gen->pt(),  jetsMHT->size(), eventWeightPU);
    //histos2D_["photon_MHTVsHT"     ]->Fill( htFillValue, mhtFillValue,      eventWeightPU);
    histos2D_["photon_MHTVsHT"     ]->Fill( htFillValue, mhtFillValue,      eventWeightPU);
    
    //histos2D_["photon_CleanMHTVsHT"]->Fill( *cleanHT, (*cleanMHT)[0].pt(), eventWeightPU);
    histos2D_["photon_CleanMHTVsHT"]->Fill( cleanHTFillValue, cleanMHTFillValue, eventWeightPU);

    ////passing acceptance
    bool acceptEta = false;
    bool acceptPt = false;
    if (fabs(gen->eta()) < photonEBMaxEta_ || 
	(fabs(gen->eta()) > photonEEMinEta_ && fabs(gen->eta()) < photonEEMaxEta_ ) )
      acceptEta = true;
    if (gen->pt() > photonMinPt_)
      acceptPt = true;

    if (acceptPt) {
      histos1D_["photonAcceptPt_Eta"    ]->Fill(gen->eta(),                 eventWeightPU);
      histos1D_["photonAcceptPt_Pt"     ]->Fill(gen->pt(),                  eventWeightPU);
      histos2D_["photonAcceptPt_EtaVsPt"]->Fill(gen->pt(), gen->eta(),      eventWeightPU);
      histos2D_["photonAcceptPt_MHTVsHT"     ]->Fill( htFillValue,      mhtFillValue,      eventWeightPU);
      histos2D_["photonAcceptPt_CleanMHTVsHT"]->Fill( cleanHTFillValue, cleanMHTFillValue, eventWeightPU);
    }

    if (acceptEta) {
      histos1D_["photonAcceptEta_Eta"    ]->Fill(gen->eta(),                 eventWeightPU);
      histos1D_["photonAcceptEta_Pt"     ]->Fill(gen->pt(),                  eventWeightPU);
      histos2D_["photonAcceptEta_EtaVsPt"]->Fill(gen->pt(), gen->eta(),      eventWeightPU);
      histos2D_["photonAcceptEta_MHTVsHT"     ]->Fill( htFillValue,      mhtFillValue,      eventWeightPU);
      histos2D_["photonAcceptEta_CleanMHTVsHT"]->Fill( cleanHTFillValue, cleanMHTFillValue, eventWeightPU);
    }
    
    if (acceptEta && acceptPt) {
      //if (gen->numberOfMothers()>0)
      //histos1D_["photonAccepted_Mother"] ->Fill(gen->mother(0)->pdgId() , eventWeightPU);
      //else 
      //histos1D_["photonAccepted_Mother"] ->Fill(0. , eventWeightPU);
      
      histos1D_["photonAccepted_Eta"    ]->Fill(gen->eta(),                  eventWeightPU);
      histos1D_["photonAccepted_Pt"     ]->Fill(gen->pt(),                   eventWeightPU);
      histos2D_["photonAccepted_EtaVsPt"]->Fill(gen->pt(),  gen->eta(),      eventWeightPU);
      //histos2D_["photonAccepted_HTJetsVsEta" ]->Fill(gen->eta(), jetsHT->size(),  eventWeightPU);
      //histos2D_["photonAccepted_HTJetsVsPt"  ]->Fill(gen->pt(),  jetsHT->size(),  eventWeightPU);
      //histos2D_["photonAccepted_MHTJetsVsEta"]->Fill(gen->eta(), jetsMHT->size(), eventWeightPU);
      //histos2D_["photonAccepted_MHTJetsVsPt" ]->Fill(gen->pt(),  jetsMHT->size(), eventWeightPU);
      histos2D_["photonAccepted_MHTVsHT"     ]->Fill( htFillValue,      mhtFillValue,      eventWeightPU);
      histos2D_["photonAccepted_CleanMHTVsHT"]->Fill( cleanHTFillValue, cleanMHTFillValue, eventWeightPU);
    }
  }//end loop over gen particles
}

// ------------ method called once each job just before starting event loop  ------------
void AcceptanceG::beginJob()
{
  edm::Service< TFileService > fs;
  
  //char title[128];
  //sprintf(title,"# of events passing HLT_%s_v*",unbiasedTrigger_.c_str());
  
  //histos1D_[ "photonAcceptPt_Eta" ] = fs->make< TH1D >( "photonAcceptPt_Eta", "#eta for #gamma passing p_{T}", photonEtaBins_.size()-1, photonEtaBins_.data());
  //histos1D_[ "photonAcceptPt_Eta" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  //histos1D_[ "photonAcceptPt_Eta" ]->SetYTitle( "# of events" );
  //histos1D_[ "photonAcceptPt_Eta" ]->Sumw2();
  //
  //histos1D_[ "photonAcceptEta_Pt" ] = fs->make< TH1D >( "photonAcceptEta_Pt", "p_T for #gamma passing #eta", photonPtBins_.size()-1, photonPtBins_.data());
  //histos1D_[ "photonAcceptEta_Pt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  //histos1D_[ "photonAcceptEta_Pt" ]->SetYTitle( "# of events" );
  //histos1D_[ "photonAcceptEta_Pt" ]->Sumw2();

  histos1D_[ "photon_Eta" ] = fs->make< TH1D >( "photon_Eta", "#eta", photonEtaBins_.size()-1, photonEtaBins_.data());
  histos1D_[ "photon_Eta" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  histos1D_[ "photon_Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "photon_Eta" ]->Sumw2();

  histos1D_[ "photon_Pt" ] = fs->make< TH1D >( "photon_Pt", "p_T", photonPtBins_.size()-1, photonPtBins_.data());
  histos1D_[ "photon_Pt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  histos1D_[ "photon_Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "photon_Pt" ]->Sumw2();

  //eta vs. pt
  histos2D_[ "photon_EtaVsPt" ] = fs->make< TH2D >( "photon_EtaVsPt", "#eta vs. p_T", photonPtBins_.size()-1, photonPtBins_.data(), photonEtaBins_.size()-1, photonEtaBins_.data());
  histos2D_[ "photon_EtaVsPt" ]->SetYTitle( "#eta^{#gamma^{GEN}}" );
  histos2D_[ "photon_EtaVsPt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  histos2D_[ "photon_EtaVsPt" ]->Sumw2();

  //mht vs. ht
  histos2D_[ "photon_MHTVsHT" ] = fs->make< TH2D >( "photon_MHTVsHT", "Event #slashH_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data());
  histos2D_[ "photon_MHTVsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos2D_[ "photon_MHTVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "photon_MHTVsHT" ]->Sumw2();

  //mht vs. ht  (cleaned for reco/iso photon)
  histos2D_[ "photon_CleanMHTVsHT" ] = fs->make< TH2D >( "photon_CleanMHTVsHT", "Event #slashH_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data());
  histos2D_[ "photon_CleanMHTVsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos2D_[ "photon_CleanMHTVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "photon_CleanMHTVsHT" ]->Sumw2();

  ////Passing Pt Cuts
  histos1D_[ "photonAcceptPt_Eta" ] = fs->make< TH1D >( "photonAcceptPt_Eta", "#eta passing #gamma p_{T} cut", photonEtaBins_.size()-1, photonEtaBins_.data());
  histos1D_[ "photonAcceptPt_Eta" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  histos1D_[ "photonAcceptPt_Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "photonAcceptPt_Eta" ]->Sumw2();

  histos1D_[ "photonAcceptPt_Pt" ] = fs->make< TH1D >( "photonAcceptPt_Pt", "p_T passing #gamma p_{T} cut", photonPtBins_.size()-1, photonPtBins_.data());
  histos1D_[ "photonAcceptPt_Pt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  histos1D_[ "photonAcceptPt_Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "photonAcceptPt_Pt" ]->Sumw2();

  //eta vs. pt
  histos2D_[ "photonAcceptPt_EtaVsPt" ] = fs->make< TH2D >( "photonAcceptPt_EtaVsPt", "#eta vs. p_T passing #gamma p_{T} cut", photonPtBins_.size()-1, photonPtBins_.data(), photonEtaBins_.size()-1, photonEtaBins_.data());
  histos2D_[ "photonAcceptPt_EtaVsPt" ]->SetYTitle( "#eta^{#gamma^{GEN}}" );
  histos2D_[ "photonAcceptPt_EtaVsPt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  histos2D_[ "photonAcceptPt_EtaVsPt" ]->Sumw2();

  //mht vs. ht  (cleaned for reco/iso photon)
  histos2D_[ "photonAcceptPt_CleanMHTVsHT" ] = fs->make< TH2D >( "photonAcceptPt_CleanMHTVsHT", "Event #slashH_{T} vs. Event H_{T} passing p_{T} cut", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data());
  histos2D_[ "photonAcceptPt_CleanMHTVsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos2D_[ "photonAcceptPt_CleanMHTVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "photonAcceptPt_CleanMHTVsHT" ]->Sumw2();

  //mht vs. ht
  histos2D_[ "photonAcceptPt_MHTVsHT" ] = fs->make< TH2D >( "photonAcceptPt_MHTVsHT", "Event #slashH_{T} vs. Event H_{T} passing p_{T} cut", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data());
  histos2D_[ "photonAcceptPt_MHTVsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos2D_[ "photonAcceptPt_MHTVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "photonAcceptPt_MHTVsHT" ]->Sumw2();

  ////Passing Eta Cuts
  histos1D_[ "photonAcceptEta_Eta" ] = fs->make< TH1D >( "photonAcceptEta_Eta", "#eta passing #gamma #eta cut", photonEtaBins_.size()-1, photonEtaBins_.data());
  histos1D_[ "photonAcceptEta_Eta" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  histos1D_[ "photonAcceptEta_Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "photonAcceptEta_Eta" ]->Sumw2();

  histos1D_[ "photonAcceptEta_Pt" ] = fs->make< TH1D >( "photonAcceptEta_Pt", "p_T passing #gamma #eta cut", photonPtBins_.size()-1, photonPtBins_.data());
  histos1D_[ "photonAcceptEta_Pt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  histos1D_[ "photonAcceptEta_Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "photonAcceptEta_Pt" ]->Sumw2();

  //eta vs. pt
  histos2D_[ "photonAcceptEta_EtaVsPt" ] = fs->make< TH2D >( "photonAcceptEta_EtaVsPt", "#eta vs. p_T passing #gamma #eta cut", photonPtBins_.size()-1, photonPtBins_.data(), photonEtaBins_.size()-1, photonEtaBins_.data());
  histos2D_[ "photonAcceptEta_EtaVsPt" ]->SetYTitle( "#eta^{#gamma^{GEN}}" );
  histos2D_[ "photonAcceptEta_EtaVsPt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  histos2D_[ "photonAcceptEta_EtaVsPt" ]->Sumw2();

  //mht vs. ht
  histos2D_[ "photonAcceptEta_MHTVsHT" ] = fs->make< TH2D >( "photonAcceptEta_MHTVsHT", "Event #slashH_{T} vs. Event H_{T} passing #eta cut", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data());
  histos2D_[ "photonAcceptEta_MHTVsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos2D_[ "photonAcceptEta_MHTVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "photonAcceptEta_MHTVsHT" ]->Sumw2();

  //mht vs. ht  (cleaned for reco/iso photon)
  histos2D_[ "photonAcceptEta_CleanMHTVsHT" ] = fs->make< TH2D >( "photonAcceptEta_CleanMHTVsHT", "Event #slashH_{T} vs. Event H_{T} passing #eta cut", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data());
  histos2D_[ "photonAcceptEta_CleanMHTVsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos2D_[ "photonAcceptEta_CleanMHTVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "photonAcceptEta_CleanMHTVsHT" ]->Sumw2();

  //histos1D_[ "photon_Mother" ] = fs->make< TH1D >( "photon_Mother", "pdgID", 61, -30.5, 30.5);
  //histos1D_[ "photon_Mother" ]->SetXTitle( "#gamma^{GEN} Mother" );
  //histos1D_[ "photon_Mother" ]->SetYTitle( "# of events" );
  //histos1D_[ "photon_Mother" ]->Sumw2();
  //
  ////vs. number of jets
  //histos2D_[ "photon_HTJetsVsEta" ] = fs->make< TH2D >( "photon_HTJetsVsEta", "#eta vs. #Jets (p_{T} > 50 GeV, |#eta| < 2.5)", photonEtaBins_.size()-1, photonEtaBins_.data(), 21, -0.5, 20.5);
  //histos2D_[ "photon_HTJetsVsEta" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  //histos2D_[ "photon_HTJetsVsEta" ]->SetYTitle( "# of jets" );
  //histos2D_[ "photon_HTJetsVsEta" ]->Sumw2();
  //
  //histos2D_[ "photon_HTJetsVsPt" ] = fs->make< TH2D >( "photon_HTJetsVsPt", "p_T vs. #Jets (p_{T} > 50 GeV, |#eta| < 2.5)", photonPtBins_.size()-1, photonPtBins_.data(), 21, -0.5, 20.5);
  //histos2D_[ "photon_HTJetsVsPt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  //histos2D_[ "photon_HTJetsVsPt" ]->SetYTitle( "# of jets" );
  //histos2D_[ "photon_HTJetsVsPt" ]->Sumw2();
  //
  //histos2D_[ "photon_MHTJetsVsEta" ] = fs->make< TH2D >( "photon_MHTJetsVsEta", "#eta vs. #Jets (p_{T} > 30 GeV, |#eta| < 5.0)", photonEtaBins_.size()-1, photonEtaBins_.data(), 21, -0.5, 20.5);
  //histos2D_[ "photon_MHTJetsVsEta" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  //histos2D_[ "photon_MHTJetsVsEta" ]->SetYTitle( "# of jets" );
  //histos2D_[ "photon_MHTJetsVsEta" ]->Sumw2();
  //
  //histos2D_[ "photon_MHTJetsVsPt" ] = fs->make< TH2D >( "photon_MHTJetsVsPt", "p_T vs. #Jets (p_{T} > 30 GeV, |#eta| < 5.0)", photonPtBins_.size()-1, photonPtBins_.data(), 21, -0.5, 20.5);
  //histos2D_[ "photon_MHTJetsVsPt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  //histos2D_[ "photon_MHTJetsVsPt" ]->SetYTitle( "# of jets" );
  //histos2D_[ "photon_MHTJetsVsPt" ]->Sumw2();
  //

  ////Accepted photons
  histos1D_[ "photonAccepted_Eta" ] = fs->make< TH1D >( "photonAccepted_Eta", "#eta", photonEtaBins_.size()-1, photonEtaBins_.data());
  histos1D_[ "photonAccepted_Eta" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  histos1D_[ "photonAccepted_Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "photonAccepted_Eta" ]->Sumw2();

  histos1D_[ "photonAccepted_Pt" ] = fs->make< TH1D >( "photonAccepted_Pt", "p_T", photonPtBins_.size()-1, photonPtBins_.data());
  histos1D_[ "photonAccepted_Pt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  histos1D_[ "photonAccepted_Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "photonAccepted_Pt" ]->Sumw2();

  //eta vs. pt
  histos2D_[ "photonAccepted_EtaVsPt" ] = fs->make< TH2D >( "photonAccepted_EtaVsPt", "#eta vs. p_T", photonPtBins_.size()-1, photonPtBins_.data(), photonEtaBins_.size()-1, photonEtaBins_.data());
  histos2D_[ "photonAccepted_EtaVsPt" ]->SetYTitle( "#eta^{#gamma^{GEN}}" );
  histos2D_[ "photonAccepted_EtaVsPt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  histos2D_[ "photonAccepted_EtaVsPt" ]->Sumw2();

  //mht vs. ht
  histos2D_[ "photonAccepted_MHTVsHT" ] = fs->make< TH2D >( "photonAccepted_MHTVsHT", "Event #slashH_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data());
  histos2D_[ "photonAccepted_MHTVsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos2D_[ "photonAccepted_MHTVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "photonAccepted_MHTVsHT" ]->Sumw2();

  //mht vs. ht  (cleaned for reco/iso photon)
  histos2D_[ "photonAccepted_CleanMHTVsHT" ] = fs->make< TH2D >( "photonAccepted_CleanMHTVsHT", "Event #slashH_{T} vs. Event H_{T}", htBins_.size()-1, htBins_.data(), mhtBins_.size()-1, mhtBins_.data());
  histos2D_[ "photonAccepted_CleanMHTVsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos2D_[ "photonAccepted_CleanMHTVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "photonAccepted_CleanMHTVsHT" ]->Sumw2();

  //histos1D_[ "photonAccepted_Mother" ] = fs->make< TH1D >( "photonAccepted_Mother", "pdgID", 61, -30.5, 30.5);
  //histos1D_[ "photonAccepted_Mother" ]->SetXTitle( "#gamma^{GEN} Mother" );
  //histos1D_[ "photonAccepted_Mother" ]->SetYTitle( "# of events" );
  //histos1D_[ "photonAccepted_Mother" ]->Sumw2();
  //
  ////vs. number of jets
  //histos2D_[ "photonAccepted_HTJetsVsEta" ] = fs->make< TH2D >( "photonAccepted_HTJetsVsEta", "#eta vs. #Jets (p_{T} > 50 GeV, |#eta| < 2.5)", photonEtaBins_.size()-1, photonEtaBins_.data(), 21, -0.5, 20.5);
  //histos2D_[ "photonAccepted_HTJetsVsEta" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  //histos2D_[ "photonAccepted_HTJetsVsEta" ]->SetYTitle( "# of jets" );
  //histos2D_[ "photonAccepted_HTJetsVsEta" ]->Sumw2();
  //
  //histos2D_[ "photonAccepted_HTJetsVsPt" ] = fs->make< TH2D >( "photonAccepted_HTJetsVsPt", "p_T vs. #Jets (p_{T} > 50 GeV, |#eta| < 2.5)", photonPtBins_.size()-1, photonPtBins_.data(), 21, -0.5, 20.5);
  //histos2D_[ "photonAccepted_HTJetsVsPt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  //histos2D_[ "photonAccepted_HTJetsVsPt" ]->SetYTitle( "# of jets" );
  //histos2D_[ "photonAccepted_HTJetsVsPt" ]->Sumw2();
  //
  //histos2D_[ "photonAccepted_MHTJetsVsEta" ] = fs->make< TH2D >( "photonAccepted_MHTJetsVsEta", "#eta vs. #Jets (p_{T} > 30 GeV, |#eta| < 5.0)", photonEtaBins_.size()-1, photonEtaBins_.data(), 21, -0.5, 20.5);
  //histos2D_[ "photonAccepted_MHTJetsVsEta" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  //histos2D_[ "photonAccepted_MHTJetsVsEta" ]->SetYTitle( "# of jets" );
  //histos2D_[ "photonAccepted_MHTJetsVsEta" ]->Sumw2();
  //
  //histos2D_[ "photonAccepted_MHTJetsVsPt" ] = fs->make< TH2D >( "photonAccepted_MHTJetsVsPt", "p_T vs. #Jets (p_{T} > 30 GeV, |#eta| < 5.0)", photonPtBins_.size()-1, photonPtBins_.data(), 21, -0.5, 20.5);
  //histos2D_[ "photonAccepted_MHTJetsVsPt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  //histos2D_[ "photonAccepted_MHTJetsVsPt" ]->SetYTitle( "# of jets" );
  //histos2D_[ "photonAccepted_MHTJetsVsPt" ]->Sumw2();
  //

}

// ------------ method called once each job just after ending the event loop  ------------
void AcceptanceG::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void AcceptanceG::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AcceptanceG);
