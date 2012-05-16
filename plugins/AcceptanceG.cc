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
// $Id$
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
AcceptanceG::AcceptanceG(const edm::ParameterSet& iConfig) :
  debug_           ( iConfig.getParameter< bool >( "debug" ) ),
  debugString_     ( iConfig.getParameter< std::string >( "debugString" ) ),
  genLabel_        ( iConfig.getParameter< edm::InputTag >( "genLabel" ) ),
  genStatus_       ( iConfig.getParameter< int >( "genStatus" ) ),

  jetLabel_( iConfig.getParameter< edm::InputTag >( "jetLabel" ) ),


  photonLabel_ ( iConfig.getParameter< edm::InputTag >( "photonLabel" ) ),
  photonMinPt_ ( iConfig.getParameter< double >( "photonMinPt" ) ),
  photonMaxEta_( iConfig.getParameter< double >( "photonMaxEta" ) )
{
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
  using namespace edm;
  //read in the gen particles
  Handle<reco::GenParticleCollection> gens;
  ev.getByLabel(genLabel_,gens);

  //get the jets
  Handle<pat::JetCollection> jets;
  ev.getByLabel(jetLabel_,jets);
  //get the photons
  Handle<pat::PhotonCollection> photons;
  ev.getByLabel(photonLabel_,photons);

  ///fill histograms

  //loop over gen particle collection
  reco::GenParticleCollection::const_iterator gen = gens->begin();
  if (debug_)
    std::cout<<debugString_<<"::looping over gen particles, searching for Ws"<<std::endl;
  for ( ; gen != gens->end(); ++gen) {
    
    if (gen->status() == genStatus_) {
      if (fabs(gen->pdgId()) == 22) {
	if (debug_) {
	  //std::cout<<"pdgID  --  status  --  #dau  --  {dau(ID,St,pT,eta),...}"<<std::endl;
	  std::cout<<"pdgID  --  status  --  #mom  --  {mom(ID,St,#mom,({pdgID,st}),pT,eta),...}"<<std::endl;
	  //std::cout<<gen->pdgId()<<"  --  "<<gen->status()<<"  --  "<<gen->numberOfDaughters()<<"  --  "<<"{";
	  std::cout<<gen->pdgId()<<"  --  "<<gen->status()<<"  --  "<<gen->numberOfMothers()<<"  --  "<<"{";
	  //for (unsigned int dau = 0; dau < gen->numberOfDaughters(); ++dau) 
	  for (unsigned int dau = 0; dau < gen->numberOfMothers(); ++dau) {
	    //std::cout<<"("<<gen->daughter(dau)->pdgId()<<","<<gen->daughter(dau)->status()<<","<<gen->daughter(dau)->pt()<<","<<gen->daughter(dau)->eta()<<") ";
	    std::cout<<"("<<gen->mother(dau)->pdgId()<<","<<gen->mother(dau)->status()<<","<<gen->mother(dau)->numberOfMothers()<<"(";
	    for (unsigned int mom = 0; mom < gen->mother(dau)->numberOfMothers(); ++mom) 
	      std::cout<<"{"<<gen->mother(dau)->mother(mom)->pdgId()<<","<<gen->mother(dau)->mother(mom)->status()<<"},";
	    std::cout<<"),";
	    std::cout<<gen->pt()<<","<<gen->eta()<<") ";
	  }
	  std::cout<<"}"<<std::endl;
	}
	if (gen->status() == genStatus_) {
	  if (fabs(gen->pdgId()) == 22) {
	    histos1D_["genPhotonEta"]->Fill(gen->eta());
	    histos1D_["genPhotonPt"] ->Fill(gen->pt());
	    if (gen->numberOfMothers()>0)
	      histos1D_["genPhotonMother"] ->Fill(gen->mother(0)->pdgId());
	    else 
	      histos1D_["genPhotonMother"] ->Fill(0);
	    histos2D_["genPhotonEtaVSPt"]->Fill(gen->eta(),gen->pt());
	    histos2D_["genPhotonEtaVSnJets"]->Fill(gen->eta(),jets->size());
	    histos2D_["genPhotonPtVSnJets"] ->Fill(gen->pt(), jets->size());
	  }
	}
      }
    }//end status check
  }//end loop over gen particles
}

// ------------ method called once each job just before starting event loop  ------------
void AcceptanceG::beginJob()
{
  edm::Service< TFileService > fs;
  
  //char title[128];
  //sprintf(title,"# of events passing HLT_%s_v*",unbiasedTrigger_.c_str());
  double photonEtaBins[5] = {-5.,-photonMaxEta_,0,photonMaxEta_,5.};
  //double photonEtaBins[5] = {-5.,-photonMaxEta_,-1.566,-1.4442_,0,1.4442,1.566,photonMaxEta_,5.};
  double photonPtBins[3] = {0.,photonMinPt_,100};

  histos1D_[ "genPhotonEta" ] = fs->make< TH1D >( "genPhotonEta", "#eta", 4, photonEtaBins);
  histos1D_[ "genPhotonEta" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  histos1D_[ "genPhotonEta" ]->SetYTitle( "# of events" );
  histos1D_[ "genPhotonEta" ]->Sumw2();

  histos1D_[ "genPhotonPt" ] = fs->make< TH1D >( "genPhotonPt", "p_T", 2, photonPtBins);
  histos1D_[ "genPhotonPt" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  histos1D_[ "genPhotonPt" ]->SetYTitle( "# of events" );
  histos1D_[ "genPhotonPt" ]->Sumw2();

  histos1D_[ "genPhotonMother" ] = fs->make< TH1D >( "genPhotonMother", "pdgID", 60, -30,30);
  histos1D_[ "genPhotonMother" ]->SetXTitle( "#gamma^{GEN} Mother" );
  histos1D_[ "genPhotonMother" ]->SetYTitle( "# of events" );
  histos1D_[ "genPhotonMother" ]->Sumw2();

  //vs. number of jets
  histos2D_[ "genPhotonEtaVSnJets" ] = fs->make< TH2D >( "genPhotonEtaVSnJets", "#eta vs. #Jets", 4, photonEtaBins, 15, 0, 15);
  histos2D_[ "genPhotonEtaVSnJets" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  histos2D_[ "genPhotonEtaVSnJets" ]->SetYTitle( "# of jets" );
  histos2D_[ "genPhotonEtaVSnJets" ]->Sumw2();

  histos2D_[ "genPhotonPtVSnJets" ] = fs->make< TH2D >( "genPhotonPtVSnJets", "p_T vs. #Jets", 2, photonPtBins, 15, 0, 15);
  histos2D_[ "genPhotonPtVSnJets" ]->SetXTitle( "p_{T}^{#gamma^{GEN}}" );
  histos2D_[ "genPhotonPtVSnJets" ]->SetYTitle( "# of jets" );
  histos2D_[ "genPhotonPtVSnJets" ]->Sumw2();

  //eta vs. pt
  histos2D_[ "genPhotonEtaVSPt" ] = fs->make< TH2D >( "genPhotonEtaVSPt", "#eta vs. p_T", 4, photonEtaBins, 2, photonPtBins);
  histos2D_[ "genPhotonEtaVSPt" ]->SetXTitle( "#eta^{#gamma^{GEN}}" );
  histos2D_[ "genPhotonEtaVSPt" ]->SetYTitle( "p_{T}^{#gamma^{GEN}}" );
  histos2D_[ "genPhotonEtaVSPt" ]->Sumw2();
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
