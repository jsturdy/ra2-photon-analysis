// -*- C++ -*-
//
// Package:    Photons
// Class:      Photons
// 
/**\class Photons GenBosonComparison.cc ZInvisibleBkgds/Photons/plugins/GenBosonComparison.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: GenBosonComparison.cc,v 1.1 2012/05/16 20:25:39 sturdy Exp $
//
//


#include "ZInvisibleBkgds/Photons/interface/GenBosonComparison.h"

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
GenBosonComparison::GenBosonComparison(const edm::ParameterSet& pset) :
  debug_      ( pset.getParameter< bool >( "debug" ) ),
  debugString_( pset.getParameter< std::string >( "debugString" ) ),
  genLabel_   ( pset.getParameter< edm::InputTag >( "genLabel" ) ),
  genStatus_  ( pset.getParameter< int >( "genStatus" ) ),
  genPDGId_   ( pset.getParameter< int >( "genPDGId" ) ),

  doPUReweight_ ( pset.getParameter< bool >( "doPUReweight" ) ),
  puWeightLabel_( pset.getParameter< edm::InputTag >( "puWeight" ) ),

  jetHTLabel_ ( pset.getParameter< edm::InputTag >( "jetHTLabel" ) ),
  jetMHTLabel_( pset.getParameter< edm::InputTag >( "jetMHTLabel" ) ),
  htLabel_    ( pset.getParameter< edm::InputTag >( "htLabel" ) ),
  mhtLabel_   ( pset.getParameter< edm::InputTag >( "mhtLabel" ) )
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


GenBosonComparison::~GenBosonComparison()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void GenBosonComparison::produce(edm::Event& ev, const edm::EventSetup& es)
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

  /*
    reconstruction/isolation efficiency in each bin found by matching the 
    gen photon to a reco photon passing isolation cuts
  */

  histos1D_["genBoson_NJets_Pt30"     ]->Fill( jetsMHT->size(), eventWeightPU);
  histos1D_["genBoson_NJets_Pt50Eta25"]->Fill( jetsHT->size(),  eventWeightPU);
  histos1D_["genBoson_HT"             ]->Fill( *ht,             eventWeightPU);
  histos1D_["genBoson_MHT"            ]->Fill( (*mht)[0].pt(),  eventWeightPU);
    
  //loop over gen particle collection
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

    if (gen->pdgId() == genPDGId_ && gen->status() == genStatus_) {
      ////////////// do i need to do some matching in the case of the photon?
      // remove reco photon from jet list, new ht/mht/dphi
      // or can i just use the cleaned reco collections and require at least one 
      // reco/iso photon?
      // i think the latter, as we will account for reco/iso eff and acceptance later
      histos1D_["genBoson_Pt" ]->Fill( gen->pt() , eventWeightPU);
      histos1D_["genBoson_Eta"]->Fill( gen->eta(), eventWeightPU);
      histos1D_["genBoson_Phi"]->Fill( gen->phi(), eventWeightPU);
      
      histos2D_["genBoson_PtVsHT" ]->Fill( *ht, gen->pt(),      eventWeightPU);
      histos2D_["genBoson_MHTVsHT"]->Fill( *ht, (*mht)[0].pt(), eventWeightPU);
    }
  }//end loop over gen particles
}

// ------------ method called once each job just before starting event loop  ------------
void GenBosonComparison::beginJob()
{
  edm::Service< TFileService > fs;
  
  //char title[128];
  //sprintf(title,"# of events passing HLT_%s_v*",unbiasedTrigger_.c_str());
  histos1D_[ "genBoson_NJets_Pt30" ] = fs->make< TH1D >( "genBoson_NJets_Pt30", "#Jets (p_{T} > 30 GeV, |#eta| < 5.0)", 25, -0.5,24.5);
  histos1D_[ "genBoson_NJets_Pt30" ]->SetXTitle( "#Jets" );
  histos1D_[ "genBoson_NJets_Pt30" ]->SetYTitle( "# of events" );
  histos1D_[ "genBoson_NJets_Pt30" ]->Sumw2();

  histos1D_[ "genBoson_NJets_Pt50Eta25" ] = fs->make< TH1D >( "genBoson_NJets_Pt50Eta25", "#Jets (p_{T} > 50 GeV, |#eta| < 2.5)", 25, -0.5,24.5);
  histos1D_[ "genBoson_NJets_Pt50Eta25" ]->SetXTitle( "#Jets" );
  histos1D_[ "genBoson_NJets_Pt50Eta25" ]->SetYTitle( "# of events" );
  histos1D_[ "genBoson_NJets_Pt50Eta25" ]->Sumw2();

  histos1D_[ "genBoson_HT" ] = fs->make< TH1D >( "genBoson_HT", "Event H_{T}", 250, 0., 5000.);
  histos1D_[ "genBoson_HT" ]->SetXTitle( "Event H_{T}" );
  histos1D_[ "genBoson_HT" ]->SetYTitle( "# of events" );
  histos1D_[ "genBoson_HT" ]->Sumw2();

  histos1D_[ "genBoson_MHT" ] = fs->make< TH1D >( "genBoson_MHT", "Event #slashH_{T}", 200, 0., 2000.);
  histos1D_[ "genBoson_MHT" ]->SetXTitle( "Event #slashH_{T}" );
  histos1D_[ "genBoson_MHT" ]->SetYTitle( "# of events" );
  histos1D_[ "genBoson_MHT" ]->Sumw2();

  histos1D_[ "genBoson_Pt" ] = fs->make< TH1D >( "genBoson_Pt", "Gen Boson p_{T}", 100, 0, 1000);
  histos1D_[ "genBoson_Pt" ]->SetXTitle( "Boson p_{T}" );
  histos1D_[ "genBoson_Pt" ]->SetYTitle( "# of events" );
  histos1D_[ "genBoson_Pt" ]->Sumw2();

  histos1D_[ "genBoson_Eta" ] = fs->make< TH1D >( "genBoson_Eta", "Gen Boson #eta", 50, -5., 5.);
  histos1D_[ "genBoson_Eta" ]->SetXTitle( "Boson #eta" );
  histos1D_[ "genBoson_Eta" ]->SetYTitle( "# of events" );
  histos1D_[ "genBoson_Eta" ]->Sumw2();

  histos1D_[ "genBoson_Phi" ] = fs->make< TH1D >( "genBoson_Phi", "Gen Boson #phi", 50, -M_PI, M_PI );
  histos1D_[ "genBoson_Phi" ]->SetXTitle( "Boson #phi" );
  histos1D_[ "genBoson_Phi" ]->SetYTitle( "# of events" );
  histos1D_[ "genBoson_Phi" ]->Sumw2();

  //pt vs. ht
  histos2D_[ "genBoson_PtVsHT" ] = fs->make< TH2D >( "genBoson_PtVsHT", "Gen Boson p_{T} vs. Event H_{T}", 250, 0., 5000., 100, 0, 1000);
  histos2D_[ "genBoson_PtVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "genBoson_PtVsHT" ]->SetYTitle( "Boson p_{T}" );
  histos2D_[ "genBoson_PtVsHT" ]->Sumw2();

  //mht vs. ht
  histos2D_[ "genBoson_MHTVsHT" ] = fs->make< TH2D >( "genBoson_MHTVsHT", "Event #slashH_{T} vs. Event H_{T}", 250, 0., 5000., 200, 0., 2000.);
  histos2D_[ "genBoson_MHTVsHT" ]->SetXTitle( "Event H_{T}" );
  histos2D_[ "genBoson_MHTVsHT" ]->SetYTitle( "Event #slashH_{T}" );
  histos2D_[ "genBoson_MHTVsHT" ]->Sumw2();

}

// ------------ method called once each job just after ending the event loop  ------------
void GenBosonComparison::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenBosonComparison::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenBosonComparison);
