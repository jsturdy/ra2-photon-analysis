// -*- C++ -*-
//
// Package:    Electrons
// Class:      Electrons
// 
/**\class Electrons AddElectronUserData.cc ZInvisibleBkgds/Photons/plugins/AddElectronUserData.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: AddElectronUserData.cc,v 1.1 2012/07/20 11:34:12 sturdy Exp $
//
//


#include "ZInvisibleBkgds/Photons/interface/AddElectronUserData.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
AddElectronUserData::AddElectronUserData(const edm::ParameterSet& pset) :
  debug_         ( pset.getParameter< bool >( "debug" ) ),
  debugString_   ( pset.getParameter< std::string >( "debugString" ) ),
  electronLabel_ ( pset.getParameter< edm::InputTag >( "electronLabel" ) ),
  floatLabels_   ( pset.getParameter< std::vector<edm::InputTag> >( "floatLabels" ) ),
  floatNames_    ( pset.getParameter< std::vector<std::string> >  ( "floatNames" ) ),
  useUserData_   ( pset.exists("userData") ),
  addConversions_( pset.getParameter< bool >  ( "embedConversionInfo" ) )
{
  using namespace pat;
  // produces vector of electrons
  produces<std::vector<pat::Electron> >();
  if ( useUserData_ ) {
    userDataHelper_ = pat::PATUserDataHelper<pat::Electron>(pset.getParameter<edm::ParameterSet>("userData"));
  }
  if ( addConversions_ ) {
    conversionsLabel_ = pset.getParameter< edm::InputTag >( "conversionsLabel" );
    beamspotLabel_    = pset.getParameter< edm::InputTag >( "beamspotLabel" );
  }  
  
}


AddElectronUserData::~AddElectronUserData()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void AddElectronUserData::produce(edm::Event& ev, const edm::EventSetup& es)
{
  using namespace pat;
  //get the electrons
  edm::Handle<edm::View<pat::Electron> > electrons;
  ev.getByLabel(electronLabel_,electrons);

  edm::Handle<reco::ConversionCollection> conversions;
  edm::Handle<reco::BeamSpot> bs;
  //const reco::BeamSpot &beamspot = *bs.product();

  if (addConversions_) {
    ev.getByLabel(conversionsLabel_,conversions);
    ev.getByLabel(beamspotLabel_, bs);
  }
  if ( floatLabels_.size()!=floatNames_.size()) {
    std::cout<<"mismatch between supplied floats and names"<<std::endl;
    std::cout<<"floatLabels_.size()="<<floatLabels_.size()<<std::endl;
    std::cout<<"floatNames_.size()="<<floatNames_.size()  <<std::endl;
    return;
  }
  std::vector<double> myFloats;
  myFloats.reserve(floatLabels_.size());
  std::vector<edm::InputTag>::const_iterator tag = floatLabels_.begin();
  for (; tag != floatLabels_.end(); ++tag) {
    edm::Handle<double> floatVal;
    ev.getByLabel(*tag,floatVal);
    myFloats.push_back(*floatVal);
  }

  std::vector<pat::Electron> * PATElectrons = new std::vector<pat::Electron>(); 
  for (edm::View<pat::Electron>::const_iterator itElectron = electrons->begin(); itElectron != electrons->end(); itElectron++) {

    pat::Electron aElectron(*itElectron);
    
    if ( useUserData_ ) {
      userDataHelper_.add( aElectron, ev, es );
    }
    std::vector<std::string>::const_iterator name = floatNames_.begin();
    for (std::vector<double>::const_iterator val = myFloats.begin(); val != myFloats.end(); ++val) {
      aElectron.addUserFloat(*name,*val);
      ++name;
    }

    if (addConversions_) {
      const reco::BeamSpot &beamspot = *bs.product();
      bool passconversionveto = !ConversionTools::hasMatchedConversion(aElectron,conversions,beamspot.position());
      aElectron.addUserInt("passConvVeto",(int)passconversionveto);
    }

    PATElectrons->push_back(aElectron);
  }

  // sort Electrons in ET
  std::sort(PATElectrons->begin(), PATElectrons->end(), eTComparator_);
  // put genEvt object in Event
  std::auto_ptr<std::vector<pat::Electron> > myElectrons(PATElectrons);
  ev.put(myElectrons);
}

// ------------ method called once each job just before starting event loop  ------------
void AddElectronUserData::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void AddElectronUserData::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void AddElectronUserData::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AddElectronUserData);
