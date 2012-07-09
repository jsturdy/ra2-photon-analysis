// -*- C++ -*-
//
// Package:    Photons
// Class:      Photons
// 
/**\class Photons AddPhotonUserData.cc ZInvisibleBkgds/Photons/plugins/AddPhotonUserData.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: AddPhotonUserData.cc,v 1.1 2012/05/16 20:25:39 sturdy Exp $
//
//


#include "ZInvisibleBkgds/Photons/interface/AddPhotonUserData.h"

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
AddPhotonUserData::AddPhotonUserData(const edm::ParameterSet& pset) :
  debug_       ( pset.getParameter< bool >( "debug" ) ),
  debugString_ ( pset.getParameter< std::string >( "debugString" ) ),
  photonLabel_ ( pset.getParameter< edm::InputTag >( "photonLabel" ) ),
  floatLabels_ ( pset.getParameter< std::vector<edm::InputTag> >( "floatLabels" ) ),
  floatNames_  ( pset.getParameter< std::vector<std::string> >  ( "floatNames" ) ),
  useUserData_ ( pset.exists("userData"))
{
  using namespace pat;
  // produces vector of photons
  produces<std::vector<pat::Photon> >();
  if ( useUserData_ ) {
    userDataHelper_ = PATUserDataHelper<Photon>(pset.getParameter<edm::ParameterSet>("userData"));
  }
  
}


AddPhotonUserData::~AddPhotonUserData()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void AddPhotonUserData::produce(edm::Event& ev, const edm::EventSetup& es)
{
  using namespace pat;
  //get the photons
  edm::Handle<edm::View<pat::Photon> > photons;
  ev.getByLabel(photonLabel_,photons);

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

  std::vector<pat::Photon> * PATPhotons = new std::vector<pat::Photon>(); 
  for (edm::View<pat::Photon>::const_iterator itPhoton = photons->begin(); itPhoton != photons->end(); itPhoton++) {

    pat::Photon aPhoton(*itPhoton);
    
    if ( useUserData_ ) {
      userDataHelper_.add( aPhoton, ev, es );
    }
    std::vector<std::string>::const_iterator name = floatNames_.begin();
    for (std::vector<double>::const_iterator val = myFloats.begin(); val != myFloats.end(); ++val) {
      aPhoton.addUserFloat(*name,*val);
      ++name;
    }
    PATPhotons->push_back(aPhoton);
  }

  // sort Photons in ET
  std::sort(PATPhotons->begin(), PATPhotons->end(), eTComparator_);
  // put genEvt object in Event
  std::auto_ptr<std::vector<pat::Photon> > myPhotons(PATPhotons);
  ev.put(myPhotons);
}

// ------------ method called once each job just before starting event loop  ------------
void AddPhotonUserData::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void AddPhotonUserData::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void AddPhotonUserData::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AddPhotonUserData);
