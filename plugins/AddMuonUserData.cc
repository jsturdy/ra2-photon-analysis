// -*- C++ -*-
//
// Package:    Muons
// Class:      Muons
// 
/**\class Muons AddMuonUserData.cc ZInvisibleBkgds/Photons/plugins/AddMuonUserData.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: AddMuonUserData.cc,v 1.1 2012/07/20 11:34:12 sturdy Exp $
//
//


#include "ZInvisibleBkgds/Photons/interface/AddMuonUserData.h"

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
AddMuonUserData::AddMuonUserData(const edm::ParameterSet& pset) :
  debug_         ( pset.getParameter< bool >( "debug" ) ),
  debugString_   ( pset.getParameter< std::string >( "debugString" ) ),
  muonLabel_     ( pset.getParameter< edm::InputTag >( "muonLabel" ) ),
  floatLabels_   ( pset.getParameter< std::vector<edm::InputTag> >( "floatLabels" ) ),
  floatNames_    ( pset.getParameter< std::vector<std::string> >  ( "floatNames" ) ),
  useUserData_   ( pset.exists("userData") )
{
  using namespace pat;
  // produces vector of muons
  produces<std::vector<pat::Muon> >();
  if ( useUserData_ ) {
    userDataHelper_ = PATUserDataHelper<Muon>(pset.getParameter<edm::ParameterSet>("userData"));
  }
}


AddMuonUserData::~AddMuonUserData()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void AddMuonUserData::produce(edm::Event& ev, const edm::EventSetup& es)
{
  using namespace pat;
  //get the muons
  edm::Handle<edm::View<pat::Muon> > muons;
  ev.getByLabel(muonLabel_,muons);

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

  std::vector<pat::Muon> * PATMuons = new std::vector<pat::Muon>(); 
  for (edm::View<pat::Muon>::const_iterator itMuon = muons->begin(); itMuon != muons->end(); itMuon++) {

    pat::Muon aMuon(*itMuon);
    
    if ( useUserData_ ) {
      userDataHelper_.add( aMuon, ev, es );
    }
    std::vector<std::string>::const_iterator name = floatNames_.begin();
    for (std::vector<double>::const_iterator val = myFloats.begin(); val != myFloats.end(); ++val) {
      aMuon.addUserFloat(*name,*val);
      ++name;
    }

    PATMuons->push_back(aMuon);
  }

  // sort Muons in ET
  std::sort(PATMuons->begin(), PATMuons->end(), eTComparator_);
  // put genEvt object in Event
  std::auto_ptr<std::vector<pat::Muon> > myMuons(PATMuons);
  ev.put(myMuons);
}

// ------------ method called once each job just before starting event loop  ------------
void AddMuonUserData::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void AddMuonUserData::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void AddMuonUserData::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AddMuonUserData);
