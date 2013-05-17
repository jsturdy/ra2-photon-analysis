// -*- C++ -*-
//
// Package:    Muons
// Class:      Muons
// 
/**\class Muons SpecialMuonCollection.cc ZInvisibleBkgds/Photons/plugins/SpecialMuonCollection.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: SpecialMuonCollection.cc,v 1.3 2013/01/19 19:36:52 sturdy Exp $
//
//


#include "ZInvisibleBkgds/Photons/interface/SpecialMuonCollection.h"

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
SpecialMuonCollection::SpecialMuonCollection(const edm::ParameterSet& pset) :
  debug_         ( pset.getParameter< bool >( "debug" ) ),
  debugString_   ( pset.getParameter< std::string >( "debugString" ) ),
  candidateLabel_( pset.getParameter< edm::InputTag >( "candidateLabel" ) ),
  mZ_            ( pset.getParameter< double >( "mZ" ) )
{
  using namespace pat;
  // produces vector of muons
  produces<std::vector<pat::Muon> >();
}


SpecialMuonCollection::~SpecialMuonCollection()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void SpecialMuonCollection::produce(edm::Event& ev, const edm::EventSetup& es)
{
  using namespace pat;
  //get the muons
  edm::Handle<std::vector<reco::CompositeCandidate> > cands;
  ev.getByLabel(candidateLabel_,cands);

  std::vector<pat::Muon> * PATMuons = new std::vector<pat::Muon>(); 
  if (!(cands->size()>0)) {
    if (debug_)
      std::cout<<"Candidate collecition empty"<<std::endl;
    edm::LogWarning("SpecialMuonCollection")<<"Candidate collecition empty";
    // put genEvt object in Event
    std::auto_ptr<std::vector<pat::Muon> > myMuons(PATMuons);
    if (debug_)
      std::cout<<"Special muon collection has size "<<myMuons->size()<<std::endl;
    edm::LogWarning("SpecialMuonCollection")<<"Special muon collection has size "<<myMuons->size();
    ev.put(myMuons);
    return;
  }
  if (debug_)
    std::cout<<"Candidate collection has size "<<cands->size()<<std::endl;
  edm::LogWarning("SpecialMuonCollection")<<"Candidate collection has size "<<cands->size();
  //
  int candIndex = -1;
  double bestDMInv = 1000.;
  int ican(0);
  for (std::vector<reco::CompositeCandidate>::const_iterator itCand = cands->begin(); itCand != cands->end(); itCand++) {
    double dMInv = fabs(itCand->mass() - mZ_);
    if (debug_)
      std::cout<<"mInv = "<<itCand->mass()<<std::endl;
    edm::LogWarning("SpecialMuonCollection")<<"mInv = "<<itCand->mass();
    if (dMInv < bestDMInv) {
      bestDMInv = dMInv;
      candIndex = ican;
    }
    ++ican;
  }
  //reco::CompositeCandidate zCand((*cands)[candIndex]);
  //const reco::Candidate* muon1 = (*cands)[candIndex].daughter("muon1");
  //const reco::Candidate* muon2 = (*cands)[candIndex].daughter("muon2");
  pat::MuonRef master1(((*cands)[candIndex].daughter("muon1"))->masterClone().castTo<pat::MuonRef>());
  pat::MuonRef master2(((*cands)[candIndex].daughter("muon2"))->masterClone().castTo<pat::MuonRef>());
  //

  pat::Muon aMuon1(*master1);
  pat::Muon aMuon2(*master2);
  PATMuons->push_back(aMuon1);
  PATMuons->push_back(aMuon2);
  if (debug_){
    std::cout<<"muon1 "<<aMuon1.pt()<<", "<<aMuon1.eta()<<std::endl;
    std::cout<<"muon2 "<<aMuon2.pt()<<", "<<aMuon2.eta()<<std::endl;
  }
  edm::LogWarning("SpecialMuonCollection") << "muon1:"<<aMuon1.pt()<<", "<<aMuon1.eta();
  edm::LogWarning("SpecialMuonCollection") << "muon2:"<<aMuon2.pt()<<", "<<aMuon2.eta();
  // sort Muons in ET
  std::sort(PATMuons->begin(), PATMuons->end(), eTComparator_);
  // put genEvt object in Event
  std::auto_ptr<std::vector<pat::Muon> > myMuons(PATMuons);
  if (debug_)
    std::cout<<"Special muon collection has size "<<myMuons->size()<<std::endl;
  edm::LogWarning("SpecialMuonCollection")<<"Special muon collection has size "<<myMuons->size();
  ev.put(myMuons);
}

// ------------ method called once each job just before starting event loop  ------------
void SpecialMuonCollection::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void SpecialMuonCollection::endJob() {
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SpecialMuonCollection::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SpecialMuonCollection);
