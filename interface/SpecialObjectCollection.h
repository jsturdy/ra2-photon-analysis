#ifndef ZInvisibleBkgds_Photons_plugins_SpecialObjectCollection_h
#define ZInvisibleBkgds_Photons_plugins_SpecialObjectCollection_h
// -*- C++ -*-
//
// Package:    Objects
// Class:      Objects
// 
/**\class Objects SpecialObjectCollection.h ZInvisibleBkgds/Photons/interface/SpecialObjectCollection.h

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: SpecialObjectCollection.h,v 1.1 2012/08/20 13:00:14 sturdy Exp $
//
//


// system include files
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <set>
#include <math.h>
#include <utility>

//ROOT includes
#include <TH1.h>
#include <TH2.h>

// Framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Run.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/UserData.h"
#include "PhysicsTools/PatAlgos/interface/PATUserDataHelper.h"
#include "CommonTools/Utils/interface/EtComparator.h"


//Used data formats
//#include "DataFormats/PatCandidates/interface/Object.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include <DataFormats/Candidate/interface/Candidate.h>

//
// class declaration
//
namespace zinvtools {

  template<class PATObjType>
    class SpecialObjectCollection : public edm::EDProducer {
  public:
    /// collection of Object objects
    typedef std::vector<PATObjType> ObjectCollection;
    /// presistent reference to a Object
    typedef edm::Ref<ObjectCollection> ObjectRef;
    ///// references to Object collection
    //typedef edm::RefProd<ObjectCollection> ObjectRefProd;
    ///// vector of references to Object objects all in the same collection
    //typedef edm::RefVector<ObjectCollection> ObjectRefVector;
    ///// iterator over a vector of references to Object objects all in the same collection
    //typedef ObjectRefVector::iterator object_iterator;
    
    explicit SpecialObjectCollection(const edm::ParameterSet& pset);
    ~SpecialObjectCollection() {};
    
    //static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
  private:
    //virtual void beginJob() ;
    //virtual void endJob() ;
    
    // ----------member data ---------------------------
  private:
    typedef GreaterByEt<PATObjType> ObjectPtComparator;
    typedef pat::PATUserDataHelper<PATObjType> ObjectUserDataHelper;

    ObjectPtComparator eTComparator_;
    ObjectUserDataHelper userDataHelper_;
    
    bool debug_;
    std::string debugString_;
    edm::InputTag candidateLabel_ ;
    double mZ_;
  };
}
/////
template <class PATObjType> 
zinvtools::SpecialObjectCollection<PATObjType>::SpecialObjectCollection(const edm::ParameterSet& pset) :
  debug_         ( pset.getParameter< bool >( "debug" ) ),
  debugString_   ( pset.getParameter< std::string >( "debugString" ) ),
  candidateLabel_( pset.getParameter< edm::InputTag >( "candidateLabel" ) ),
  mZ_            ( pset.getParameter< double >( "mZ" ) )
{
  using namespace pat;
  // produces vector of muons
  produces<std::vector<PATObjType> >();
}



//
// member functions
//

// ------------ method called to produce the data  ------------
template <class PATObjType>
void
zinvtools::SpecialObjectCollection<PATObjType>::produce(edm::Event& ev, const edm::EventSetup& es)
{
  using namespace pat;
  //get the muons
  //edm::Handle<std::vector<reco::CompositeCandidate> > cands;
  edm::Handle<reco::CompositeCandidateCollection > cands;
  ev.getByLabel(candidateLabel_,cands);
  
  if (!(cands->size()>0)) {
    std::cout<<"Candidate collecition empty"<<std::endl;
    return;
  }
  if (debug_)
    std::cout<<"Candidate collection has size "<<cands->size()<<std::endl;
  //
  int candIndex = -1;
  double bestDMInv = 1000.;
  int ican(0);
  for (std::vector<reco::CompositeCandidate>::const_iterator itCand = cands->begin(); itCand != cands->end(); itCand++) {
    double dMInv = fabs(itCand->mass() - mZ_);
    if (debug_)
      std::cout<<"mInv = "<<itCand->mass()<<std::endl;
    if (dMInv < bestDMInv) {
      bestDMInv = dMInv;
      candIndex = ican;
    }
    ++ican;
  }
  //std::vector<PATObjType> * PATObjects = new std::vector<PATObjType>(); 
  ObjectCollection * PATObjects = new std::vector<PATObjType>(); 
  //const_iterator_imp_specific 
  for (unsigned int dau = 0; dau < (*cands)[candIndex].numberOfDaughters(); ++dau){
    ObjectRef master(((*cands)[candIndex].daughter(dau))->masterClone().castTo<ObjectRef>());
    PATObjType aObject(*master);
    PATObjects->push_back(aObject);
    if (debug_)
      std::cout<<"object "<<aObject.pt()<<", "<<aObject.eta()<<std::endl;
  }
  //
  
  // sort Objects in ET
  std::sort(PATObjects->begin(), PATObjects->end(), eTComparator_);
  // put object in Event
  std::auto_ptr<std::vector<PATObjType> > myObjects(PATObjects);
  if (debug_)
    std::cout<<"Special object collection has size "<<myObjects->size()<<std::endl;
  ev.put(myObjects);
}

//// ------------ method called once each job just before starting event loop  ------------
//template <class PATObjType>
//void 
//zinvtools::SpecialObjectCollection<PATObjType>::beginJob()
//{
//}
//
//// ------------ method called once each job just after ending the event loop  ------------
//template <class PATObjType>
//void
//zinvtools::SpecialObjectCollection<PATObjType>::endJob() 
//{
//}
//
//// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//template <class PATObjType>
//void
//zinvtools::SpecialObjectCollection<PATObjType>::fillDescriptions(edm::ConfigurationDescriptions& descriptions) 
//{
//  //The following says we do not know what parameters are allowed so do no validation
//  // Please change this to state exactly what you do use, even if it is no parameters
//  edm::ParameterSetDescription desc;
//  desc.setUnknown();
//  descriptions.addDefault(desc);
//}

#endif
