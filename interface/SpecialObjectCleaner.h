#ifndef ZInvisibleBkgds_Photons_plugins_SpecialObjectCleaner_h
#define ZInvisibleBkgds_Photons_plugins_SpecialObjectCleaner_h
// -*- C++ -*-
//
// Package:    Objects
// Class:      Objects
// 
/**\class Objects SpecialObjectCleaner.h ZInvisibleBkgds/Photons/interface/SpecialObjectCleaner.h

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: SpecialObjectCleaner.h,v 1.1 2012/07/20 11:34:50 sturdy Exp $
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
    class SpecialObjectCleaner : public edm::EDProducer {
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
    
    explicit SpecialObjectCleaner(const edm::ParameterSet& pset);
    ~SpecialObjectCleaner() {};
    
    //static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    virtual void produce(edm::Event&, const edm::EventSetup&);
    
  private:
    //virtual void beginJob() ;
    //virtual void endJob() ;
    
    // ----------member data ---------------------------
  private:
    GreaterByEt<pat::Jet>            eTComparator_;
    pat::PATUserDataHelper<pat::Jet> userDataHelper_;
    
    bool debug_;
    std::string debugString_;
    edm::InputTag objectLabel_ ;
    edm::InputTag jetLabel_ ;
    double maxDR_;
    //std::string arbitration_;
    bool arbitration_;
  };
}
/////
template <class PATObjType>
zinvtools::SpecialObjectCleaner<PATObjType>::SpecialObjectCleaner(const edm::ParameterSet& pset) :
  debug_      ( pset.getParameter< bool >( "debug" ) ),
  debugString_( pset.getParameter< std::string >( "debugString" ) ),
  objectLabel_( pset.getParameter< edm::InputTag >( "objectLabel" ) ),
  jetLabel_   ( pset.getParameter< edm::InputTag >( "jetLabel" ) ),
  maxDR_      ( pset.getParameter< double >( "maxDR" ) ),
  //arbitration_( pset.getParameter< std::string >( "arbitration" ) )
  arbitration_( pset.getParameter< bool >( "arbitration" ) )
{
  using namespace pat;
  // produces vector of muons
  //if (arbitration_ == "bestDR") {
  //}
  //else if (arbitration_ == "highPt") {
  //}
  //else if (arbitration_ == "all") {
  //}
  produces<std::vector<pat::Jet> >();
}



//
// member functions
//

// ------------ method called to produce the data  ------------
template <class PATObjType>
void
zinvtools::SpecialObjectCleaner<PATObjType>::produce(edm::Event& ev, const edm::EventSetup& es)
{
  using namespace pat;
  //get the muons
  edm::Handle<std::vector<PATObjType> > objects;
  ev.getByLabel(objectLabel_,objects);

  edm::Handle<pat::JetCollection > jets;
  ev.getByLabel(jetLabel_,jets);
  
  if (!(objects->size()>0)) {
    std::cout<<"Object collecition empty"<<std::endl;
    return;
  }
  if (debug_) {
    std::cout<<"Object collection has size "<<objects->size()<<std::endl;
    std::cout<<"Jet collection has size "<<jets->size()<<std::endl;
  }

  //remove only closest DR object, or all objects within minDR?
  int candIndex = -1;
  double bestDR = 1000.;
  int ican(0);
  for (pat::JetCollection::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {
    if (!arbitration_) {
      double dR = reco::deltaR((*objects)[0].eta(),(*objects)[0].phi(),jet->eta(),jet->phi());
      
      if (dR < bestDR) {
	bestDR = dR;
	candIndex = ican;
      }
      ++ican;
    }
    if (arbitration_) {
      for (typename ObjectCollection::const_iterator itObj = objects->begin(); itObj != objects->end(); ++itObj) {
	double dR = reco::deltaR(itObj->eta(),itObj->phi(),jet->eta(),jet->phi());
	if (debug_)
	  std::cout<<"mInv = "<<itObj->mass()<<std::endl;
	if (dR < bestDR) {
	  bestDR = dR;
	  candIndex = ican;
	}
	++ican;
      }
    }
  }


  std::vector<pat::Jet> * PATCleanedJets = new std::vector<pat::Jet>(); 
  ican = 0;
  for (pat::JetCollection::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {
    
    if (ican == candIndex) {
      if (bestDR > maxDR_)
	PATCleanedJets->push_back(*jet);
    }
    else
      PATCleanedJets->push_back(*jet);
    
    ++ican;
  }
  
  // sort Objects in ET
  std::sort(PATCleanedJets->begin(), PATCleanedJets->end(), eTComparator_);
  std::auto_ptr<std::vector<pat::Jet> > myObjects(PATCleanedJets);
  // put object in Event
  if (debug_)
    std::cout<<"Cleaned jet collection has size "<<myObjects->size()<<std::endl;
  ev.put(myObjects);
}

#endif
