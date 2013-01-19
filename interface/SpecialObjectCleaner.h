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
// $Id: SpecialObjectCleaner.h,v 1.2 2012/09/03 10:33:29 sturdy Exp $
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
    /// references to Object collection
    typedef edm::RefProd<ObjectCollection> ObjectRefProd;
    /// vector of references to Object objects all in the same collection
    typedef edm::RefVector<ObjectCollection> ObjectRefVector;
    /// iterator over a vector of references to Object objects all in the same collection
    typedef typename ObjectRefVector::iterator object_iterator;
    
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
    if (debug_) 
      std::cout<<"Object collecition empty"<<std::endl;
    std::vector<pat::Jet> * PATCleanedJets = new std::vector<pat::Jet>(); 
    for (pat::JetCollection::const_iterator jet = jets->begin(); jet != jets->end(); ++jet)
      PATCleanedJets->push_back(*jet);
    
    std::sort(PATCleanedJets->begin(), PATCleanedJets->end(), eTComparator_);
    std::auto_ptr<std::vector<pat::Jet> > myObjects(PATCleanedJets);
    // put object in Event
    if (debug_)
      std::cout<<"Cleaned jet collection has size "<<myObjects->size()<<std::endl;
    ev.put(myObjects);
    return;
  }
  if (debug_) {
    std::cout<<"Object collection has size "<<objects->size()<<std::endl;
    printf("Object      pt      eta      phi        M\n");
    for (unsigned int i = 0;i < objects->size();++i) 
      printf("%d  %2.4f    %2.2f    %2.2f    %2.2f\n",
	     i,(*objects)[i].pt(),(*objects)[i].eta(),(*objects)[i].phi(),(*objects)[i].mass());
    std::cout<<"Jet collection has size "<<jets->size()<<std::endl;
    printf("jet      pt      eta      phi        chf   nhf   cemf   nemf  phef   elef   muef  hfhf  hfemf\n");
    for (unsigned int i = 0;i < jets->size();++i) 
      if (i < 10)
	printf("%d  %2.4f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f    %2.2f\n",
	       i,(*jets)[i].pt(),(*jets)[i].eta(),(*jets)[i].phi(),
	       (*jets)[i].chargedHadronEnergyFraction(),
	       (*jets)[i].neutralHadronEnergyFraction(),
	       (*jets)[i].chargedEmEnergyFraction(),
	       (*jets)[i].neutralEmEnergyFraction(),
	       (*jets)[i].photonEnergyFraction(),
	       (*jets)[i].electronEnergyFraction(),
	       (*jets)[i].muonEnergyFraction(),
	       (*jets)[i].HFHadronEnergyFraction(),
	       (*jets)[i].HFEMEnergyFraction()
	       );
  }

  //remove only closest DR object, or all objects within minDR?
  int candIndex = -1;
  int objIndex  = -1;
  double bestDR = 1000.;
  int ican(0);

  for (pat::JetCollection::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {
    if (!arbitration_) {
      double dR = reco::deltaR((*objects)[0].eta(),(*objects)[0].phi(),jet->eta(),jet->phi());
      if (dR < bestDR) {
	bestDR = dR;
	candIndex = ican;
	objIndex  = 0;
      }
      //++ican;
    }
    if (arbitration_) {
      int iobj(0);
      for (typename ObjectCollection::const_iterator itObj = objects->begin(); itObj != objects->end(); ++itObj) {
	double dR = reco::deltaR(itObj->eta(),itObj->phi(),jet->eta(),jet->phi());
	if (debug_)
	  std::cout<<"mInv = "<<itObj->mass()<<std::endl;
	if (dR < bestDR) {
	  bestDR = dR;
	  candIndex = ican;
	  objIndex = iobj;
	}
	++iobj;
      }
    }
    ++ican;
  }


  std::vector<pat::Jet> * PATCleanedJets = new std::vector<pat::Jet>(); 
  ican = 0;
  for (pat::JetCollection::const_iterator jet = jets->begin(); jet != jets->end(); ++jet) {
    
    if (ican == candIndex) {
      if (bestDR > maxDR_)
	PATCleanedJets->push_back(*jet);
      else 
	if (debug_)
	  std::cout<<"Removing jet "<<candIndex
		   <<" - pt("<<jet->pt()
		   <<"), eta("<<jet->eta()
		   <<"), phi("<<jet->phi()
		   <<")"<<std::endl<<"matched DR("<<bestDR<<") to object "<<objIndex
		   <<" - pt("<<(*objects)[objIndex].pt()
		   <<"), eta("<<(*objects)[objIndex].eta()
		   <<"), phi("<<(*objects)[objIndex].phi()
		   <<")"
		   <<std::endl;
    }
    else
      PATCleanedJets->push_back(*jet);
    
    ++ican;
  }
  
  // sort Objects in ET
  std::sort(PATCleanedJets->begin(), PATCleanedJets->end(), eTComparator_);
  std::auto_ptr<std::vector<pat::Jet> > myObjects(PATCleanedJets);
  // put object in Event
  if (debug_) {
    std::cout<<"Object collection has size "<<objects->size()<<std::endl;
    std::cout<<"Cleaned jet collection has size "<<myObjects->size()<<std::endl;
    std::cout<<"removed "<<jets->size()-myObjects->size()<<" objects"<<std::endl;
  }
  ev.put(myObjects);
}

#endif
