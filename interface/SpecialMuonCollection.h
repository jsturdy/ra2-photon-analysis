// -*- C++ -*-
//
// Package:    Muons
// Class:      Muons
// 
/**\class Muons SpecialMuonCollection.h ZInvisibleBkgds/Photons/interface/SpecialMuonCollection.h

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: SpecialMuonCollection.h,v 1.1 2012/07/20 11:34:50 sturdy Exp $
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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <DataFormats/Candidate/interface/Candidate.h>

//
// class declaration
//

class SpecialMuonCollection : public edm::EDProducer {
public:
  explicit SpecialMuonCollection(const edm::ParameterSet&);
  ~SpecialMuonCollection();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
private:
  GreaterByEt<pat::Muon> eTComparator_;

  bool debug_;
  std::string debugString_;
  edm::InputTag candidateLabel_ ;
  double mZ_;
  pat::PATUserDataHelper<pat::Muon>      userDataHelper_;
};

