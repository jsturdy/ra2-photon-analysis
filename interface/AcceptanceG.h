// -*- C++ -*-
//
// Package:    Photons
// Class:      Photons
// 
/**\class Photons AcceptanceG.h ZInvisibleBkgds/Photons/interface/AcceptanceG.h

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: AcceptanceG.h,v 1.1 2012/05/16 20:25:39 sturdy Exp $
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
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//Used data formats
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

//For MC truth information
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

//
// class declaration
//

class AcceptanceG : public edm::EDProducer {
public:
  explicit AcceptanceG(const edm::ParameterSet&);
  ~AcceptanceG();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
private:
  bool debug_;
  std::string debugString_;
  
  edm::InputTag genLabel_;
  int genStatus_;

  bool doPUReweight_;
  edm::InputTag puWeightLabel_;

  edm::InputTag jetHTLabel_;
  edm::InputTag jetMHTLabel_;
  edm::InputTag htLabel_;
  edm::InputTag mhtLabel_;
  std::vector<double> htBins_ ;
  std::vector<double> mhtBins_;

  edm::InputTag cleanJetHTLabel_;
  edm::InputTag cleanJetMHTLabel_;
  edm::InputTag cleanHTLabel_;
  edm::InputTag cleanMHTLabel_;
  
  edm::InputTag photonLabel_ ;
  std::vector<double> photonPtBins_ ;
  std::vector<double> photonEtaBins_;
  double photonMinPt_        ;
  double photonEBMaxEta_       ;
  double photonEEMinEta_       ;
  double photonEEMaxEta_       ;

  //double photonPtBinsTmp_[7]  ;// = photonPtBins_.data();
  //double photonEtaBinsTmp_[11];// = photonEtaBins_.data();
  std::map<std::string, TH1D*> histos1D_;
  std::map<std::string, TH2D*> histos2D_;
};

