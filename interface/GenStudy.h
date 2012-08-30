// -*- C++ -*-
//
// Package:    Photons
// Class:      Photons
// 
/**\class Photons GenStudy.h ZInvisibleBkgds/Photons/interface/GenStudy.h

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: GenStudy.h,v 1.1 2012/07/09 14:31:26 sturdy Exp $
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
#include <TH3.h>

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

class GenStudy : public edm::EDProducer {
public:
  explicit GenStudy(const edm::ParameterSet&);
  ~GenStudy();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
public:
  bool debug_;
  std::string debugString_;
  
  edm::InputTag genLabel_;
  int pdgId_;
  int genStatus_;
  int mom_pdgId_;
  int mom_genStatus_;
  edm::InputTag genJetLabel_;
  edm::InputTag recoJetLabel_;
  edm::InputTag recoPhotonLabel_;

  bool doPUReweight_;
  edm::InputTag puWeightLabel_;
  
  std::vector<double> nJetBins_ ;
  std::vector<double> htBins_ ;
  std::vector<double> mhtBins_;
  double minHT_, minMHT_ ;

  std::vector<double> bosonPtBins_ ;
  std::vector<double> bosonEtaBins_;

  bool studyAcc_, studyRecoIso_, removePhot_;
  double bosonMinPt_        ;
  double bosonEBMaxEta_       ;
  double bosonEEMinEta_       ;
  double bosonEEMaxEta_       ;

  std::map<std::string, TH2D*> histos2D_;
  std::map<std::string, TH3D*> histos3D_;
//  std::map<std::string, TH1D*> histos1DNJets_[4];
//  std::map<std::string, TH3D*> histos3DNJets_[4];
};

