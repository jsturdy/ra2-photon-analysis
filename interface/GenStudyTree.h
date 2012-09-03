// -*- C++ -*-
//
// Package:    Photons
// Class:      Photons
// 
/**\class Photons GenStudyTree.h ZInvisibleBkgds/Photons/interface/GenStudyTree.h

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: GenStudyTree.h,v 1.1 2012/08/30 09:44:42 sturdy Exp $
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
#include <TTree.h>

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

class GenStudyTree : public edm::EDProducer {
public:
  explicit GenStudyTree(const edm::ParameterSet&);
  ~GenStudyTree();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  // ----------member data ---------------------------
public:
  bool debug_;
  double scale_;
  std::string debugString_;
  
  edm::InputTag genSrc_;
  edm::InputTag genJetSrc_;
  edm::InputTag vertexSrc_, recoPhotonSrc_;
  edm::InputTag recoJetSrc_,htJetSrc_, bJetSrc_,
    htSrc_, mhtSrc_, htNoBosonSrc_, mhtNoBosonSrc_;

  bool doPUReweight_;
  edm::InputTag puWeightSrc_;
  
  edm::Service<TFileService> fs;
  TTree *reducedValues;

  bool   studyAcc_, studyRecoIso_, removePhot_;
  double bosonMinPt_   ;
  double bosonEBMaxEta_;
  double bosonEEMinEta_;
  double bosonEEMaxEta_;

  int m_nJetsGenPt30Eta50, m_nJetsGenPt30Eta24, 
    m_bJetsGenPt30Eta24,   m_nJetsGenPt50Eta25,
    m_Vertices,            m_genBosons;
  int m_nJetsPt30Eta50,    m_nJetsPt30Eta24,
    m_bJetsPt30Eta24,      m_nJetsPt50Eta25,
    m_nBosons;

  bool m_genPassAcc, m_genPassRecoIso;

  double m_genHT,  m_genMHT,       m_genHTNoBoson, m_genMHTNoBoson,
    m_genBoson1Pt, m_genBoson1Eta, m_genBoson1M,   
    m_genDPhi1,    m_genDPhi2,     m_genDPhi3,
    m_genDPhiNoBoson1,    m_genDPhiNoBoson2,     m_genDPhiNoBoson3,
    m_EventWt,     m_PUWt;
  double m_genJet1Pt, m_genJet1Eta, 
    m_genJet2Pt,      m_genJet2Eta,
    m_genJet3Pt,      m_genJet3Eta,
    m_genBoson2Pt,    m_genBoson2Eta, m_genBoson2M;

  double m_HT,  m_MHT,       m_HTNoBoson, m_MHTNoBoson,
    m_boson1Pt, m_boson1Eta, m_boson1M,   
    m_dPhi1,    m_dPhi2,     m_dPhi3;
  double m_Jet1Pt, m_Jet1Eta, 
    m_Jet2Pt,      m_Jet2Eta,
    m_Jet3Pt,      m_Jet3Eta,
    m_boson2Pt,    m_boson2Eta, m_boson2M;
  
};

