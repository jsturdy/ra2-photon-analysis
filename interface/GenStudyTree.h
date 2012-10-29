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
// $Id: GenStudyTree.h,v 1.1 2012/09/03 10:30:38 sturdy Exp $
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
  edm::InputTag genJetSrc_, genMETSrc_;
  edm::InputTag vertexSrc_, recoPhotonSrc_;
  edm::InputTag recoJetSrc_,htJetSrc_, bJetSrc_,
    htSrc_, mhtSrc_, metSrc_, htNoBosonSrc_, mhtNoBosonSrc_;

  bool doPUReweight_;
  edm::InputTag puWeightSrc_;
  edm::InputTag eventWeightSrc_;
  
  edm::Service<TFileService> fs;
  TTree *reducedValues;

  bool   studyAcc_, studyRecoIso_, removePhot_;
  double bosonMinPt_   ;
  double bosonEBMaxEta_;
  double bosonEEMinEta_;
  double bosonEEMaxEta_;

  int m_nJetsGenPt30Eta50, m_nJetsGenPt30Eta24, 
    m_bJetsGenPt30Eta24,   m_nJetsGenPt50Eta25,
    m_nJetsGenPt50Eta25MInv,
    m_Vertices,            m_genBosons;
  int m_nJetsPt30Eta50,    m_nJetsPt30Eta24,
    m_bJetsPt30Eta24,      m_nJetsPt50Eta25,
    m_nJetsPt50Eta25MInv,  m_nBosons;

  bool m_genPassAcc, m_genPassRecoIso, m_gen1PassRecoIso;

  double m_genHT,  m_genMHT,       m_genHTMInv,  m_genMET,  m_genMETNoBoson,
    m_genBoson1Pt, m_genBoson1Eta, m_genBoson1M, m_genBoson1MinDR,
    m_genDPhiMHT1,    m_genDPhiMHT2,     m_genDPhiMHT3,   m_genDPhiMHT4,
    m_genDPhiMET1,    m_genDPhiMET2,     m_genDPhiMET3,   m_genDPhiMET4,
    m_genDPhiMETNoBoson1,    m_genDPhiMETNoBoson2,     m_genDPhiMETNoBoson3,     m_genDPhiMETNoBoson4,
    m_EventWt,     m_PUWt;
  double m_genJet1Pt, m_genJet1Eta, 
    m_genJet2Pt,      m_genJet2Eta,
    m_genJet3Pt,      m_genJet3Eta,
    m_genJet4Pt,      m_genJet4Eta,
    m_genBoson2Pt,    m_genBoson2Eta, m_genBoson2M;

  double m_HT,  m_MHT,       m_HTMInv,    m_MET,    m_METNoBoson,
    m_boson1Pt, m_boson1Eta, m_boson1M, m_boson1MinDR,
    m_dPhiMHT1,    m_dPhiMHT2,     m_dPhiMHT3,     m_dPhiMHT4,
    m_dPhiMET1,    m_dPhiMET2,     m_dPhiMET3,     m_dPhiMET4,
    m_dPhiMETNoBoson1,    m_dPhiMETNoBoson2,     m_dPhiMETNoBoson3,     m_dPhiMETNoBoson4;

  double m_daughter1Pt, m_daughter1Eta,    m_daughter1M,   m_daughter1ID,   
    m_daughter2Pt,      m_daughter2Eta,    m_daughter2M,   m_daughter2ID,   
    m_combDaughterPt,   m_combDaughterEta, m_combDaughterM;
  double m_Jet1Pt, m_Jet1Eta, 
    m_Jet2Pt,      m_Jet2Eta,
    m_Jet3Pt,      m_Jet3Eta,
    m_Jet4Pt,      m_Jet4Eta,
    m_boson2Pt,    m_boson2Eta, m_boson2M;
  
};

