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
// $Id: GenStudyTree.h,v 1.5 2013/01/19 19:36:51 sturdy Exp $
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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

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
  //void studyRecoPhotons(photons,jets,gen,genjets,ht,genht,mht,genmht,met,genmet);
  //void studyRecoMuons();
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
  edm::InputTag vertexSrc_, recoPhotonSrc_, recoMuonSrc_;
  edm::InputTag recoJetSrc_,htJetSrc_, bJetSrc_,
    htSrc_, mhtSrc_, metSrc_, htNoBosonSrc_, mhtNoBosonSrc_, metNoBosonSrc_;
  edm::InputTag electronVetoSrc_, muonVetoSrc_, isoTrkVetoSrc_,
    ra2ElectronSrc_, ra2MuonSrc_;

  bool doPUReweight_, storeExtraVetos_;
  edm::InputTag puWeightSrc_;
  edm::InputTag eventWeightSrc_;
  
  edm::Service<TFileService> fs;
  TTree *reducedValues;

  bool studyAcc_, studyRecoIso_, removePhot_;
  bool studyPhotons_, studyMuons_;
  int nParticles_;

  double bosonMinPt_   ;
  double bosonEBMaxEta_;
  double bosonEEMinEta_;
  double bosonEEMaxEta_;

  int m_nJetsGenPt30Eta50, m_nJetsGenPt30Eta24, 
    m_bJetsGenPt30Eta24,   m_nJetsGenPt50Eta25,
    m_nJetsGenPt50Eta25MInv,
    m_Vertices,            m_genBosons;
  int m_nJetsPt30Eta24,m_nJetsPt50Eta24,m_nJetsPt70Eta24,
    m_nJetsCSVM, m_nJetsCSVT,
    m_nJetsPt30Eta50, m_nJetsPt50Eta25,
    m_nJetsPt50Eta25MInv,  m_nBosons;

  int m_event, m_run, m_lumi;
  bool m_genPassAcc,
    m_genMatchRecoID,      m_gen1MatchRecoID,      m_reco1MatchRecoID, 
    m_genMatchRecoTightID, m_gen1MatchRecoTightID, m_reco1MatchRecoTightID,
    m_genMatchRecoIDPixV,  m_gen1MatchRecoIDPixV,  m_reco1MatchRecoIDPixV,
    m_genMatchRecoIDCSEV,  m_gen1MatchRecoIDCSEV,  m_reco1MatchRecoIDCSEV,
    m_genMatchRecoIDIso,   m_gen1MatchRecoIDIso,   m_reco1MatchRecoIDIso;

  double m_genHT,  m_genMHT,       m_genHTMInv,  m_genMET,  m_genMETNoBoson,
    m_genBoson1Pt, m_genBoson1Eta, m_genBoson1M, m_genBoson1MinDR,
    m_genDPhiMHT1, m_genDPhiMHT2,  m_genDPhiMHT3,m_genDPhiMHT4,
    m_genDPhiMET1, m_genDPhiMET2,  m_genDPhiMET3,m_genDPhiMET4,
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

  bool m_boson1PassTight, m_boson1PassPixV, m_boson1PassCSEV, m_boson1PassIso;
  double m_daughter1Pt,   m_daughter1Eta,    m_daughter1M,    m_daughter1MinDR,   m_daughter1ID,   
         m_daughter2Pt,   m_daughter2Eta,    m_daughter2M,    m_daughter2MinDR,   m_daughter2ID,   
         m_combDaughterPt,m_combDaughterEta, m_combDaughterM, m_combDaughterMinDR;
  double m_Jet1Pt, m_Jet1Eta, 
    m_Jet2Pt,      m_Jet2Eta,
    m_Jet3Pt,      m_Jet3Eta,
    m_Jet4Pt,      m_Jet4Eta,
    m_boson2Pt,    m_boson2Eta, m_boson2M, m_boson2MinDR;
  
  bool m_passDirIsoElVeto, m_passDirIsoMuVeto, m_passIsoTrkVeto,
    m_passRA2ElVeto, m_passRA2MuVeto;
};

