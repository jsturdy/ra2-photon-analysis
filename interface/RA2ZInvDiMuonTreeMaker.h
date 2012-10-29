// system include files
#include <memory>
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/Candidate/interface/Candidate.h>
#include <DataFormats/PatCandidates/interface/Jet.h>

// Trigger includes
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "FWCore/Common/interface/TriggerNames.h"

// TFile Service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TH1.h"

class RA2ZInvDiMuonTreeMaker : public edm::EDAnalyzer {
public:
  explicit RA2ZInvDiMuonTreeMaker(const edm::ParameterSet&);
  ~RA2ZInvDiMuonTreeMaker();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  void BookTree();

  bool debug_;
  bool data_;
  double scale_;
  edm::InputTag muonSrc_;
  edm::InputTag vertexSrc_, jetSrc_, htJetSrc_, bJetSrc_, htSrc_, mhtSrc_, metSrc_;
  std::string topTaggerSrc_;
  bool          doPUReWeight_;
  edm::InputTag puWeightSrc_, eventWeightSrc_;

  //Trigger information
  bool getHLTfromConfig_; 
  bool checkedProcess_; 
  bool getL1Info_;
  edm::InputTag hlTriggerResults_; 
  std::string processName_;
  edm::TriggerNames triggerNames_;     // TriggerNames class
  HLTConfigProvider hltConfig;

  //Output stuff
  edm::Service<TFileService> fs;
  TTree *reducedValues;

  int m_nMuonsIso, m_nJetsPt30Eta50, m_nJetsPt30Eta24,
    m_bJetsPt30Eta24, m_nJetsPt50Eta25, m_nJetsPt50Eta25MInv,
    m_Vertices;
  double m_HT, m_HTMInv, m_MHT, m_MET, m_modMET,
    m_dPhiMHT1, m_dPhiMHT2, m_dPhiMHT3, m_dPhiMHT4, m_dPhiMHTMin, m_dPhiMHTMinB,
    m_dPhiMET1, m_dPhiMET2, m_dPhiMET3, m_dPhiMET4, m_dPhiMETMin, m_dPhiMETMinB,
    m_dPhiModMET1, m_dPhiModMET2, m_dPhiModMET3, m_dPhiModMET4, m_dPhiModMETMin, m_dPhiModMETMinB,
    m_EventWt, m_PUWt, 
    m_Muon1Pt, m_Muon1Eta, m_Muon2Pt, m_Muon2Eta,
    m_DiMuonInvM, m_DiMuonPt;
  double m_Jet1Pt, m_Jet1Eta,
    m_Jet2Pt, m_Jet2Eta,
    m_Jet3Pt, m_Jet3Eta,
    m_Jet4Pt, m_Jet4Eta;

  bool m_Mu13_Mu8, m_Mu17_Mu8, m_DoubleMu5_IsoMu5;
  
  ///Top tagger variables
  double m_bestTopJetMass, m_TbJet, m_TbestTopJet, m_MT2;
};
