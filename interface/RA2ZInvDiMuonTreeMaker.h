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
  edm::InputTag vertexSrc_, jetSrc_, htJetSrc_, bJetSrc_, htSrc_, mhtSrc_;
  bool          doPUReWeight_;
  edm::InputTag puWeightSrc_;

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

  int m_nMuonsIso, m_nJetsPt30Eta50, m_bJetsPt30Eta24, m_nJetsPt50Eta25, m_Vertices;
  double m_HT, m_MHT, m_dPhi1, m_dPhi2, m_dPhi3, m_EventWt, m_PUWt, 
    m_Muon1Pt, m_Muon1Eta, m_Muon2Pt, m_Muon2Eta, m_DiMuonInvM, m_DiMuonPt;
  double m_Jet1Pt, m_Jet1Eta, m_Jet2Pt, m_Jet2Eta, m_Jet3Pt, m_Jet3Eta;

  bool m_DoubleMu8_Mass8_PFHT175, m_DoubleMu8_Mass8_PFHT225, 
    m_DoubleMu8_Mass8_PFNoPUHT175, m_DoubleMu8_Mass8_PFNoPUHT225, 
    m_DoubleRelIso1p0Mu5_Mass8_PFHT175, m_DoubleRelIso1p0Mu5_Mass8_PFHT225, 
    m_DoubleRelIso1p0Mu5_Mass8_PFNoPUHT175, m_DoubleRelIso1p0Mu5_Mass8_PFNoPUHT225,
    m_Mu13_Mu8, m_Mu17_Mu8, m_DoubleMu5_IsoMu5;
};
