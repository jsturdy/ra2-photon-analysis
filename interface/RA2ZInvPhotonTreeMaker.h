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
#include <DataFormats/PatCandidates/interface/Photon.h>
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

class RA2ZInvPhotonTreeMaker : public edm::EDAnalyzer {
public:
  explicit RA2ZInvPhotonTreeMaker(const edm::ParameterSet&);
  ~RA2ZInvPhotonTreeMaker();
  
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
  edm::InputTag photonSrc_;
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

  int m_nPhotonsIso, m_nJetsPt30Eta50, m_bJetsPt30Eta24, m_nJetsPt50Eta25, m_Vertices;
  double m_HT, m_MHT, m_dPhi1, m_dPhi2, m_dPhi3, m_EventWt, m_PUWt, m_Photon1Pt, m_Photon1Eta;
  double m_Jet1Pt, m_Jet1Eta, m_Jet2Pt, m_Jet2Eta, m_Jet3Pt, m_Jet3Eta, m_Photon2Pt, m_Photon2Eta;
  bool m_Photon70PFMET100, m_Photon70PFHT400, m_Photon70PFNoPUHT400, m_Photon135, m_Photon150;
};
