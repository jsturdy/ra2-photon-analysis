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

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/MET.h"
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

#include "ZInvisibleBkgds/Photons/interface/RA2ZInvTreeMakerFunctions.h"

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
  std::string debugString_;
  bool data_;
  double scale_;
  edm::InputTag photonSrc_, tightPhotonSrc_;
  edm::InputTag electronVetoSrc_, muonVetoSrc_, isoTrkVetoSrc_,
    ra2ElectronSrc_, ra2MuonSrc_;
  edm::InputTag vertexSrc_, jetSrc_, htJetSrc_,
    bJetSrc_, htSrc_, mhtSrc_, metSrc_,
    ra2HTSrc_, ra2MHTSrc_, ra2METSrc_;
  bool computeMET_, doPUReWeight_, storeExtraVetos_;
  edm::InputTag puWeightSrc_, eventWeightSrc_;

  //Trigger information
  bool getHLTfromConfig_; 
  bool checkedProcess_; 
  bool getL1Info_;
  edm::InputTag triggerResults_, hlTriggerResults_;
  std::string processName_;
  edm::TriggerNames triggerNames_;     // TriggerNames class
  HLTConfigProvider hltConfig;

  //Output stuff
  edm::Service<TFileService> fs;
  TTree *reducedValues;

  ///GEN level information
  bool runGenStudy_;
  edm::InputTag genSrc_, genJetSrc_, genMETSrc_, ra2JetSrc_;
  double maxDR_;
  int  m_genBosons, m_gen_nJetsPt30Eta50, m_gen_nJetsPt50Eta25,
    m_gen_nGenJetsPt30Eta50, m_gen_nGenJetsPt50Eta25;
  double m_genBoson1Pt, m_genBoson1Eta, m_genBoson1M, m_genBoson1MinDR, m_genBoson1DRJet1,
    m_genBoson2Pt, m_genBoson2Eta, m_genBoson2M, m_genBoson2MinDR, m_genBoson2DRJet1;
  bool m_genPassAcc,       m_genMatchRecoID,      m_gen1MatchRecoID,      m_reco1MatchRecoID, 
    m_genMatchRecoIDPixV,  m_gen1MatchRecoIDPixV,  m_reco1MatchRecoIDPixV,
    m_genMatchRecoIDCSEV,  m_gen1MatchRecoIDCSEV,  m_reco1MatchRecoIDCSEV,
    m_genMatchRecoIDIso,   m_gen1MatchRecoIDIso,   m_reco1MatchRecoIDIso;
  double m_gen_HT, m_gen_MHT, m_gen_GenHT, m_gen_GenMHT,
    m_gen_dPhiMHT1, m_gen_dPhiMHT2, m_gen_dPhiMHT3, m_gen_dPhiMHT4, m_gen_dPhiMHTMin,
    m_gen_genDPhiMHT1, m_gen_genDPhiMHT2, m_gen_genDPhiMHT3, m_gen_genDPhiMHT4, m_gen_genDPhiMHTMin,
    m_gen_genDPhiMET1, m_gen_genDPhiMET2, m_gen_genDPhiMET3, m_gen_genDPhiMET4, m_gen_genDPhiMETMin,
    m_gen_MET, m_gen_GenMET, 
    m_gen_dPhiMET1, m_gen_dPhiMET2, m_gen_dPhiMET3, m_gen_dPhiMET4, m_gen_dPhiMETMin;
  double m_genJet1Pt, m_genJet1Eta,
    m_genJet2Pt, m_genJet2Eta,
    m_genJet3Pt, m_genJet3Eta,
    m_genJet4Pt, m_genJet4Eta;

  ///RECO level information
  int m_nPhotonsID, m_nPhotonsTightIso, m_Photon1PDGID, m_Photon1Status, m_Photon2PDGID, m_Photon2Status,
    m_nJetsPt30Eta24, m_nJetsPt50Eta24, m_nJetsPt70Eta24,
    m_nJetsPt30Eta50, m_nJetsPt50Eta25, m_nJetsPt50Eta25MInv,
    m_nJetsCSVM, m_nJetsCSVT,
    m_Vertices,  m_event, m_run, m_lumi;
  double m_HT, m_HTMInv, m_MHT, m_ra2_HT, m_ra2_MHT,
    m_dPhiMHT1, m_dPhiMHT2, m_dPhiMHT3, m_dPhiMHT4, m_dPhiMHTMin,
    m_MET,m_ra2_MET,
    m_dPhiMET1, m_dPhiMET2, m_dPhiMET3, m_dPhiMET4, m_dPhiMETMin,
    m_EventWt, m_PUWt;
  double m_Photon1Pt, m_Photon1Eta, m_Photon1MinDR, m_Photon1DRJet1,
    m_Photon1Phi, m_Photon1SigmaIetaIeta, m_Photon1HadTowOverEm, 
    m_Photon1pfCH, m_Photon1pfNU, m_Photon1pfGA;
  bool  m_Photon1EConvVeto, m_Photon1PixelVeto, m_Photon1IsTightPFIso;
  double m_Photon2Pt, m_Photon2Eta, m_Photon2MinDR, m_Photon2DRJet1,
    m_Photon2Phi, m_Photon2SigmaIetaIeta, m_Photon2HadTowOverEm, 
    m_Photon2pfCH, m_Photon2pfNU, m_Photon2pfGA;
  bool  m_Photon2EConvVeto, m_Photon2PixelVeto, m_Photon2IsTightPFIso;

  double m_Jet1Pt, m_Jet1Eta,
    m_Jet2Pt, m_Jet2Eta,
    m_Jet3Pt, m_Jet3Eta,
    m_Jet4Pt, m_Jet4Eta;
  double m_JetCSVM1Pt, m_JetCSVM1Eta,
    m_JetCSVM2Pt, m_JetCSVM2Eta,
    m_JetCSVT1Pt, m_JetCSVT1Eta,
    m_JetCSVT2Pt, m_JetCSVT2Eta;
  bool m_Photon70PFMET100, m_Photon70PFHT400, m_Photon70PFNoPUHT400, m_Photon135, m_Photon150;

  //Lepton vetos
  bool m_passDirIsoElVeto, m_passDirIsoMuVeto, m_passIsoTrkVeto,
    m_passRA2ElVeto, m_passRA2MuVeto;
};
