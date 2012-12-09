// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <DataFormats/PatCandidates/interface/Jet.h>

// TFile Service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TH1.h"

class RA2ZInvTreeMaker : public edm::EDAnalyzer {
public:
  explicit RA2ZInvTreeMaker(const edm::ParameterSet&);
  ~RA2ZInvTreeMaker();
  
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
  edm::InputTag genLabel_, electronVetoSrc_, muonVetoSrc_, tauVetoSrc_, isoTrkVetoSrc_;
  edm::InputTag vertexSrc_, jetSrc_, htJetSrc_, bJetSrc_, htSrc_, mhtSrc_, metSrc_;
  std::string looseTopTaggerSrc_, nominalTopTaggerSrc_;
  bool          doPUReWeight_, runTopTagger_, storeExtraVetos_;
  edm::InputTag puWeightSrc_, eventWeightSrc_;

  //Output stuff
  edm::Service<TFileService> fs;
  TTree *reducedValues;

  int m_nBosons, m_nJetsPt30Eta50, m_nJetsPt30Eta24,
    m_nJetsPt50Eta25, m_nJetsPt50Eta25MInv,
    m_nJetsCSVM, m_nJetsCSVT, 
    m_Vertices,  m_event, m_run, m_lumi;
  double m_HT, m_HTMInv, m_MHT, m_MET,
    m_dPhiMHT1, m_dPhiMHT2, m_dPhiMHT3, m_dPhiMHT4, m_dPhiMHTMin, m_dPhiMHTMinBCSVM, m_dPhiMHTMinBCSVT,
    m_dPhiMET1, m_dPhiMET2, m_dPhiMET3, m_dPhiMET4, m_dPhiMETMin, m_dPhiMETMinBCSVM, m_dPhiMETMinBCSVT,
    m_EventWt, m_PUWt,
    m_boson1Pt, m_boson1Eta, m_boson1M, m_boson1MinDR, m_boson1DRJet1;
  double m_Jet1Pt, m_Jet1Eta,
    m_Jet2Pt, m_Jet2Eta,
    m_Jet3Pt, m_Jet3Eta,
    m_Jet4Pt, m_Jet4Eta;
  //m_boson2Pt, m_boson2Eta, m_boson2M;
  double m_JetCSVM1Pt, m_JetCSVM1Eta,
    m_JetCSVM2Pt, m_JetCSVM2Eta,
    m_JetCSVT1Pt, m_JetCSVT1Eta,
    m_JetCSVT2Pt, m_JetCSVT2Eta;

  ///Top tagger variables
  int m_loose_bestTopJetIdx, m_loose_pickedRemainingCombfatJetIdx;
  bool m_loose_remainPassCSVS, m_passLooseTopTagger,    m_passLooseTopJetIdx,
    m_passLooseTopMassCut,    m_passLooseCSVCut,    m_passLooseRemainingSystem,
    m_passLooseMT2Cuts;
  double m_loose_bestTopJetMass, m_loose_MTbJet, m_loose_MTbestTopJet, m_loose_MT2,
    m_loose_MTbestWJet, m_loose_MTbestbJet, m_loose_MTremainingTopJet, m_loose_linearCombMTbJetPlusMTbestTopJet;

  int m_nominal_bestTopJetIdx, m_nominal_pickedRemainingCombfatJetIdx;
  bool m_nominal_remainPassCSVS, m_passNominalTopTagger,    m_passNominalTopJetIdx,
    m_passNominalTopMassCut,    m_passNominalCSVCut,    m_passNominalRemainingSystem,
    m_passNominalMT2Cuts;
  double m_nominal_bestTopJetMass, m_nominal_MTbJet, m_nominal_MTbestTopJet, m_nominal_MT2,
    m_nominal_MTbestWJet, m_nominal_MTbestbJet, m_nominal_MTremainingTopJet, m_nominal_linearCombMTbJetPlusMTbestTopJet;

  //Lepton vetos
  bool m_passElVeto, m_passMuVeto, m_passTauVeto, m_passIsoTrkVeto;
};
