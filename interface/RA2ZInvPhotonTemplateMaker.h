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

class RA2ZInvPhotonTemplateMaker : public edm::EDAnalyzer {
public:
  explicit RA2ZInvPhotonTemplateMaker(const edm::ParameterSet&);
  ~RA2ZInvPhotonTemplateMaker();
  
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
  edm::InputTag vertexSrc_, jetSrc_, htJetSrc_, htSrc_, mhtSrc_, metSrc_;
  bool          doPUReWeight_, runTopTagger_;
  edm::InputTag puWeightSrc_, eventWeightSrc_;

  //Output stuff
  edm::Service<TFileService> fs;
  TTree *reducedValues;

  int m_nPhotonsID, m_nPhotonsTightIso,
    m_Photon1PDGID,
    m_nJetsPt30Eta50, m_nJetsPt50Eta25,
    m_Vertices,  m_event, m_run, m_lumi;
  double m_HT, m_HTMInv, m_MHT, m_MET,
    m_dPhiMHT1, m_dPhiMHT2, m_dPhiMHT3, m_dPhiMHT4, m_dPhiMHTMin,
    m_dPhiMET1, m_dPhiMET2, m_dPhiMET3, m_dPhiMET4, m_dPhiMETMin,
    m_EventWt, m_PUWt,
    m_Photon1Pt, m_Photon1Eta, m_Photon1MinDR, m_Photon1DRJet1,
    m_Photon1Phi, 
    m_Photon1SigmaIetaIeta, m_Photon1HadTowOverEm, 
    m_Photon1pfCH, m_Photon1pfNU, m_Photon1pfGA;
  bool  m_Photon1EConvVeto, m_Photon1PixelVeto, m_Photon1IsTightID, m_Photon1IsTightPFIso,
    m_Photon1PassPFCh, m_Photon1PassPFNu, m_Photon1PassPFGa;

  bool m_Photon70PFMET100, m_Photon70PFHT400, m_Photon70PFNoPUHT400, m_Photon135, m_Photon150;
};
