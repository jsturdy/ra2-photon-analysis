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
  edm::InputTag genLabel_;
  edm::InputTag vertexSrc_, jetSrc_, htJetSrc_, bJetSrc_, htSrc_, mhtSrc_, metSrc_;
  std::string topTaggerSrc_;
  bool          doPUReWeight_;
  edm::InputTag puWeightSrc_, eventWeightSrc_;

  //Output stuff
  edm::Service<TFileService> fs;
  TTree *reducedValues;

  int m_nBosons, m_nJetsPt30Eta50, m_bJetsPt30Eta24, 
    m_nJetsPt30Eta24, m_nJetsPt50Eta25, m_nJetsPt50Eta25MInv,
    m_Vertices;
  double m_HT, m_HTMInv, m_MHT, m_MET,
    m_dPhiMHT1, m_dPhiMHT2, m_dPhiMHT3, m_dPhiMHT4, m_dPhiMHTMin, m_dPhiMHTMinB,
    m_dPhiMET1, m_dPhiMET2, m_dPhiMET3, m_dPhiMET4, m_dPhiMETMin, m_dPhiMETMinB,
    m_EventWt, m_PUWt,
    m_boson1Pt, m_boson1Eta, m_boson1M, m_boson1MinDR;
  double m_Jet1Pt, m_Jet1Eta,
    m_Jet2Pt, m_Jet2Eta,
    m_Jet3Pt, m_Jet3Eta,
    m_Jet4Pt, m_Jet4Eta,
    m_boson2Pt, m_boson2Eta, m_boson2M;
  ///Top tagger variables
  double m_bestTopJetMass, m_TbJet, m_TbestTopJet, m_MT2;
};
