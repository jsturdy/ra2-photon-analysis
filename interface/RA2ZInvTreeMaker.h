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
  edm::InputTag vertexSrc_, jetSrc_, bJetSrc_, jetHTSrc_, genLabel_;
  bool          doPUReWeight_;
  edm::InputTag puWeightSrc_;

  //int           ra2NJets_;
  //double        ra2HT_;
  //double        ra2MHT_;
  //bool          ra2ApplyDphiCuts_;

  //Output stuff
  edm::Service<TFileService> fs;
  TTree *reducedValues;

  int m_nJetsPt30Eta50, m_bJetsPt30Eta24, m_nJetsPt50Eta25, m_Vertices, m_nBosons;
  double m_HT, m_MHT, m_bosonPt, m_bosonEta, m_bosonM, m_dPhi1, m_dPhi2, m_dPhi3, m_EventWt, m_PUWt;
};