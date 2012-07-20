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
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>

// TFile Service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"

class RA2ZInvAnalyzer : public edm::EDAnalyzer {
public:
  explicit RA2ZInvAnalyzer(const edm::ParameterSet&);
  ~RA2ZInvAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  //80.0, 100.0, 120.0, 200,300,400,500 and >500  
  //0.0, 0.9, 1.442, 1.566, 2.1, 2.5, 5.0    
  static const int NPhotPtBins  = 8;
  static const int NPhotEtaBins = 7;

private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  void doGenAnalysis(edm::Handle<reco::GenParticleCollection>& genParticles, 
		     edm::Handle<edm::View<pat::Jet> >& jetsHT,
		     edm::Handle<edm::View<pat::Jet> >& jetsMHT,
		     std::vector<const pat::Photon*> IsoPhotons,
		     //edm::Handle<edm::View<pat::Photon> >& photonsIso,
		     double pu_event_wt);

  void BookHistograms();

  bool debug_;
  bool data_;
  edm::InputTag jetSrc_, jetNoPhotSrc_, jetHTSrc_;
  edm::InputTag genParticles_;
  bool          doOptimizePlots_;
  bool          doPUReWeight_;
  edm::InputTag puWeightSrc_;
  edm::InputTag pfRhoSrc_;

  int           ra2NJets_;
  double        ra2HT_;
  double        ra2MHT_;
  bool          ra2ApplyDphiCuts_;
  bool          doGenAnalysis_;

  edm::Service<TFileService> fs;
  TH1F *h_puWeight;
  TH1F *h_Vertices;
  TH1F *h_VerticesReWeighted;
  TH1F *h_pfRho[60];

  TH1F *h_PreRA2_Vertices;
  TH1F *h_PreRA2_VerticesReWeighted;
  TH1F *h_PreRA2_NJets_Pt30, *h_PreRA2_NJets_Pt50Eta25;
  TH1F *h_PreRA2_HT, *h_PreRA2_MHT, *h_PreRA2_MEff;
  TH1F *h_PreRA2_DPhiMHTJet1, *h_PreRA2_DPhiMHTJet2, *h_PreRA2_DPhiMHTJet3;

  TH1F *h_PreDPhi_Vertices;
  TH1F *h_PreDPhi_VerticesReWeighted;
  TH1F *h_PreDPhi_NJets_Pt30, *h_PreDPhi_NJets_Pt50Eta25;
  TH1F *h_PreDPhi_HT, *h_PreDPhi_MHT, *h_PreDPhi_MEff;
  TH1F *h_PreDPhi_DPhiMHTJet1, *h_PreDPhi_DPhiMHTJet2, *h_PreDPhi_DPhiMHTJet3;

  TH1F *h_RA2_Vertices;
  TH1F *h_RA2_VerticesReWeighted;
  TH1F *h_RA2_NJets_Pt30, *h_RA2_NJets_Pt50Eta25;
  TH1F *h_RA2_HT, *h_RA2_MHT, *h_RA2_MEff;
  TH1F *h_RA2_DPhiMHTJet1, *h_RA2_DPhiMHTJet2, *h_RA2_DPhiMHTJet3;

  TH1F *h_GenBoson_preRA2_dPhiMHTJet1, *h_GenBoson_preRA2_dPhiMHTJet2, *h_GenBoson_preRA2_dPhiMHTJet3;

  //MHT=200, 350, 500, 600,  800             (600-800 means>600)
  //HT =350, 500, 800, 1000, 1200, 1400, 1600 (1400-1600 means>1400)
  static const int NRA2MHTBins = 5;
  static const int NRA2HTBins  = 7;
  double RA2MHTVal[NRA2MHTBins] ;
  double RA2HTVal [NRA2HTBins];
  TH2F *h_RA2_HTMHT_Excl;
  TH2F *h_RA2_GenBosonPtEta_Excl[NRA2HTBins-1][NRA2MHTBins-1];


  // 
  static const int NMHTBins = 8;
  static const int NHTBins  = 8;
  double MHTVal[NMHTBins] ;
  double HTVal [NHTBins];

  TH2F *h_RA2_HTMHT_Incl;
  TH1F *h_RA2_NJets_HT350_MHT200;
  TH1F *h_RA2_NJets_HT800_MHT200;
  TH1F *h_RA2_NJets_HT800_MHT500;
  TH1F *h_RA2_NJets_HT500_MHT350;
  TH1F *h_RA2_NJets_HT500_MHT200;
  TH1F *h_RA2_NJets_HT1000_MHT600;
  TH1F *h_RA2_NJets_HT1200_MHT400;

  TH1F *h_RA2_HT_HT350_MHT200,  *h_RA2_MHT_HT350_MHT200,  *h_RA2_MEff_HT350_MHT200;
  TH1F *h_RA2_HT_HT800_MHT200,  *h_RA2_MHT_HT800_MHT200,  *h_RA2_MEff_HT800_MHT200;
  TH1F *h_RA2_HT_HT800_MHT500,  *h_RA2_MHT_HT800_MHT500,  *h_RA2_MEff_HT800_MHT500;
  TH1F *h_RA2_HT_HT500_MHT350,  *h_RA2_MHT_HT500_MHT350,  *h_RA2_MEff_HT500_MHT350;
  TH1F *h_RA2_HT_HT500_MHT200,  *h_RA2_MHT_HT500_MHT200,  *h_RA2_MEff_HT500_MHT200;
  TH1F *h_RA2_HT_HT1000_MHT600, *h_RA2_MHT_HT1000_MHT600, *h_RA2_MEff_HT1000_MHT600;
  TH1F *h_RA2_HT_HT1200_MHT400, *h_RA2_MHT_HT1200_MHT400, *h_RA2_MEff_HT1200_MHT400;

  TH2F *h_RA2_GenBosonPtEta_HT350_MHT200;
  TH2F *h_RA2_GenBosonPtEta_HT800_MHT200;
  TH2F *h_RA2_GenBosonPtEta_HT800_MHT500;
  TH2F *h_RA2_GenBosonPtEta_HT500_MHT350;
  TH2F *h_RA2_GenBosonPtEta_HT500_MHT200;
  TH2F *h_RA2_GenBosonPtEta_HT1000_MHT600;
  TH2F *h_RA2_GenBosonPtEta_HT1200_MHT400;
 

  //double GenBosonPtVal[NGenBosonPtBins];
  //double GenBosonEtaVal[NGenBosonEtaBins];
  TH2F *h_PreRA2_GenBosonPtEta;
  TH2F *h_RA2_GenBosonPtEta;

  // gen histograms
  TH1F *h_GenBoson_NBosons;
  TH1F *h_GenBoson_NJets_Pt30, *h_GenBoson_NJets_Pt50Eta25;
  TH1F *h_GenBoson_Pt, *h_GenBoson_Eta, *h_GenBoson_Phi;
  TH1F *h_GenBoson_HT;
  TH1F *h_GenBoson_MHT;
  TH2F *h_GenBoson_PtVsHT;
  TH2F *h_GenBoson_MHTVsHT;
  
  TH1F *h_GenBoson_RA2_Pt, *h_GenBoson_RA2_Eta, *h_GenBoson_RA2_Phi;
  TH1F *h_GenBoson_RA2_HT;
  TH1F *h_GenBoson_RA2_MHT;
  TH2F *h_GenBosonTotal_RA2_HTMHT_Excl;
  TH2F *h_GenBosonTotal_RA2_HTMHT_Incl;
  TH2F *h_GenBosonKineAcc_RA2_HTMHT_Excl;
  TH2F *h_GenBosonKineAcc_RA2_HTMHT_Incl;
  TH2F *h_GenBosonMatched_RA2_HTMHT_Excl;
  TH2F *h_GenBosonMatched_RA2_HTMHT_Incl;
};
