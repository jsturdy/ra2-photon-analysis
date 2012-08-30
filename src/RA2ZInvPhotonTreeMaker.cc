// -*- C++ -*-
//
// Package:    RA2ZInvPhotonTreeMaker
// Class:      RA2ZInvPhotonTreeMaker
// 
/**\class RA2ZInvPhotonTreeMaker RA2ZInvPhotonTreeMaker.cc SusyAnalysis/RA2ZInvPhotonTreeMaker/src/RA2ZInvPhotonTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seema Sharma
//         Created:  Mon Jun 20 12:58:08 CDT 2011
// $Id: RA2ZInvPhotonTreeMaker.cc,v 1.1 2012/07/20 11:35:34 sturdy Exp $
//
//


// system include files
#include <memory>
#include <iomanip>
#include <cmath>

#include "ZInvisibleBkgds/Photons/interface/RA2ZInvPhotonTreeMaker.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include <DataFormats/METReco/interface/MET.h>


RA2ZInvPhotonTreeMaker::RA2ZInvPhotonTreeMaker(const edm::ParameterSet& pset) {

  // read parameters from config file
  debug_          = pset.getParameter<bool>("Debug");
  data_           = pset.getParameter<bool>("Data");
  scale_          = pset.getParameter<double>("ScaleFactor");
  photonSrc_      = pset.getParameter<edm::InputTag>("PhotonSrc");
  vertexSrc_      = pset.getParameter<edm::InputTag>("VertexSrc");
  jetSrc_         = pset.getParameter<edm::InputTag>("JetSrc");
  bJetSrc_        = pset.getParameter<edm::InputTag>("bJetSrc");
  jetHTSrc_       = pset.getParameter<edm::InputTag>("JetHTSource");
  doPUReWeight_   = pset.getParameter<bool>("DoPUReweight");
  puWeightSrc_    = pset.getParameter<edm::InputTag>("PUWeightSource");
  
  //getHLTfromConfig_   = pset.getUntrackedParameter<bool>("getHLTfromConfig",false);
  //if (getHLTfromConfig_)
  //  hlTriggerResults_ = pset.getUntrackedParameter<edm::InputTag>("hlTriggerResults");

  //key to help getting the hlt process from event provenance
  checkedProcess_ = false;
  processName_    = "";

}

RA2ZInvPhotonTreeMaker::~RA2ZInvPhotonTreeMaker() {
  delete reducedValues;
  reducedValues = 0; 
}

void RA2ZInvPhotonTreeMaker::analyze(const edm::Event& ev, const edm::EventSetup& es) {

  using namespace edm;

  // get event-id
  unsigned int event = (ev.id()).event();
  unsigned int run   = (ev.id()).run();
  unsigned int lumi  =  ev.luminosityBlock();

  // reject event if there an isolated electron or muon present
  edm::Handle< std::vector<pat::Electron> > patElectrons;
  ev.getByLabel("patElectronsIDIso", patElectrons);

  edm::Handle< std::vector<pat::Muon> > patMuons;
  ev.getByLabel("patMuonsIDIso", patMuons);
  
  if (patElectrons->size() != 0 || patMuons->size() != 0) { 
    std::cout << "Isolated Lepton found : Event Rejected : ( run, event, lumi ) " 
    << run << " " << event << " " << lumi << std::endl;
    return;
  }
  
  // get photons 
  edm::Handle< std::vector<pat::Photon> > patPhotons;
  ev.getByLabel(photonSrc_, patPhotons); 

  edm::Handle< reco::VertexCollection > vertices;
  ev.getByLabel(vertexSrc_, vertices); 
  
  // get jetcollection
  edm::Handle< std::vector<pat::Jet> > recoJets;
  ev.getByLabel(jetSrc_, recoJets); 
  
  edm::Handle<edm::View<pat::Jet> > jets;
  ev.getByLabel(jetSrc_, jets);

  // get jetcollection
  edm::Handle< std::vector<pat::Jet> > recobJets;
  ev.getByLabel(bJetSrc_, recobJets); 

  edm::Handle<edm::View<pat::Jet> > bJets;
  ev.getByLabel(bJetSrc_, bJets);

  edm::Handle<edm::View<pat::Jet> > jetsHT;
  ev.getByLabel(jetHTSrc_, jetsHT);

  // if MC, do PU reweighting
  double pu_event_wt = 1.0;
  edm::Handle<double> puweight;
  if( doPUReWeight_ ) {
    ev.getByLabel(puWeightSrc_, puweight);
    pu_event_wt = *puweight;
  }

  m_EventWt = scale_;
  m_PUWt    = pu_event_wt;
  m_Vertices = vertices->size();
  std::vector<const pat::Photon*> IsoPhotons;
  std::vector<const pat::Photon*> RA2Photons;
  double photon_ptcut  = 30.0;
  
  for(unsigned int iphot=0; iphot<patPhotons->size(); ++iphot) {
    
    const pat::Photon *p1 = &((*patPhotons)[iphot]);
    
    double eta          = p1->eta();
    double sigEtaEta    = p1->sigmaIetaIeta();

    double hadOverEm2012 = p1->hadTowOverEm();

    double pfChargedIso = p1->userFloat("pfChargedIsoAlt");
    double pfNeutralIso = p1->userFloat("pfNeutralIsoAlt");
    double pfGammaIso   = p1->userFloat("pfGammaIsoAlt");
    
    double pfChargedIsoPURel = p1->userFloat("pfChargedPURel");
    double pfNeutralIsoPURel = p1->userFloat("pfNeutralPURel");
    double pfGammaIsoPURel   = p1->userFloat("pfGammaPURel");
    
    double pfChargedIsoPU = p1->userFloat("pfChargedPU");
    double pfNeutralIsoPU = p1->userFloat("pfNeutralPU");
    double pfGammaIsoPU   = p1->userFloat("pfGammaPU");
    
    bool   isPixel      = p1->hasPixelSeed();
    bool   passElecVeto = p1->userInt("passElectronConvVeto");
    bool   hOverE2012  = ( hadOverEm2012 < p1->userFloat("hadTowOverEmTightCut"));
    bool   kineAcc     = ( p1->et() > photon_ptcut && 
			 ((p1->isEE() && std::fabs(eta) > 1.566 && std::fabs(eta) < 2.5) || 
			  (p1->isEB() && std::fabs(eta) < 1.4442) ) ) ;
    bool   showerShape = sigEtaEta < p1->userFloat("showerShapeTightCut");

    bool   isPhotonPFChargedIso = (pfChargedIsoPU < p1->userFloat("pfChargedTightCut"));
    bool   isPhotonPFNeutralIso = (pfNeutralIsoPU < p1->userFloat("pfNeutralTightCut"));
    bool   isPhotonPFGammaIso   = (pfGammaIsoPU   < p1->userFloat("pfGammaTightCut"));
    bool   isPhotonPFChargedIsoRel = (pfChargedIsoPURel < p1->userFloat("pfChargedRelTightCut"));
    bool   isPhotonPFNeutralIsoRel = (pfNeutralIsoPURel < p1->userFloat("pfNeutralRelTightCut"));
    bool   isPhotonPFGammaIsoRel   = (pfGammaIsoPURel   < p1->userFloat("pfGammaRelTightCut"));
    bool   isPhotonPFIso = ( isPhotonPFChargedIso && isPhotonPFNeutralIso && isPhotonPFGammaIso );
    //bool   isPhoton2012PF  = ( kineAcc && !isPixel && hOverE2012 && showerShape && isPhotonPFIso );
    bool   isPhoton2012PF  = ( kineAcc && passElecVeto && hOverE2012 && showerShape && isPhotonPFIso );

    if (debug_) {
      std::cout<<
	"photon("     <<iphot             <<"):("<< 
	p1->pt()<<","<<p1->et()<<","<<p1->eta()<<","<<p1->phi()<<") - (";
      if (p1->genPhoton()) {
	std::cout<<
	  p1->genPhoton()->pdgId()<<","<<
	  p1->genPhoton()->status()<<","<<
	  p1->genPhoton()->mother()->pdgId()<<") - ";
      }
      else 
	std::cout<<") - ";
      std::cout<<"isPhoton2012PF("<<isPhoton2012PF      <<") - ";
      //<< //std::endl<<
      if (!kineAcc)
	std::cout<<"kineAcc("     <<kineAcc               <<") - ";//<< 
      if (isPixel)
	std::cout<<"!isPixel("    <<!isPixel              <<") - ";//<< 
      if (!passElecVeto)
	std::cout<<"passElecVeto("<<passElecVeto          <<") - ";//<< 
      if (!hOverE2012)
	std::cout<<"hOverE2012("  <<p1->userFloat("hadTowOverEmTightCut")<<","<<hOverE2012            <<","<<hadOverEm2012<<") - ";//<< 
      if (!showerShape)
	std::cout<<"showerShape(" <<p1->userFloat("showerShapeTightCut")<<","<<showerShape           <<","<<sigEtaEta    <<") - ";//<<// std::endl<<
      if (!isPhotonPFChargedIso || !isPhotonPFChargedIsoRel){
	std::cout<<"pfChargedIso("<<pfChargedIso<<","<<p1->userFloat("pfChargedTightCut")<<","<<p1->userFloat("pfChargedRelTightCut")<<") - "<<
	  "pfChargedIso("<<isPhotonPFChargedIso<<","<<pfChargedIsoPU<<") - "<<
	  "pfChargedIsoRel("<<isPhotonPFChargedIsoRel<<","<<pfChargedIsoPURel<<") - ";//<<
      }
      if (!isPhotonPFNeutralIso || !isPhotonPFNeutralIsoRel){
	std::cout<<"pfNeutralIso("<<pfNeutralIso<<","<<p1->userFloat("pfNeutralTightCut")<<","<<p1->userFloat("pfNeutralRelTightCut")<<") - "<<
	  "pfNeutralIso("<<isPhotonPFNeutralIso<<","<<pfNeutralIsoPU<<") - "<<
	  "pfNeutralIsoRel("<<isPhotonPFNeutralIsoRel<<","<<pfNeutralIsoPURel<<") - ";//<<
      }
      if (!isPhotonPFGammaIso || !isPhotonPFGammaIsoRel) {
	std::cout<<"pfGammaIso("  <<pfGammaIso  <<","<<p1->userFloat("pfGammaTightCut")  <<","<<p1->userFloat("pfGammaRelTightCut")  <<") - "<<
	  "pfGammaIso("  <<isPhotonPFGammaIso  <<","<<pfGammaIsoPU  <<") - "<<
	  "pfGammaIsoRel("  <<isPhotonPFGammaIsoRel  <<","<<pfGammaIsoPURel  <<")";
      }
      //"isPhotonPFIso       ("<<isPhotonPFIso       <<") - "<<
      std::cout<<std::endl;
    }
   
    if( isPhoton2012PF ) IsoPhotons.push_back(p1);
  } // loop over patPhotons

  // require atleast one isolated photon else return
  if (IsoPhotons.size()<1) {
    std::cout << "No isolated photons found : Event Rejected : ( run, event, lumi ) " 
    << run << " " << event << " " << lumi << std::endl;
    return;
    }

  double photon_eta = IsoPhotons[0]->eta();
  double photon_phi = IsoPhotons[0]->phi();

  double photon2_eta = -999.;
  double photon2_phi = -999.;
  
  // clean jet collection

  std::vector<const pat::Jet*> Jets; // clean jet collection
  std::vector<const pat::Jet*> BJets; // clean jet collection

  int    bestDRPhot = -1;
  int    bestDRJet  = -1;
  double bestDRMin  =999.0;

  double dRMin=999.0;
  int    dRMinJetIdx = -1;
  double dRMin2=999.0;
  int    dRMinJetIdx2 = -1;
  for(unsigned int i=0; i<recoJets->size(); i++) {
    const pat::Jet *r = &((*recoJets)[i]);
   
    double jet_pt  = r->pt();
    double jet_eta = r->eta();
    double jet_phi = r->phi();
    
    if( jet_pt<30.0 ) continue;
    
    double dR = reco::deltaR(photon_eta, photon_phi, jet_eta, jet_phi);
    if(dR<dRMin) { 
      dRMin       = dR;
      dRMinJetIdx = i;
    }

    if (IsoPhotons.size() > 1) {
      double dR2 = reco::deltaR(photon2_eta, photon2_phi, jet_eta, jet_phi);
      if(dR2<dRMin2) { 
	dRMin2       = dR2;
	dRMinJetIdx2 = i;
      }
    }
  }
  
  //original
  if (dRMin < dRMin2) {
    bestDRMin = dRMin;
    bestDRPhot = 0;
    bestDRJet = dRMinJetIdx;
  }
  else if (dRMin2 < dRMin) {
    bestDRMin = dRMin2;
    bestDRPhot = 1;
    bestDRJet = dRMinJetIdx2;
  }

  //////
  if(debug_ && IsoPhotons.size() > 0) {
    std::cout << "Isolated photons : " << std::endl;
    for( unsigned int iphot=0; iphot<IsoPhotons.size(); iphot++) {
      std::cout << iphot << " " <<IsoPhotons[iphot]->pt() 
		<< " " << IsoPhotons[iphot]->eta()
		<< " " << IsoPhotons[iphot]->phi() 
		<< std::endl;
    }
    std::cout << "Jets("<<recoJets->size()<<") : " << std::endl;
    for(unsigned int i=0; i<recoJets->size(); i++) {
      const pat::Jet *r = &((*recoJets)[i]);
      std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
    }
    //std::cout << "Cleaned Jets("<<Jets.size()<<") : " << std::endl;
    //for(unsigned int i=0; i<Jets.size(); ++i) {
    //  const pat::Jet *r = Jets[i];
    //  std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
    //}
    std::cout << "bJets("<<recobJets->size()<<") : " << std::endl;
    for(unsigned int i=0; i<recobJets->size(); i++) {
      const pat::Jet *r = &((*recobJets)[i]);
      std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
    }
    //std::cout << "Cleaned bJets("<<BJets.size()<<") : " << std::endl;
    //for(unsigned int i=0; i<BJets.size(); ++i) {
    //  const pat::Jet *r = BJets[i];
    //  std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
    //}
  }
  //////

  // if no matching jet found in dRMin<0.1 then reject this event
  //if(dRMinJetIdx<0 || dRMin>0.1 ) return;
  if(bestDRJet<0 || bestDRMin>0.1 ) {
    if (debug_) {
      std::cout << "No jet matched to isolated photon found : Event Rejected : ( run, event, lumi ) " 
		<< run << " " << event << " " << lumi << std::endl;
    }
    return;
  }

  //original
  if (dRMin < dRMin2) {
    bestDRMin = dRMin;
    bestDRPhot = 0;
    bestDRJet = dRMinJetIdx;
  }
  else if (dRMin2 < dRMin) {
    bestDRMin = dRMin2;
    bestDRPhot = 1;
    bestDRJet = dRMinJetIdx2;
  }
  

  if ( bestDRPhot == 0 ) {
    RA2Photons.push_back(IsoPhotons[0]);
    if (IsoPhotons.size() > 1)
      RA2Photons.push_back(IsoPhotons[1]);
  }
  if ( bestDRPhot == 1 ) {
    RA2Photons.push_back(IsoPhotons[1]);
    RA2Photons.push_back(IsoPhotons[0]);
    }

  m_nPhotonsIso = RA2Photons.size();
  m_PhotonPt  = RA2Photons[0]->pt();
  m_PhotonEta = RA2Photons[0]->eta();
  for(unsigned int i=0; i<recoJets->size(); ++i) {
    const pat::Jet *r = &((*recoJets)[i]);
    int jj = (int) i;
    if( jj != bestDRJet ) {
      Jets.push_back(r);
      if (r->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898)
	//if (r->bDiscriminator("combinedSecondaryVertexMVABJetTags") > 0.898)
	BJets.push_back(r);
    }
  }
  /*
  // clean b-jet collection
  std::vector<const pat::Jet*> BJets; // clean jet collection

  bestDRPhot = -1;
  bestDRJet  = -1;
  bestDRMin  =999.0;

  dRMin=999.0;
  dRMinJetIdx = -1;
  dRMin2=999.0;
  dRMinJetIdx2 = -1;
  for(unsigned int i=0; i<recobJets->size(); i++) {
    const pat::Jet *r = &((*recobJets)[i]);
   
    double jet_pt  = r->pt();
    double jet_eta = r->eta();
    double jet_phi = r->phi();
    
    if( jet_pt<30.0 ) continue;
    
    double dR = reco::deltaR(photon_eta, photon_phi, jet_eta, jet_phi);
    if(dR<dRMin) { 
      dRMin       = dR;
      dRMinJetIdx = i;
    }

    if (IsoPhotons.size() > 1) {
      double dR2 = reco::deltaR(photon2_eta, photon2_phi, jet_eta, jet_phi);
      if(dR2<dRMin2) { 
	dRMin2       = dR2;
	dRMinJetIdx2 = i;
      }
    }
  }
  
  for(unsigned int i=0; i<recobJets->size(); ++i) {
    const pat::Jet *r = &((*recobJets)[i]);
    int jj = (int) i;
    if( jj != bestDRJet ) BJets.push_back(r);
  }
  */
  //if(debug_) {
  //  std::cout << "Isolated photons : " << std::endl;
  //  for( unsigned int iphot=0; iphot<IsoPhotons.size(); iphot++) {
  //    std::cout << iphot << " " <<IsoPhotons[iphot]->pt() 
  //            << " " << IsoPhotons[iphot]->eta()
  //		<< " " << IsoPhotons[iphot]->phi() 
  //		<< std::endl;
  //  }
  //  std::cout << "Jets("<<recoJets->size()<<") : " << std::endl;
  //  for(unsigned int i=0; i<recoJets->size(); i++) {
  //    const pat::Jet *r = &((*recoJets)[i]);
  //    std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
  //  }
  //  std::cout << "Cleaned Jets("<<Jets.size()<<") : " << std::endl;
  //  for(unsigned int i=0; i<Jets.size(); ++i) {
  //    const pat::Jet *r = Jets[i];
  //    std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
  //  }
  //  std::cout << "bJets("<<recobJets->size()<<") : " << std::endl;
  //  for(unsigned int i=0; i<recobJets->size(); i++) {
  //    const pat::Jet *r = &((*recobJets)[i]);
  //    std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
  //  }
  //  std::cout << "Cleaned bJets("<<BJets.size()<<") : " << std::endl;
  //  for(unsigned int i=0; i<BJets.size(); ++i) {
  //    const pat::Jet *r = BJets[i];
  //    std::cout << i << " " << r->pt()<<" "<<r->eta()<<" "<<r->phi()<<std::endl;
  //  }
  //}

  
  //compare the sizes of the collections computed within the code an with the external module

  m_nJetsPt30Eta50 = 0;
  m_nJetsPt50Eta25 = 0;
  m_bJetsPt30Eta24 = 0;
  m_HT  = 0.0;
  reco::MET::LorentzVector mht(0,0,0,0);
  for(unsigned int i=0; i<Jets.size(); ++i) {
    //const pat::Jet *r = Jets[i];
    if(Jets[i]->pt() > 50.0 && fabs(Jets[i]->eta()) < 2.50) {
      m_nJetsPt50Eta25++;
      m_HT  += Jets[i]->pt();
    }
    if(Jets[i]->pt() > 30.0 && fabs(Jets[i]->eta()) < 5.0) { 
      m_nJetsPt30Eta50++;
      mht  -= Jets[i]->p4();
    }
  }
  reco::MET MHT = reco::MET(mht, reco::MET::Point());
  m_MHT =  MHT.pt();

  //count the photon cleaned bJets
  for(unsigned int i=0; i<BJets.size(); ++i) {
    //const pat::Jet *r = BJets[i];
    if(BJets[i]->pt() > 30.0 && fabs(BJets[i]->eta()) < 2.40) {
      m_bJetsPt30Eta24++;
    }
  }
  const pat::Jet  *r1, *r2, *r3;
  m_dPhi1 = -1.0;
  m_dPhi2 = -1.0;
  m_dPhi3 = -1.0;
  if(m_nJetsPt30Eta50 >= 1) {
    r1 = Jets[0];
    m_dPhi1 = fabs(reco::deltaPhi(r1->phi(),MHT.phi()));
    if(m_nJetsPt30Eta50 >= 2) {
      r2 = Jets[1];
      m_dPhi2 = fabs(reco::deltaPhi(r2->phi(),MHT.phi()));
      
      if(m_nJetsPt30Eta50 >= 3) {
	r3 = Jets[2];
	m_dPhi3 = fabs(reco::deltaPhi(r3->phi(),MHT.phi()));
      }
    }
  }

  /////Trigger information
  m_Photon70PFMET100 = true;
  m_Photon70PFHT400 = true;
  m_Photon70PFNoPUHT400 = true;
  m_Photon135 = true;
  m_Photon150 = true;

  /******************************************************************
   * Here we do all the HLT related trigger stuff
   *
   *
   ******************************************************************/
  // Get the HLT results and check validity
  
  if (data_) {
    m_Photon70PFMET100 = false;
    m_Photon70PFHT400 = false;
    m_Photon70PFNoPUHT400 = false;
    m_Photon135 = false;
    m_Photon150 = false;

    if (!getHLTfromConfig_) 
      if (processName_=="") {
	Handle<trigger::TriggerEvent> hltEventHandle;
	ev.getByLabel("hltTriggerSummaryAOD", hltEventHandle);
	processName_ = hltEventHandle.provenance()->processName();
	if (debug_)
	  std::cout<<processName_<<std::endl;
      }
    hlTriggerResults_ = InputTag("TriggerResults","",processName_);
    
    edm::LogInfo("HLTEventSelector") << "Using trigger results for InputTag " << hlTriggerResults_;
    
    edm::Handle<edm::TriggerResults> hltHandle;
    ev.getByLabel(hlTriggerResults_, hltHandle);
    
    if ( !hltHandle.isValid() ) {
      edm::LogWarning("HLTEventSelector") << "No trigger results for InputTag " << hlTriggerResults_;
      if (debug_)
	std::cout<<"HLT results not valid"<<std::endl;
      return;
    }
    
    const edm::TriggerNames& trgNames = ev.triggerNames(*hltHandle);
    
    int          prescaleSet = hltConfig.prescaleSet(ev,es);
    
    if (debug_)
      std::cout<<"Prescale set is: "<<prescaleSet<<std::endl;
    for (unsigned int hltnum = 0; hltnum < trgNames.size(); ++hltnum) {
      std::string  tmpName     = trgNames.triggerName(hltnum);
      unsigned int trgIndex    = trgNames.triggerIndex(tmpName);
      int          trgResult   = hltHandle->accept(trgIndex);
      
      if (trgResult > 0) {
	if (tmpName.rfind("HLT_Photon70_CaloIdXL_PFHT400_v") != std::string::npos)
	  m_Photon70PFHT400 = true;
	else if (tmpName.rfind("HLT_Photon70_CaloIdXL_PFNoPUHT400_v") != std::string::npos)
	  m_Photon70PFNoPUHT400 = true;
	else if (tmpName.rfind("HLT_Photon70_CaloIdXL_PFMET100_v") != std::string::npos)
	  m_Photon70PFMET100 = true;
	else if (tmpName.rfind("HLT_Photon135_v") != std::string::npos)
	  m_Photon135 = true;
	else if (tmpName.rfind("HLT_Photon150_v") != std::string::npos)
	  m_Photon150 = true;
      }
    }
  }
  //if (reducedValues)
  reducedValues->Fill();
 
}

void RA2ZInvPhotonTreeMaker::beginJob() {
  //book trees
  BookTree();
}

void RA2ZInvPhotonTreeMaker::endJob() {

}


void RA2ZInvPhotonTreeMaker::BookTree() {

  edm::Service<TFileService> fs;

  reducedValues = fs->make<TTree>( "RA2Values", "Variables for reduced studies" );

  reducedValues->Branch("ra2_HT",  &m_HT,  "m_HT/D" );
  reducedValues->Branch("ra2_MHT", &m_MHT, "m_MHT/D");
  reducedValues->Branch("ra2_Vertices", &m_Vertices, "m_Vertices/I");

  reducedValues->Branch("ra2_dPhi1", &m_dPhi1, "m_dPhi1/D");
  reducedValues->Branch("ra2_dPhi2", &m_dPhi2, "m_dPhi2/D");
  reducedValues->Branch("ra2_dPhi3", &m_dPhi3, "m_dPhi3/D");

  reducedValues->Branch("ra2_PUWt",    &m_PUWt,    "m_PUWt/D");
  reducedValues->Branch("ra2_EventWt", &m_EventWt, "m_EventWt/D");

  reducedValues->Branch("ra2_nPhotonsIso", &m_nPhotonsIso, "m_nPhotonsIso/I");
  reducedValues->Branch("ra2_PhotonPt",    &m_PhotonPt,    "m_PhotonPt/D" );
  reducedValues->Branch("ra2_PhotonEta",   &m_PhotonEta,   "m_PhotonEta/D");

  reducedValues->Branch("ra2_Photon70PFMET100",    &m_Photon70PFMET100,    "m_Photon70PFMET100/O" );
  reducedValues->Branch("ra2_Photon70PFHT400",     &m_Photon70PFHT400,     "m_Photon70PFHT400/O" );
  reducedValues->Branch("ra2_Photon70PFNoPUHT400", &m_Photon70PFNoPUHT400, "m_Photon70PFNoPUHT400/O" );
  reducedValues->Branch("ra2_Photon135",           &m_Photon135,           "m_Photon135/O" );
  reducedValues->Branch("ra2_Photon150",           &m_Photon150,           "m_Photon150/O" );

  reducedValues->Branch("ra2_nJetsPt30Eta50", &m_nJetsPt30Eta50, "m_nJetsPt30Eta50/I" );
  reducedValues->Branch("ra2_bJetsPt30Eta24", &m_bJetsPt30Eta24, "m_bJetsPt30Eta24/I");
  reducedValues->Branch("ra2_nJetsPt50Eta25", &m_nJetsPt50Eta25, "m_nJetsPt50Eta25/I" );

  reducedValues->SetAutoSave(1);
}


void  RA2ZInvPhotonTreeMaker::beginRun(edm::Run const& run, edm::EventSetup const& es) {
  bool changed = false;
  if (data_) {
    if (hltConfig.init(run,es,"HLT",changed)) {
      if (changed) {
	edm::LogWarning("RA2ZInvPhotonTreeMaker") << "beginRun: The HLT config has changed!";
      }
    }
    else {
      edm::LogError("TriggerEvent") << " HLT config extraction failure";
    }
  }
}

void RA2ZInvPhotonTreeMaker::endRun(edm::Run const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonTreeMaker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonTreeMaker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvPhotonTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RA2ZInvPhotonTreeMaker);

