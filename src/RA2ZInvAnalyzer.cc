// -*- C++ -*-
//
// Package:    RA2ZInvAnalyzer
// Class:      RA2ZInvAnalyzer
// 
/**\class RA2ZInvAnalyzer RA2ZInvAnalyzer.cc SusyAnalysis/RA2ZInvAnalyzer/src/RA2ZInvAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Seema Sharma
//         Created:  Mon Jun 20 12:58:08 CDT 2011
// $Id$
//
//


// system include files
#include <memory>
#include <iomanip>
#include <cmath>

#include "ZInvisibleBkgds/Photons/interface/RA2ZInvAnalyzer.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include <DataFormats/METReco/interface/MET.h>


RA2ZInvAnalyzer::RA2ZInvAnalyzer(const edm::ParameterSet& pset) {

  // read parameters from config file
  debug_          = pset.getParameter<bool>("Debug");
  doGenAnalysis_  = pset.getParameter<bool>("DoGenAnalysis");
  genParticles_   = pset.getParameter<edm::InputTag>("GenParticles");  
  data_           = pset.getParameter<bool>("Data");
  jetSrc_         = pset.getParameter<edm::InputTag>("JetSrc");
  jetHTSrc_       = pset.getParameter<edm::InputTag>("JetHTSource");
  doPUReWeight_   = pset.getParameter<bool>("DoPUReweight");
  puWeightSrc_    = pset.getParameter<edm::InputTag>("PUWeightSource");
  pfRhoSrc_       = pset.getParameter<edm::InputTag>("PFRhoSource");
  
  //minPhoPt_   = pset.getParameter<double>("minPhoPt");
  //ebMaxEta_   = pset.getParameter<double>("ebMaxEta");
  //eeMinEta_   = pset.getParameter<double>("eeMinEta");
  //eeMaxEta_   = pset.getParameter<double>("eeMaxEta");
  //ebSigIeIe_  = pset.getParameter<double>("ebSigIeIe");
  //eeSigIeIe_  = pset.getParameter<double>("eeSigIeIe");
  //recoHoverE_ = pset.getParameter<double>("recoHoverE");

  ra2NJets_         = pset.getParameter<unsigned int>("RA2NJets");
  ra2HT_            = pset.getParameter<double>("RA2HT");
  ra2MHT_           = pset.getParameter<double>("RA2MHT");
  ra2ApplyDphiCuts_ = pset.getParameter<bool>("RA2ApplyDphiCuts");
}

RA2ZInvAnalyzer::~RA2ZInvAnalyzer() {

}

void RA2ZInvAnalyzer::analyze(const edm::Event& ev, const edm::EventSetup& es) {

  using namespace edm;

  // get event-id
  unsigned int event = (ev.id()).event();
  unsigned int run   = (ev.id()).run();
  unsigned int lumi  =  ev.luminosityBlock();

  // reject event if there an isolated electron or muon present
  edm::Handle< std::vector<pat::Electron> > patElectrons;
  ev.getByLabel("patElectronsPFIDIso", patElectrons);

  edm::Handle< std::vector<pat::Muon> > patMuons;
  ev.getByLabel("patMuonsPFIDIso", patMuons);
  
  if (patElectrons->size() != 0 || patMuons->size() != 0) { 
    //std::cout << "Isolated Lepton found : Event Rejected : ( run, event, lumi ) " 
    //<< run << " " << event << " " << lumi << std::endl;
    return;
  }
  
  // get vertices
  edm::Handle< std::vector<reco::Vertex> > Vertices;
  ev.getByLabel("offlinePrimaryVertices", Vertices);
  int nVertices = Vertices->size();
  
  // get jetcollection
  edm::Handle< std::vector<pat::Jet> > recoJets;
  ev.getByLabel(jetSrc_, recoJets); 
  
  edm::Handle<edm::View<pat::Jet> > jets;
  ev.getByLabel(jetSrc_, jets);

  edm::Handle<edm::View<pat::Jet> > jetsHT;
  ev.getByLabel(jetHTSrc_, jetsHT);

  // if MC, do PU reweighting
  double pu_event_wt = 1.0;
  edm::Handle<double> puweight;
  if( doPUReWeight_ ) {
    ev.getByLabel(puWeightSrc_, puweight);
    pu_event_wt = *puweight;
  }
  h_puWeight ->Fill(pu_event_wt);
  h_Vertices ->Fill(nVertices);
  h_VerticesReWeighted ->Fill(nVertices, pu_event_wt);

  // get pfRho for pileup correction to isolation
  edm::Handle<double> pfRho;
  ev.getByLabel(pfRhoSrc_, pfRho);
  double pf_event_rho = *pfRho;
  
  if( nVertices<60 )  h_pfRho[nVertices]->Fill(pf_event_rho);
  
  double areaR03       = 0.3 * 0.3 * 4.0*atan(1.0);
  double areaR04       = 0.4 * 0.4 * 4.0*atan(1.0);
  
  // do the gen analysis if MC events
  //if(! ev.isRealData()) {
  //  edm::Handle<reco::GenParticleCollection> genParticles;
  //  ev.getByLabel(genParticles_, genParticles);
  //  doGenAnalysis(genParticles, jetsHT, jets, IsoPhotons, pu_event_wt);
  //}

  // require atleast one GenBoson else return
  //if (IsoPhotons.size()<1) return;

  int n_jets_pt30 = 0, n_jets_pt50eta25 = 0;
  double ht  = 0.0;
  reco::MET::LorentzVector mht(0,0,0,0);
  for(unsigned int i=0; i<jets->size(); ++i) {
    const pat::Jet r = (*jets)[i];
    if((*jets)[i].pt() > 50.0 && fabs((*jets)[i].eta()) < 2.50) {
      n_jets_pt50eta25++;
      ht  += (*jets)[i].pt();
    }
    if((*jets)[i].pt() > 30.0 && fabs((*jets)[i].eta()) < 5.0) { 
      n_jets_pt30++;
      mht  -= (*jets)[i].p4();
    }
  }
  reco::MET MHT = reco::MET(mht, reco::MET::Point());
  double mht_val =  MHT.pt();

  double dphi_mht_j1 = 0, dphi_mht_j2 = 0, dphi_mht_j3 = 0;
  if(n_jets_pt30 >= 3) {
    dphi_mht_j1 = fabs(reco::deltaPhi((*jets)[0].phi(),MHT.phi()));
    dphi_mht_j2 = fabs(reco::deltaPhi((*jets)[1].phi(),MHT.phi()));
    dphi_mht_j3 = fabs(reco::deltaPhi((*jets)[2].phi(),MHT.phi()));
  }

  // do RA2 analysis here
  //if(n_jets_pt50eta25 >= 3) {
  h_PreRA2_Vertices ->Fill(nVertices);
  h_PreRA2_VerticesReWeighted ->Fill(nVertices, pu_event_wt);

  h_PreRA2_NJets_Pt30     ->Fill(n_jets_pt30, pu_event_wt);
  h_PreRA2_NJets_Pt50Eta25->Fill(n_jets_pt50eta25, pu_event_wt);
  if(n_jets_pt30 >= 3) {
    h_PreRA2_DPhiMHTJet1->Fill(dphi_mht_j1, pu_event_wt);
    h_PreRA2_DPhiMHTJet2->Fill(dphi_mht_j2, pu_event_wt);
    h_PreRA2_DPhiMHTJet3->Fill(dphi_mht_j3, pu_event_wt);
  }
  h_PreRA2_HT ->Fill(ht , pu_event_wt);
  h_PreRA2_MHT->Fill(mht_val, pu_event_wt);
  h_PreRA2_MEff->Fill(ht+mht_val, pu_event_wt);
  
  if(n_jets_pt50eta25 >= ra2NJets_) {
    
    //if( ht> 500.0) {
    if( ht> ra2HT_) {
      
      //if(mht_val > 200.0) {
      if(mht_val > ra2MHT_) {
	h_PreDPhi_Vertices ->Fill(nVertices);
	h_PreDPhi_VerticesReWeighted ->Fill(nVertices, pu_event_wt);

	h_PreDPhi_NJets_Pt30     ->Fill(n_jets_pt30, pu_event_wt);
	h_PreDPhi_NJets_Pt50Eta25->Fill(n_jets_pt50eta25, pu_event_wt);
	if(n_jets_pt30 >= 3) {
	  h_PreDPhi_DPhiMHTJet1->Fill(dphi_mht_j1, pu_event_wt);
	  h_PreDPhi_DPhiMHTJet2->Fill(dphi_mht_j2, pu_event_wt);
	  h_PreDPhi_DPhiMHTJet3->Fill(dphi_mht_j3, pu_event_wt);
	}
	h_PreDPhi_HT ->Fill(ht , pu_event_wt);
	h_PreDPhi_MHT->Fill(mht_val, pu_event_wt);
	h_PreDPhi_MEff->Fill(ht+mht_val, pu_event_wt);
	
	if( !ra2ApplyDphiCuts_ || (ra2ApplyDphiCuts_ && (dphi_mht_j1 > 0.5 && dphi_mht_j2 > 0.5 && dphi_mht_j3 > 0.3)) ) {
	  
	  h_RA2_Vertices ->Fill(nVertices);
	  h_RA2_VerticesReWeighted ->Fill(nVertices, pu_event_wt);
	  
	  if( debug_) {
	    std::cout << "RA2 " <<std::setw(9)<<run<<":"<<lumi<<":"<<event
		      << std::setw(5)<<" == MHT::"<<mht_val<<std::setw(5)<<" == HT::"<<ht
		      << std::setprecision(4)
		      << std::endl;
	  }
	  
	  //	  for(int ipt=0; ipt<NPhotPtBins-1; ipt++) {
	  //	    for(int ieta=0; ieta<NPhotEtaBins-1; ieta++) {
	  //	      if( photon_pt>PhotPtVal[ipt] && photon_pt<PhotPtVal[ipt+1] ) {
	  //		if( std::fabs(photon_eta)>PhotEtaVal[ieta] && std::fabs(photon_eta)<PhotEtaVal[ieta+1] ) {
	  //		  NPhotonPtEta_RA2[ipt][ieta] += pu_event_wt;
	  //		}
	  //	      }
	  //	    }
	  //	  }

	  // get the event yield in inclusive bins
	  if( ht>350.0 && mht_val>200.0 ) {
	    h_RA2_NJets_HT350_MHT200 ->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT350_MHT200    ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT350_MHT200   ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT350_MHT200  ->Fill(ht+mht_val,       pu_event_wt);
	  }
	  if( ht>800.0 && mht_val>200.0 ) {
	    h_RA2_NJets_HT800_MHT200 ->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT800_MHT200    ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT800_MHT200   ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT800_MHT200  ->Fill(ht+mht_val,       pu_event_wt);
	  }
	  if( ht>800.0 && mht_val>500.0 ) {
	    h_RA2_NJets_HT800_MHT500 ->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT800_MHT500    ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT800_MHT500   ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT800_MHT500  ->Fill(ht+mht_val,       pu_event_wt);
	  }
	  if( ht>500.0 && mht_val>350.0 ) {
	    h_RA2_NJets_HT500_MHT350 ->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT500_MHT350    ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT500_MHT350   ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT500_MHT350  ->Fill(ht+mht_val,       pu_event_wt);
	  }
	  if( ht>500.0 && mht_val>200.0 ) {
	    h_RA2_NJets_HT500_MHT200 ->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT500_MHT200    ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT500_MHT200   ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT500_MHT200  ->Fill(ht+mht_val,       pu_event_wt);
	  }
	  if( ht>1000.0 && mht_val>600.0 ) {
	    h_RA2_NJets_HT1000_MHT600->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT1000_MHT600   ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT1000_MHT600  ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT1000_MHT600 ->Fill(ht+mht_val,       pu_event_wt);
	  }
	  if( ht>1200.0 && mht_val>400.0 ) {
	    h_RA2_NJets_HT1200_MHT400->Fill(n_jets_pt50eta25, pu_event_wt);
	    h_RA2_HT_HT1200_MHT400   ->Fill(ht,               pu_event_wt);
	    h_RA2_MHT_HT1200_MHT400  ->Fill(mht_val,          pu_event_wt);
	    h_RA2_MEff_HT1200_MHT400 ->Fill(ht+mht_val,       pu_event_wt);
	  }
	  
	  
	  // get event yields in exclusive bins
	  h_RA2_HTMHT_Excl->Fill( ht, mht_val, pu_event_wt);
	  int ibin = h_RA2_HTMHT_Excl->FindBin(ht, mht_val);
	  int ihtBin, imhtBin, iZ;
	  h_RA2_HTMHT_Excl->GetBinXYZ(ibin, ihtBin, imhtBin, iZ);
	  if( ihtBin >h_RA2_HTMHT_Excl->GetNbinsX() ) ihtBin =h_RA2_HTMHT_Excl->GetNbinsX();
	  if( imhtBin>h_RA2_HTMHT_Excl->GetNbinsY() ) imhtBin=h_RA2_HTMHT_Excl->GetNbinsY();
	  //std::cout << "ibin "<< ibin << " " << ihtBin << " " << imhtBin << std::endl;


	  // get event yield in inclusive bins
	  int ibinInc = h_RA2_HTMHT_Incl->FindBin(ht, mht_val);
	  int ihtBinInc, imhtBinInc, iZInc;
	  h_RA2_HTMHT_Incl->GetBinXYZ(ibinInc, ihtBinInc, imhtBinInc, iZInc);
	  if( ihtBinInc >h_RA2_HTMHT_Incl->GetNbinsX() ) ihtBinInc =h_RA2_HTMHT_Incl->GetNbinsX();
	  if( imhtBinInc>h_RA2_HTMHT_Incl->GetNbinsY() ) imhtBinInc=h_RA2_HTMHT_Incl->GetNbinsY();
	  for(int iht=1; iht<=ihtBinInc; iht++) {
	    for(int imht=1; imht<=imhtBinInc; imht++) {
	      double binHT  = h_RA2_HTMHT_Incl->GetXaxis()->GetBinCenter(iht);
	      double binMHT = h_RA2_HTMHT_Incl->GetYaxis()->GetBinCenter(imht);
	      h_RA2_HTMHT_Incl->Fill(binHT, binMHT, pu_event_wt);
	    }
	  }


	  // just some monitoring plots
	  h_RA2_NJets_Pt30     ->Fill(n_jets_pt30, pu_event_wt);
	  h_RA2_NJets_Pt50Eta25->Fill(n_jets_pt50eta25, pu_event_wt);
	  if(n_jets_pt30 >= 3) {
	    h_RA2_DPhiMHTJet1->Fill(dphi_mht_j1, pu_event_wt);
	    h_RA2_DPhiMHTJet2->Fill(dphi_mht_j2, pu_event_wt);
	    h_RA2_DPhiMHTJet3->Fill(dphi_mht_j3, pu_event_wt);
	  }
	  h_RA2_HT ->Fill(ht , pu_event_wt);
	  h_RA2_MHT->Fill(mht_val, pu_event_wt);
	  h_RA2_MEff->Fill(ht+mht_val, pu_event_wt);
	  
	} 
	
      }

    }
  }
}

void RA2ZInvAnalyzer::doGenAnalysis(edm::Handle<reco::GenParticleCollection>& genParticles, 
					  edm::Handle<edm::View<pat::Jet> >& jetsHT,
					  edm::Handle<edm::View<pat::Jet> >& jetsMHT,
					  std::vector<const pat::Photon*>    IsoPhotons,
					  double pu_event_wt) {
  pu_event_wt = 1.0;

  bool electronMistag= false;
  unsigned int ipart = 0;

  for(reco::GenParticleCollection::const_iterator part = genParticles->begin(); part<genParticles->end(); ++part){

    const reco::GenParticle *mother = (reco::GenParticle *) part->mother();

    if( part->pt()>10.0 ) {

      // definition of direct photon 
      bool directPhoton    = ( std::fabs(part->pdgId()) == 22 && part->status() == 3 &&
			       (std::fabs(mother->pdgId())<10 || mother->pdgId()==21 )  ) ;

      // definition of secondary photon or decay photons
      bool secondaryPhoton = ( std::fabs(part->pdgId()) == 22 && part->status() == 1 && std::fabs(mother->pdgId())>100 );

      // definition of fragmentation photon
      bool fragmenPhoton   = ( std::fabs(part->pdgId()) == 22 && part->status() == 1 && (std::fabs(mother->pdgId())<10 || mother->pdgId()==21 )  );
      
      // definition of electron mistags
      // to be done
      
      // consider only direct photons here
      if( directPhoton ) {
	
	double phot_pt  = part->pt();
	double phot_eta = part->eta();
	double phot_phi = part->phi();
	
	// check if genPhoton is within acceptance 
	//why these |eta| cuts and why this pT cut?
	bool kineAccPhot = ( phot_pt>10.0 && 
			     ((std::fabs(phot_eta)<1.4442) || (std::fabs(phot_eta)>1.566 && 
							      std::fabs(phot_eta)<2.5)) );
	///these are our nominal cuts
	//bool kineAccPhot = ( phot_pt>100.0 && 
	//		     ((std::fabs(phot_eta)<1.4442) || (std::fabs(phot_eta)>1.566 && 
	//						     std::fabs(phot_eta)<2.5)) );
	
	// find a matching recoPhoton to check efficiency
	bool   foundPhotonMatch = false;
	double dRMinPhot=99.0;
	int    dRMinPhotIdx = -1;
	for( unsigned int iphot=0; iphot<IsoPhotons.size(); iphot++) {	      
	  double rphot_pt  = IsoPhotons[iphot]->pt();
	  double rphot_eta = IsoPhotons[iphot]->eta();
	  double rphot_phi = IsoPhotons[iphot]->phi();
	  double dR = reco::deltaR(rphot_eta, rphot_phi, phot_eta, phot_phi);
	  if( dR<dRMinPhot) {
	    dRMinPhot    = dR;
	    dRMinPhotIdx = iphot;
	  }
	}
	if( dRMinPhotIdx>=0 && dRMinPhot<0.1 ) {
	  foundPhotonMatch = true;
	} else {
	  //std::cout << "not Matching to leading Phot " << " dRMin " << dRMin << " NRecoPhotons " << IsoPhotons.size() << std::endl;	    
	}

	// find & remove closest reconstructed jet
	double dRMinJet    = 99.0;
	int    dRMinJetIdx = -1;
	for(unsigned int ijet=0; ijet<jetsMHT->size(); ijet++) {

	  const pat::Jet *r = &((*jetsMHT)[ijet]);
	  double jet_pt  = r->pt();
	  double jet_eta = r->eta();
	  double jet_phi = r->phi();	  

	  if( jet_pt>30.0 ) {
	    double dR = reco::deltaR(phot_eta, phot_phi, jet_eta, jet_phi);
	    if(dR<dRMinJet) { 
	      dRMinJet    = dR;
	      dRMinJetIdx = ijet;
	    }
	  }
	}
	
	// clean jet collection	
	std::vector<const pat::Jet*> jetsClean; 
	if(dRMinJetIdx>=0 && dRMinJet<0.1 ) {
	  for(unsigned int i=0; i<jetsMHT->size(); ++i) {
	    const pat::Jet *r = &((*jetsMHT)[i]);
	    int jj = (int) i;
	    if( jj != dRMinJetIdx ) jetsClean.push_back(r);
	  }
	}
	
	// calculate HT/MHT here
	int n_jets_pt50eta25 = 0;
	double ht  = 0.0;
	double ht_jet30 = 0.0;
	reco::MET::LorentzVector mhtV(0,0,0,0);

	bool goodJetId = true;
	for(unsigned int i=0; i<jetsClean.size(); ++i) {
	  const pat::Jet *r = jetsClean[i];
	  if(jetsClean[i]->pt() > 50.0 && std::fabs(jetsClean[i]->eta())<2.50) {
	    n_jets_pt50eta25++;
	    ht  += jetsClean[i]->pt();
	  }
	  if(jetsClean[i]->pt() > 30.0 && std::fabs(jetsClean[i]->eta())<5.0) { 
	    if( jetsClean[i]->neutralHadronEnergyFraction()>0.90 || jetsClean[i]->photonEnergyFraction()/jetsClean[i]->jecFactor(0)>0.95 ) goodJetId=false; 
	    mhtV     -= jetsClean[i]->p4();
	    ht_jet30 += jetsClean[i]->pt();
	  }
	}
	double mht    = mhtV.pt();
	double mhtPhi = mhtV.phi();
	
	h_GenBoson_NJets_Pt30     ->Fill( jetsMHT->size() );
	h_GenBoson_NJets_Pt50Eta25->Fill( n_jets_pt50eta25);
	h_GenBoson_PtVsHT         ->Fill( ht_jet30, part->pt() );
	h_GenBoson_MHTVsHT        ->Fill( ht_jet30, mht );
	h_GenBoson_Pt             ->Fill( part->pt() );
	h_GenBoson_Eta            ->Fill( part->eta() );
	h_GenBoson_Phi            ->Fill( part->phi() );
	h_GenBoson_HT             ->Fill( ht );
	h_GenBoson_MHT            ->Fill( mht );



	// do RA2 analysis here
	double dphi_mht_j1 = 0, dphi_mht_j2 = 0, dphi_mht_j3 = 0;
	if( jetsClean.size()> 2 ){
	  dphi_mht_j1 = std::fabs(reco::deltaPhi(jetsClean[0]->phi(),mhtPhi));
	  dphi_mht_j2 = std::fabs(reco::deltaPhi(jetsClean[1]->phi(),mhtPhi));
	  dphi_mht_j3 = std::fabs(reco::deltaPhi(jetsClean[2]->phi(),mhtPhi));

	}
	if(goodJetId && n_jets_pt50eta25 >= ra2NJets_) {
   
	  if( ht> ra2HT_) {
      
	    if(mht > ra2MHT_) {

	      h_GenBoson_preRA2_dPhiMHTJet1->Fill(dphi_mht_j1);
	      h_GenBoson_preRA2_dPhiMHTJet2->Fill(dphi_mht_j2);
	      h_GenBoson_preRA2_dPhiMHTJet3->Fill(dphi_mht_j3);
	
	      if( !ra2ApplyDphiCuts_ || (ra2ApplyDphiCuts_ && (dphi_mht_j1 > 0.5 && dphi_mht_j2 > 0.5 && dphi_mht_j3 > 0.3)) ) {
		h_GenBoson_RA2_Pt             ->Fill( part->pt() );
		h_GenBoson_RA2_Eta            ->Fill( part->eta() );
		h_GenBoson_RA2_Phi            ->Fill( part->phi() );
		h_GenBoson_RA2_HT             ->Fill( ht );
		h_GenBoson_RA2_MHT            ->Fill( mht );
		
		// get event yields in exclusive bins
		h_GenBosonTotal_RA2_HTMHT_Excl->Fill( ht, mht);

		int ibinInc = h_GenBosonTotal_RA2_HTMHT_Incl->FindBin(ht, mht);
		int ihtBinInc, imhtBinInc, iZInc;
		h_GenBosonTotal_RA2_HTMHT_Incl->GetBinXYZ(ibinInc, ihtBinInc, imhtBinInc, iZInc);
		if( ihtBinInc >h_GenBosonTotal_RA2_HTMHT_Incl->GetNbinsX() ) 
		  ihtBinInc =h_GenBosonTotal_RA2_HTMHT_Incl->GetNbinsX();
		if( imhtBinInc>h_GenBosonTotal_RA2_HTMHT_Incl->GetNbinsY() ) 
		  imhtBinInc=h_GenBosonTotal_RA2_HTMHT_Incl->GetNbinsY();
		for(int iht=1; iht<=ihtBinInc; iht++) {
		  for(int imht=1; imht<=imhtBinInc; imht++) {
		    double binHT  = h_GenBosonTotal_RA2_HTMHT_Incl->GetXaxis()->GetBinCenter(iht);
		    double binMHT = h_GenBosonTotal_RA2_HTMHT_Incl->GetYaxis()->GetBinCenter(imht);
		    h_GenBosonTotal_RA2_HTMHT_Incl->Fill(binHT, binMHT, pu_event_wt);
		  }
		}
		if( kineAccPhot ) {
		  h_GenBosonKineAcc_RA2_HTMHT_Excl->Fill( ht, mht);
		  for(int iht=1; iht<=ihtBinInc; iht++) {
		    for(int imht=1; imht<=imhtBinInc; imht++) {
		      double binHT  = h_GenBosonKineAcc_RA2_HTMHT_Incl->GetXaxis()->GetBinCenter(iht);
		      double binMHT = h_GenBosonKineAcc_RA2_HTMHT_Incl->GetYaxis()->GetBinCenter(imht);
		      h_GenBosonKineAcc_RA2_HTMHT_Incl->Fill(binHT, binMHT, pu_event_wt);
		    }
		  }
		}
		if( foundPhotonMatch ) {
		  h_GenBosonMatched_RA2_HTMHT_Excl->Fill( ht, mht);
		  for(int iht=1; iht<=ihtBinInc; iht++) {
		    for(int imht=1; imht<=imhtBinInc; imht++) {
		      double binHT  = h_GenBosonMatched_RA2_HTMHT_Incl->GetXaxis()->GetBinCenter(iht);
		      double binMHT = h_GenBosonMatched_RA2_HTMHT_Incl->GetYaxis()->GetBinCenter(imht);
		      h_GenBosonMatched_RA2_HTMHT_Incl->Fill(binHT, binMHT, pu_event_wt);
		    }
		  }
		}
	      }
	    }
	  }
	}
     


      } // if direct photon

      /*
      if( debug_ ) {
	if( directPhoton )    std::cout << "It's a Direct Photon "        << std::endl;
	if( secondaryPhoton ) std::cout << "It's a Secondary Photon "    << std::endl;
	if( fragmenPhoton )   std::cout << "It's a Fragmentation Photon " << std::endl;
	if( directPhoton || secondaryPhoton || fragmenPhoton ) {
	  std::cout << ipart << " " << part->pdgId() << " status " << part->status() 
		    << " pt " << part->pt() << " eta " << part->eta() << " phi " << part->phi()
		    << " mother pt " << mother->pt() << " mother Status " << mother->status()
		    << " mother pdgId " << mother->pdgId()
		    << std::endl;
	}
      }
      */
    } // if a parton with pt>10
  
    ipart++;
  } // looop over gen particles

}







void RA2ZInvAnalyzer::beginJob() {

  double RA2MHTValTemp[NRA2MHTBins] = {200, 350, 500, 600,  800};
  double RA2HTValTemp [NRA2HTBins ] = {350, 500, 800, 1000, 1200, 1400, 1600};
  for(int imht=0; imht<NRA2MHTBins; imht++) RA2MHTVal[imht] = RA2MHTValTemp[imht];
  for(int iht=0;  iht<NRA2HTBins;   iht++)  RA2HTVal [iht]  = RA2HTValTemp[iht];
  
  for(int iht=0;  iht<NRA2HTBins;   iht++){
    for(int imht=0; imht<NRA2MHTBins; imht++) {
      std::cout << "(" << RA2HTVal[iht] << "," << RA2MHTVal[imht] << ")" <<"  ";
    }
    std::cout << std::endl;
  }

  
  double MHTValTemp[NMHTBins] = {200, 250, 300, 350, 400,  500,  600,  800};
  double HTValTemp [NHTBins]  = {350, 400, 500, 600, 800, 1000, 1200, 1400};
  for(int imht=0; imht<NMHTBins; imht++) MHTVal[imht] = MHTValTemp[imht];
  for(int iht=0;  iht<NHTBins;   iht++)  HTVal [iht]  = HTValTemp [iht];
  for(int iht=0;  iht<NHTBins;   iht++){
    for(int imht=0; imht<NMHTBins; imht++) {
      std::cout << "(" << HTVal[iht] << "," << MHTVal[imht] << ")" <<"  ";
    }
    std::cout << std::endl;
  }
 

  //book histograms
  BookHistograms();
}

void RA2ZInvAnalyzer::endJob() {

  // treat the overflows properly 
  for(int imht=1; imht<=h_RA2_HTMHT_Excl->GetNbinsY(); imht++) {
    int    overbin    = h_RA2_HTMHT_Excl->GetNbinsX()+1;
    double overflow   = h_RA2_HTMHT_Excl->GetBinContent(overbin, imht);
    int    lastBinMHT = h_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(imht);
    int    lastBinHT  = h_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(overbin-1);
    h_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);

    int    thisBin    = h_RA2_HTMHT_Excl->GetBin(overbin, imht);
    h_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
  }

  for(int iht=1; iht<=h_RA2_HTMHT_Excl->GetNbinsX(); iht++){
    int    overbin    = h_RA2_HTMHT_Excl->GetNbinsY()+1;
    double overflow   = h_RA2_HTMHT_Excl->GetBinContent(iht, overbin);
    double lastBinMHT = h_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(overbin-1);
    double lastBinHT  = h_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(iht);
    h_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);

    int    thisBin  = h_RA2_HTMHT_Excl->GetBin(iht, overbin);
    h_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
  }


  for(int imht=1; imht<=h_RA2_HTMHT_Excl->GetNbinsY(); imht++) {
    int    overbin    = h_RA2_HTMHT_Excl->GetNbinsX()+1;
    double overflow   = h_RA2_HTMHT_Excl->GetBinContent(overbin, imht);
    int    lastBinMHT = h_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(imht);
    int    lastBinHT  = h_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(overbin-1);
    h_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);

    int    thisBin    = h_RA2_HTMHT_Excl->GetBin(overbin, imht);
    h_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
  }

  if( doGenAnalysis_) { 
    for(int iht=1; iht<=h_GenBosonTotal_RA2_HTMHT_Excl->GetNbinsX(); iht++){
      int    overbin    = h_GenBosonTotal_RA2_HTMHT_Excl->GetNbinsY()+1;
      double overflow   = h_GenBosonTotal_RA2_HTMHT_Excl->GetBinContent(iht, overbin);
      double lastBinMHT = h_GenBosonTotal_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(overbin-1);
      double lastBinHT  = h_GenBosonTotal_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(iht);
      h_GenBosonTotal_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);
      
      int    thisBin  = h_GenBosonTotal_RA2_HTMHT_Excl->GetBin(iht, overbin);
      h_GenBosonTotal_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
    }

    for(int iht=1; iht<=h_GenBosonKineAcc_RA2_HTMHT_Excl->GetNbinsX(); iht++){
      int    overbin    = h_GenBosonKineAcc_RA2_HTMHT_Excl->GetNbinsY()+1;
      double overflow   = h_GenBosonKineAcc_RA2_HTMHT_Excl->GetBinContent(iht, overbin);
      double lastBinMHT = h_GenBosonKineAcc_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(overbin-1);
      double lastBinHT  = h_GenBosonKineAcc_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(iht);
      h_GenBosonKineAcc_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);
      
      int    thisBin  = h_GenBosonKineAcc_RA2_HTMHT_Excl->GetBin(iht, overbin);
      h_GenBosonKineAcc_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
    }

    for(int iht=1; iht<=h_GenBosonMatched_RA2_HTMHT_Excl->GetNbinsX(); iht++){
      int    overbin    = h_GenBosonMatched_RA2_HTMHT_Excl->GetNbinsY()+1;
      double overflow   = h_GenBosonMatched_RA2_HTMHT_Excl->GetBinContent(iht, overbin);
      double lastBinMHT = h_GenBosonMatched_RA2_HTMHT_Excl->GetYaxis()->GetBinCenter(overbin-1);
      double lastBinHT  = h_GenBosonMatched_RA2_HTMHT_Excl->GetXaxis()->GetBinCenter(iht);
      h_GenBosonMatched_RA2_HTMHT_Excl->Fill(lastBinHT, lastBinMHT, overflow);
      
      int    thisBin  = h_GenBosonMatched_RA2_HTMHT_Excl->GetBin(iht, overbin);
      h_GenBosonMatched_RA2_HTMHT_Excl->SetBinContent(thisBin, 0.0);
    }

  }

}


void RA2ZInvAnalyzer::BookHistograms() {

  char hname[200], htit[200];
  h_puWeight = fs->make<TH1F>("h_puWeight", "h_puWeight", 20, 0.0, 2.0);
  h_Vertices = fs->make<TH1F>("h_Vertices", "h_Vertices", 60, 0.0, 60.0);
  h_VerticesReWeighted = fs->make<TH1F>("h_VerticesReWeighted", "h_VerticesReWeighted", 60, 0.0, 60.0);

  for(int i=0; i<60; i++){
    sprintf(hname, "h_pfRho_nVtx_%i", i+1);
    sprintf(htit, "pfRho nVtx==%i", i+1);
    h_pfRho[i] = fs->make<TH1F>(hname, htit, 200,0.0,20.0);

    //sprintf(hname, "h_totalIsolation_nVtx_%i", i+1);
    //sprintf(hname, "totalIsolation nVtx==%i", i+1);
    //h_totalIsolation[i] = fs->make<TH1F>(hname, htit, 100,0.0,30.0);
    //
    //sprintf(hname, "h_totalIsolationRho_nVtx_%i", i+1);
    //sprintf(hname, "totalIsolationRho nVtx==%i", i+1);
    //h_totalIsolationRho[i] = fs->make<TH1F>(hname, htit, 100,0.0,30.0);

  }

  h_PreRA2_Vertices = fs->make<TH1F>("h_PreRA2_Vertices", "h_PreRA2_Vertices", 60, 0.0, 60.0);
  h_PreRA2_VerticesReWeighted = fs->make<TH1F>("h_PreRA2_VerticesReWeighted", "h_PreRA2_VerticesReWeighted", 60, 0.0, 60.0);

  h_PreRA2_NJets_Pt30       = fs->make<TH1F>("h_PreRA2_NJets_Pt30", "h_PreRA2_NJets_Pt30", 20, -0.5, 19.5);
  h_PreRA2_NJets_Pt50Eta25       = fs->make<TH1F>("h_PreRA2_NJets_Pt50Eta25", "h_PreRA2_NJets_Pt50Eta25", 20, -0.5, 19.5);
  h_PreRA2_HT          = fs->make<TH1F>("h_PreRA2_HT",  "h_PreRA2_HT",  100, 0, 5000.0);
  h_PreRA2_MHT         = fs->make<TH1F>("h_PreRA2_MHT", "h_PreRA2_MHT", 150, 0, 1500.0);
  h_PreRA2_MEff        = fs->make<TH1F>("h_PreRA2_MEff","h_PreRA2_MEff",100, 0, 5000.0);
  h_PreRA2_DPhiMHTJet1 = fs->make<TH1F>("h_PreRA2_DPhiMHTJet1", "h_PreRA2_DPhiMHTJet1", 64, -3.2, 3.2);
  h_PreRA2_DPhiMHTJet2 = fs->make<TH1F>("h_PreRA2_DPhiMHTJet2", "h_PreRA2_DPhiMHTJet2", 64, -3.2, 3.2);
  h_PreRA2_DPhiMHTJet3 = fs->make<TH1F>("h_PreRA2_DPhiMHTJet3", "h_PreRA2_DPhiMHTJet3", 64, -3.2, 3.2);

  h_PreDPhi_Vertices = fs->make<TH1F>("h_PreDPhi_Vertices", "h_PreDPhi_Vertices", 60, 0.0, 60.0);
  h_PreDPhi_VerticesReWeighted = fs->make<TH1F>("h_PreDPhi_VerticesReWeighted", "h_PreDPhi_VerticesReWeighted", 60, 0.0, 60.0);

  h_PreDPhi_NJets_Pt30       = fs->make<TH1F>("h_PreDPhi_NJets_Pt30", "h_PreDPhi_NJets_Pt30", 20, -0.5, 19.5);
  h_PreDPhi_NJets_Pt50Eta25       = fs->make<TH1F>("h_PreDPhi_NJets_Pt50Eta25", "h_PreDPhi_NJets_Pt50Eta25", 20, -0.5, 19.5);
  h_PreDPhi_HT          = fs->make<TH1F>("h_PreDPhi_HT",  "h_PreDPhi_HT",  100, 0, 5000.0);
  h_PreDPhi_MHT         = fs->make<TH1F>("h_PreDPhi_MHT", "h_PreDPhi_MHT", 150, 0, 1500.0);
  h_PreDPhi_MEff        = fs->make<TH1F>("h_PreDPhi_MEff","h_PreDPhi_MEff",100, 0, 5000.0);
  h_PreDPhi_DPhiMHTJet1 = fs->make<TH1F>("h_PreDPhi_DPhiMHTJet1", "h_PreDPhi_DPhiMHTJet1", 64, -3.2, 3.2);
  h_PreDPhi_DPhiMHTJet2 = fs->make<TH1F>("h_PreDPhi_DPhiMHTJet2", "h_PreDPhi_DPhiMHTJet2", 64, -3.2, 3.2);
  h_PreDPhi_DPhiMHTJet3 = fs->make<TH1F>("h_PreDPhi_DPhiMHTJet3", "h_PreDPhi_DPhiMHTJet3", 64, -3.2, 3.2);

  // RA2 cuts in addition to identified & isolated photons 
  h_RA2_Vertices = fs->make<TH1F>("h_RA2_Vertices", "h_RA2_Vertices", 60, 0.0, 60.0);
  h_RA2_VerticesReWeighted = fs->make<TH1F>("h_RA2_VerticesReWeighted", "h_RA2_VerticesReWeighted", 60, 0.0, 60.0);

  h_RA2_NJets_Pt30       = fs->make<TH1F>("h_RA2_NJets_Pt30", "h_RA2_NJets_Pt30", 20, -0.5, 19.5);
  h_RA2_NJets_Pt50Eta25  = fs->make<TH1F>("h_RA2_NJets_Pt50Eta25", "h_RA2_NJets_Pt50Eta25", 20, -0.5, 19.5);
  h_RA2_HT          = fs->make<TH1F>("h_RA2_HT",  "h_RA2_HT",  100, 0, 5000.0);
  h_RA2_MHT         = fs->make<TH1F>("h_RA2_MHT", "h_RA2_MHT", 150, 0, 1500.0);
  h_RA2_MEff        = fs->make<TH1F>("h_RA2_MEff","h_RA2_MEff",100, 0, 5000.0);
  h_RA2_DPhiMHTJet1 = fs->make<TH1F>("h_RA2_DPhiMHTJet1", "h_RA2_DPhiMHTJet1", 64, -3.2, 3.2);
  h_RA2_DPhiMHTJet2 = fs->make<TH1F>("h_RA2_DPhiMHTJet2", "h_RA2_DPhiMHTJet2", 64, -3.2, 3.2);
  h_RA2_DPhiMHTJet3 = fs->make<TH1F>("h_RA2_DPhiMHTJet3", "h_RA2_DPhiMHTJet3", 64, -3.2, 3.2);

  // event yield in exclusive bins
  h_RA2_HTMHT_Excl = fs->make<TH2F>("h_RA2_HTMHT_Excl", "h_RA2_HTMHT_Excl", NRA2HTBins-1, RA2HTVal, NRA2MHTBins-1, RA2MHTVal);
  h_RA2_HTMHT_Excl ->Sumw2();
  h_RA2_HTMHT_Excl->GetXaxis()->SetTitle("HT");
  h_RA2_HTMHT_Excl->GetYaxis()->SetTitle("MHT");

//  for(int iht=1; iht<=h_RA2_HTMHT_Excl->GetNbinsX(); iht++){
//    for(int imht=1; imht<=h_RA2_HTMHT_Excl->GetNbinsY(); imht++) {
//      int htxL  = (int)h_RA2_HTMHT_Excl->GetXaxis()->GetBinLowEdge(iht);
//      int htxU  = (int)h_RA2_HTMHT_Excl->GetXaxis()->GetBinUpEdge(iht);
//
//      int mhtyL = (int)h_RA2_HTMHT_Excl->GetYaxis()->GetBinLowEdge(imht);     
//      int mhtyU = (int)h_RA2_HTMHT_Excl->GetYaxis()->GetBinUpEdge(imht);     
//
//      sprintf(hname, "h_RA2_PhotPtEta_Excl_%iHT%i_%iMHT%i", htxL,htxU, mhtyL,mhtyU);
//      //std::cout <<iht <<" "<< imht << " "<<  hname << std::endl;
//      sprintf(htit,  "h_RA2_PhotPtEta_Excl_%iHT%i_%iMHT%i", htxL,htxU, mhtyL,mhtyU);
//      h_RA2_PhotPtEta_Excl[iht-1][imht-1]  = fs->make<TH2F>(hname, htit, NPhotEtaBins-1, PhotEtaVal, NPhotPtBins-1, PhotPtVal);
//    }
//  }
//
  // event yield in inclusive bins
  h_RA2_HTMHT_Incl = fs->make<TH2F>("h_RA2_HTMHT_Incl", "h_RA2_HTMHT_Incl", NHTBins-1, HTVal, NMHTBins-1, MHTVal);
  h_RA2_HTMHT_Incl ->Sumw2();
  h_RA2_HTMHT_Incl->GetXaxis()->SetTitle("HT");
  h_RA2_HTMHT_Incl->GetYaxis()->SetTitle("MHT");

  sprintf(hname, "h_RA2_NJets_HT350_MHT200");
  sprintf(htit,  "NJets: RA2 Baseline (HT>350 & MHT>200 GeV)");
  h_RA2_NJets_HT350_MHT200 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT350_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_NJets_HT800_MHT200");
  sprintf(htit,  "NJets: RA2 High HT (HT>800 & MHT>200 GeV)");
  h_RA2_NJets_HT800_MHT200 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT800_MHT200 ->Sumw2();

  sprintf(hname, "h_RA2_NJets_HT800_MHT500");
  sprintf(htit,  "NJets: RA2 High MHT (HT>800 & MHT>500 GeV)");
  h_RA2_NJets_HT800_MHT500 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT800_MHT500 ->Sumw2();

  sprintf(hname, "h_RA2_NJets_HT500_MHT350");
  sprintf(htit,  "NJets: RA2 Medium HT,MHT (HT>500 & MHT>350 GeV)");
  h_RA2_NJets_HT500_MHT350 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT500_MHT350 ->Sumw2();

  sprintf(hname, "h_RA2_NJets_HT500_MHT200");
  sprintf(htit,  "NJets: RA2 New Baseline (HT>500 & MHT>200 GeV)");
  h_RA2_NJets_HT500_MHT200 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT500_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_NJets_HT1000_MHT600");
  sprintf(htit,  "NJets: RA2 New High MHT (HT>1000 & MHT>600 GeV)");
  h_RA2_NJets_HT1000_MHT600 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT1000_MHT600 ->Sumw2();

  sprintf(hname, "h_RA2_NJets_HT1200_MHT400");
  sprintf(htit,  "NJets: RA2 New High HT (HT>1200 & MHT>400 GeV)");
  h_RA2_NJets_HT1200_MHT400 = fs->make<TH1F>(hname, htit,  10, 0, 10);
  h_RA2_NJets_HT1200_MHT400 ->Sumw2();


  sprintf(hname, "h_RA2_HT_HT350_MHT200");
  sprintf(htit,  "HT: RA2 Baseline (HT>350 & MHT>200 GeV)");
  h_RA2_HT_HT350_MHT200 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT350_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_HT_HT800_MHT200");
  sprintf(htit,  "HT: RA2 High HT (HT>800 & MHT>200 GeV)");
  h_RA2_HT_HT800_MHT200 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT800_MHT200 ->Sumw2();

  sprintf(hname, "h_RA2_HT_HT800_MHT500");
  sprintf(htit,  "HT: RA2 High MHT (HT>800 & MHT>500 GeV)");
  h_RA2_HT_HT800_MHT500 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT800_MHT500 ->Sumw2();

  sprintf(hname, "h_RA2_HT_HT500_MHT350");
  sprintf(htit,  "HT: RA2 Medium HT,MHT (HT>500 & MHT>350 GeV)");
  h_RA2_HT_HT500_MHT350 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT500_MHT350 ->Sumw2();

  sprintf(hname, "h_RA2_HT_HT500_MHT200");
  sprintf(htit,  "HT: RA2 New Baseline (HT>500 & MHT>200 GeV)");
  h_RA2_HT_HT500_MHT200 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT500_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_HT_HT1000_MHT600");
  sprintf(htit,  "HT: RA2 New High MHT (HT>1000 & MHT>600 GeV)");
  h_RA2_HT_HT1000_MHT600 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT1000_MHT600 ->Sumw2();

  sprintf(hname, "h_RA2_HT_HT1200_MHT400");
  sprintf(htit,  "HT: RA2 New High HT (HT>1200 & MHT>400 GeV)");
  h_RA2_HT_HT1200_MHT400 = fs->make<TH1F>(hname, htit,  500, 0, 5000.0);
  h_RA2_HT_HT1200_MHT400 ->Sumw2();

  sprintf(hname, "h_RA2_MHT_HT350_MHT200");
  sprintf(htit,  "MHT: RA2 Baseline (HT>350 & MHT>200 GeV)");
  h_RA2_MHT_HT350_MHT200 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT350_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_MHT_HT800_MHT200");
  sprintf(htit,  "MHT: RA2 High HT (HT>800 & MHT>200 GeV)");
  h_RA2_MHT_HT800_MHT200 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT800_MHT200 ->Sumw2();

  sprintf(hname, "h_RA2_MHT_HT800_MHT500");
  sprintf(htit,  "MHT: RA2 High MHT (HT>800 & MHT>500 GeV)");
  h_RA2_MHT_HT800_MHT500 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT800_MHT500 ->Sumw2();

  sprintf(hname, "h_RA2_MHT_HT500_MHT350");
  sprintf(htit,  "MHT: RA2 Medium HT,MHT (HT>500 & MHT>350 GeV)");
  h_RA2_MHT_HT500_MHT350 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT500_MHT350 ->Sumw2();

  sprintf(hname, "h_RA2_MHT_HT500_MHT200");
  sprintf(htit,  "MHT: RA2 New Baseline (HT>500 & MHT>200 GeV)");
  h_RA2_MHT_HT500_MHT200 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT500_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_MHT_HT1000_MHT600");
  sprintf(htit,  "MHT: RA2 New High MHT (HT>1000 & MHT>600 GeV)");
  h_RA2_MHT_HT1000_MHT600 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT1000_MHT600 ->Sumw2();

  sprintf(hname, "h_RA2_MHT_HT1200_MHT400");
  sprintf(htit,  "MHT: RA2 New High HT (HT>1200 & MHT>400 GeV)");
  h_RA2_MHT_HT1200_MHT400 = fs->make<TH1F>(hname, htit,  100, 0, 1000.0);
  h_RA2_MHT_HT1200_MHT400 ->Sumw2();

  sprintf(hname, "h_RA2_MEff_HT350_MHT200");
  sprintf(htit,  "MEff: RA2 Baseline (HT>350 & MHT>200 GeV)");
  h_RA2_MEff_HT350_MHT200 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT350_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_MEff_HT800_MHT200");
  sprintf(htit,  "MEff: RA2 High HT (HT>800 & MHT>200 GeV)");
  h_RA2_MEff_HT800_MHT200 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT800_MHT200 ->Sumw2();

  sprintf(hname, "h_RA2_MEff_HT800_MHT500");
  sprintf(htit,  "MEff: RA2 High MHT (HT>800 & MHT>500 GeV)");
  h_RA2_MEff_HT800_MHT500 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT800_MHT500 ->Sumw2();

  sprintf(hname, "h_RA2_MEff_HT500_MHT350");
  sprintf(htit,  "MEff: RA2 Medium HT,MHT (HT>500 & MHT>350 GeV)");
  h_RA2_MEff_HT500_MHT350 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT500_MHT350 ->Sumw2();

  sprintf(hname, "h_RA2_MEff_HT500_MHT200");
  sprintf(htit,  "MEff: RA2 New Baseline (HT>500 & MHT>200 GeV)");
  h_RA2_MEff_HT500_MHT200 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT500_MHT200 ->Sumw2(); 

  sprintf(hname, "h_RA2_MEff_HT1000_MHT600");
  sprintf(htit,  "MEff: RA2 New High MHT (HT>1000 & MHT>600 GeV)");
  h_RA2_MEff_HT1000_MHT600 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT1000_MHT600 ->Sumw2();

  sprintf(hname, "h_RA2_MEff_HT1200_MHT400");
  sprintf(htit,  "MEff: RA2 New High HT (HT>1200 & MHT>400 GeV)");
  h_RA2_MEff_HT1200_MHT400 = fs->make<TH1F>(hname, htit,  150, 0.0, 7500.0);
  h_RA2_MEff_HT1200_MHT400 ->Sumw2();

  if(doGenAnalysis_) {
    h_GenBoson_preRA2_dPhiMHTJet1 = fs->make<TH1F>("h_GenBoson_preRA2_dPhiMHTJet1", "h_GenBoson_preRA2_dPhiMHTJet1", 128, -3.2, 3.2);
    h_GenBoson_preRA2_dPhiMHTJet2 = fs->make<TH1F>("h_GenBoson_preRA2_dPhiMHTJet2", "h_GenBoson_preRA2_dPhiMHTJet2", 128, -3.2, 3.2);
    h_GenBoson_preRA2_dPhiMHTJet3 = fs->make<TH1F>("h_GenBoson_preRA2_dPhiMHTJet3", "h_GenBoson_preRA2_dPhiMHTJet3", 128, -3.2, 3.2);

    
    h_GenBoson_PtVsHT  = fs->make<TH2F>("h_GenBoson_PtVsHT",  "GenBoson_Pt  Vs HTJet30",  500, 0.0, 5000.0, 100, 0.0, 1000.0);
    h_GenBoson_PtVsHT->Sumw2();
    h_GenBoson_MHTVsHT = fs->make<TH2F>("h_GenBoson_MHTVsHT", "GenBoson_MHT Vs HTJet30",  500, 0.0, 5000.0, 100, 0.0, 1000.0);
    h_GenBoson_MHTVsHT->Sumw2();
    h_GenBoson_NJets_Pt30       = fs->make<TH1F>("h_GenBoson_NJets_Pt30",      "h_GenBoson_NJets_Pt30", 25, -0.5, 24.5);
    h_GenBoson_NJets_Pt30->Sumw2();
    h_GenBoson_NJets_Pt50Eta25  = fs->make<TH1F>("h_GenBoson_NJets_Pt50Eta25", "h_GenBoson_NJets_Pt50Eta25", 25, -0.5, 24.5);
    h_GenBoson_NJets_Pt50Eta25->Sumw2();
    h_GenBoson_Pt        = fs->make<TH1F>("h_GenBoson_Pt",         "Pt of Gen Boson",            100, 0.0, 1000.0);
    h_GenBoson_Pt->Sumw2();
    h_GenBoson_Eta       = fs->make<TH1F>("h_GenBoson_Eta",        "Eta of Gen Boson",          100, -5.0, 5.0);
    h_GenBoson_Eta->Sumw2();
    h_GenBoson_Phi       = fs->make<TH1F>("h_GenBoson_Phi",        "Phi of Gen Boson",           64, -3.2, 3.2);
    h_GenBoson_Phi->Sumw2();
    h_GenBoson_HT        = fs->make<TH1F>("h_GenBoson_HT",         "HT of Gen Boson",            500, 0.0, 5000.0);
    h_GenBoson_HT->Sumw2();
    h_GenBoson_MHT        = fs->make<TH1F>("h_GenBoson_MHT",       "MHT of Gen Boson",           100, 0.0, 1000.0);
    h_GenBoson_MHT->Sumw2();

    h_GenBoson_RA2_Pt        = fs->make<TH1F>("h_GenBoson_RA2_Pt",         "Pt of Gen Boson",            100, 0.0, 1000.0);
    h_GenBoson_RA2_Pt->Sumw2();
    h_GenBoson_RA2_Eta       = fs->make<TH1F>("h_GenBoson_RA2_Eta",        "Eta of Gen Boson",          100, -5.0, 5.0);
    h_GenBoson_RA2_Eta->Sumw2();
    h_GenBoson_RA2_Phi       = fs->make<TH1F>("h_GenBoson_RA2_Phi",        "Phi of Gen Boson",           64, -3.2, 3.2);
    h_GenBoson_RA2_Phi->Sumw2();
    h_GenBoson_RA2_HT        = fs->make<TH1F>("h_GenBoson_RA2_HT",         "HT of Gen Boson",            500, 0.0, 5000.0);
    h_GenBoson_RA2_HT->Sumw2();
    h_GenBoson_RA2_MHT        = fs->make<TH1F>("h_GenBoson_RA2_MHT",       "MHT of Gen Boson",           100, 0.0, 1000.0);
    h_GenBoson_RA2_MHT->Sumw2();

    // event yield in exclusive bins
    h_GenBosonTotal_RA2_HTMHT_Excl = fs->make<TH2F>("h_GenBosonTotal_RA2_HTMHT_Excl", "h_GenBosonTotal_RA2_HTMHT_Excl", NRA2HTBins-1, RA2HTVal, NRA2MHTBins-1, RA2MHTVal);
    h_GenBosonTotal_RA2_HTMHT_Excl ->Sumw2();
    h_GenBosonTotal_RA2_HTMHT_Excl->GetXaxis()->SetTitle("HT");
    h_GenBosonTotal_RA2_HTMHT_Excl->GetYaxis()->SetTitle("MHT");

    h_GenBosonKineAcc_RA2_HTMHT_Excl = fs->make<TH2F>("h_GenBosonKineAcc_RA2_HTMHT_Excl", "h_GenBosonKineAcc_RA2_HTMHT_Excl", NRA2HTBins-1, RA2HTVal, NRA2MHTBins-1, RA2MHTVal);
    h_GenBosonKineAcc_RA2_HTMHT_Excl ->Sumw2();
    h_GenBosonKineAcc_RA2_HTMHT_Excl->GetXaxis()->SetTitle("HT");
    h_GenBosonKineAcc_RA2_HTMHT_Excl->GetYaxis()->SetTitle("MHT");

    h_GenBosonMatched_RA2_HTMHT_Excl = fs->make<TH2F>("h_GenBosonMatched_RA2_HTMHT_Excl", "h_GenBosonMatched_RA2_HTMHT_Excl", NRA2HTBins-1, RA2HTVal, NRA2MHTBins-1, RA2MHTVal);
    h_GenBosonMatched_RA2_HTMHT_Excl ->Sumw2();
    h_GenBosonMatched_RA2_HTMHT_Excl->GetXaxis()->SetTitle("HT");
    h_GenBosonMatched_RA2_HTMHT_Excl->GetYaxis()->SetTitle("MHT");
    
    // event yield in inclusive bins
    h_GenBosonTotal_RA2_HTMHT_Incl = fs->make<TH2F>("h_GenBosonTotal_RA2_HTMHT_Incl", "h_GenBosonTotal_RA2_HTMHT_Incl", NHTBins-1, HTVal, NMHTBins-1, MHTVal);
    h_GenBosonTotal_RA2_HTMHT_Incl ->Sumw2();
    h_GenBosonTotal_RA2_HTMHT_Incl->GetXaxis()->SetTitle("HT");
    h_GenBosonTotal_RA2_HTMHT_Incl->GetYaxis()->SetTitle("MHT");

    h_GenBosonKineAcc_RA2_HTMHT_Incl = fs->make<TH2F>("h_GenBosonKineAcc_RA2_HTMHT_Incl", "h_GenBosonKineAcc_RA2_HTMHT_Incl", NHTBins-1, HTVal, NMHTBins-1, MHTVal);
    h_GenBosonKineAcc_RA2_HTMHT_Incl ->Sumw2();
    h_GenBosonKineAcc_RA2_HTMHT_Incl->GetXaxis()->SetTitle("HT");
    h_GenBosonKineAcc_RA2_HTMHT_Incl->GetYaxis()->SetTitle("MHT");

    h_GenBosonMatched_RA2_HTMHT_Incl = fs->make<TH2F>("h_GenBosonMatched_RA2_HTMHT_Incl", "h_GenBosonMatched_RA2_HTMHT_Incl", NHTBins-1, HTVal, NMHTBins-1, MHTVal);
    h_GenBosonMatched_RA2_HTMHT_Incl ->Sumw2();
    h_GenBosonMatched_RA2_HTMHT_Incl->GetXaxis()->SetTitle("HT");
    h_GenBosonMatched_RA2_HTMHT_Incl->GetYaxis()->SetTitle("MHT");
  }

}


void  RA2ZInvAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&) {
}

void RA2ZInvAnalyzer::endRun(edm::Run const&, edm::EventSetup const&) {
}

void RA2ZInvAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void RA2ZInvAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RA2ZInvAnalyzer);

