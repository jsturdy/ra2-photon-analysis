#ifndef ZINVISIBLEBKGDS_PHOTONS_RA2ZINVTREEMAKERFUNCTIONS_H
#define ZINVISIBLEBKGDS_PHOTONS_RA2ZINVTREEMAKERFUNCTIONS_H

// system include files
#include <memory>
#include <vector>
#include <string>

//// user include files
//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
//#include "FWCore/Framework/interface/Event.h"
//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include <DataFormats/PatCandidates/interface/Photon.h>
#include <DataFormats/PatCandidates/interface/Jet.h>

//// Trigger includes
//#include "DataFormats/Common/interface/TriggerResults.h"
//#include "DataFormats/HLTReco/interface/TriggerEvent.h"
//#include "DataFormats/HLTReco/interface/TriggerObject.h"
//#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
//#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
//#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
//#include "FWCore/Common/interface/TriggerNames.h"
//
//// TFile Service
//#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "CommonTools/UtilAlgos/interface/TFileService.h"
//#include "TTree.h"
//#include "TH1.h"

namespace RA2ZInvTreeMakerFunctions {
  
  //  // public:
  //  void computeBJetVars(const edm::View<pat::Jet> bJets, 
  //		       const reco::MET mhtVec,
  //		       const reco::MET metVec,
  //		       int& nTightBJets, 
  //		       std::vector<double*>* ptValsCSVM,
  //		       std::vector<double*>* etaValsCSVM,
  //		       std::vector<double*>* ptValsCSVT,
  //		       std::vector<double*>* etaValsCSVT,
  //		       double& dphiMHTMinBCSVM, double& dphiMHTMinBCSVT,
  //		       double& dphiMETMinBCSVM, double& dphiMETMinBCSVT
  //		       );
  //
  //  void computeAllJetVars(const edm::View<pat::Jet> allJets,
  //			 const reco::MET mhtVec,
  //			 const reco::MET metVec,
  //			 //int& nTightBJets, 
  //			 std::vector<double*>* ptVals,
  //			 std::vector<double*>* etaVals,
  //			 std::vector<double*>* dphiValsMHT,
  //			 std::vector<double*>* dphiValsMET,
  //			 double& dphiMHTMin, double& dphiMETMin,
  //			 int& nJetsPt30Eta24, int& nJetsPt50Eta24, int& nJetsPt70Eta24
  //			 );
  //  
  //};


  inline void computeBJetVars(const edm::View<pat::Jet> bJets, 
			      const reco::MET mhtVec,
			      const reco::MET metVec,
			      int& nTightBJets, 
			      std::vector<double*>* ptValsCSVM,
			      std::vector<double*>* etaValsCSVM,
			      std::vector<double*>* ptValsCSVT,
			      std::vector<double*>* etaValsCSVT,
			      double& dphiMHTMinBCSVM, double& dphiMHTMinBCSVT,
			      double& dphiMETMinBCSVM, double& dphiMETMinBCSVT
						
			      ) {
    const pat::Jet  *b;
    edm::View<pat::Jet>::const_iterator bjet = bJets.begin();
    for (unsigned int ii = 0; ii < ptValsCSVM->size(); ++ii){
      if (bJets.size()>ii){
	b = &((bJets)[ii]);
	*(ptValsCSVM->at(ii))  = b->pt();
	*(etaValsCSVM->at(ii)) = b->eta();
	if (b->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898) {
	  *(ptValsCSVT->at(ii))  = b->pt();
	  *(etaValsCSVT->at(ii)) = b->eta();
	}
      }
    }
    for (; bjet!= bJets.end(); ++bjet) {
      double tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),mhtVec.phi()));
      if (tmpDPhiB < dphiMHTMinBCSVM)
	dphiMHTMinBCSVM = tmpDPhiB;
      
      tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),metVec.phi()));
      if (tmpDPhiB < dphiMETMinBCSVM)
	dphiMETMinBCSVM = tmpDPhiB;
      if (bjet->bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898) {
	++nTightBJets;
	tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),mhtVec.phi()));
	if (tmpDPhiB < dphiMHTMinBCSVT)
	  dphiMHTMinBCSVT = tmpDPhiB;
      
	tmpDPhiB = fabs(reco::deltaPhi(bjet->phi(),metVec.phi()));
	if (tmpDPhiB < dphiMETMinBCSVT)
	  dphiMETMinBCSVT = tmpDPhiB;
      }
    }
  
  }
  
  inline void computeAllJetVars(const edm::View<pat::Jet> allJets,
				const reco::MET mhtVec,
				const reco::MET metVec,
				//int& nTightBJets, 
				std::vector<double*>* ptVals,
				std::vector<double*>* etaVals,
				std::vector<double*>* dphiValsMHT,
				std::vector<double*>* dphiValsMET,
				double& dphiMHTMin, double& dphiMETMin,
				int& nJetsPt30Eta24, int& nJetsPt50Eta24, int& nJetsPt70Eta24
				) {
    const pat::Jet  *j;
    edm::View<pat::Jet>::const_iterator jet = allJets.begin();
    for (unsigned int ii = 0; ii < ptVals->size(); ++ii){
      if (allJets.size()>ii){
	j = &((allJets)[ii]);
	*(ptVals->at(ii))  = j->pt();
	*(etaVals->at(ii)) = j->eta();
	*(dphiValsMHT->at(ii))  = fabs(reco::deltaPhi(j->phi(),mhtVec.phi()));
	
	if (dphiValsMET->size())
	  *(dphiValsMET->at(ii))  = fabs(reco::deltaPhi(j->phi(),metVec.phi()));
	
      }
    }
    for (; jet!= allJets.end(); ++jet) {
      //++nTightJets;
      if (fabs(jet->eta()) < 2.4) {
	if (jet->pt() > 30)
	  ++nJetsPt30Eta24;
	if (jet->pt() > 50)
	  ++nJetsPt50Eta24;
	if (jet->pt() > 70)
	  ++nJetsPt70Eta24;
      }
      
      double tmpDPhi = fabs(reco::deltaPhi(jet->phi(),
					   mhtVec.phi()));
      if (tmpDPhi < dphiMHTMin)
	dphiMHTMin = tmpDPhi;
      
      if (dphiValsMET->size()) {
	tmpDPhi = fabs(reco::deltaPhi(jet->phi(),metVec.phi()));
	if (tmpDPhi < dphiMETMin)
	  dphiMETMin = tmpDPhi;
      }
    }
  }

  inline void computeGenJetVars(const edm::View<reco::GenJet> genJets,
				const reco::GenMET metVec,
				double& genHT,
				double& genMHT,
				std::vector<double*>* ptVals,
				std::vector<double*>* etaVals,
				std::vector<double*>* dphiValsMHT,
				std::vector<double*>* dphiValsMET,
				double& dphiMHTMin, double& dphiMETMin,
				int& nJetsPt30Eta50, int& nJetsPt50Eta25
				) {

    //loop over the gen jets, compute genHT, genMHT and nGenJets
    const reco::GenJet *j;
    edm::View<reco::GenJet>::const_iterator gjet = genJets.begin();
    //reco::GenJetCollection::const_iterator gjet = genJets.begin();
    double htGenJets(0.);
    reco::MET::LorentzVector mhtGen(0,0,0,0);
    std::vector<const reco::Jet*> gen_mhtJets;
    for (; gjet != genJets.end(); ++gjet) {
      if (gjet->pt() > 30) 
	if (fabs(gjet->eta()) < 5.0) {
	  ++nJetsPt30Eta50 ;
	  gen_mhtJets.push_back(&(*gjet));
	  mhtGen -= gjet->p4();
	  if (gjet->pt() > 50. && fabs(gjet->eta()) < 2.5) {
	    ++nJetsPt50Eta25;
	    htGenJets += gjet->pt();
	  }
	}
    }
    reco::MET genMHTV = reco::MET(mhtGen, reco::MET::Point());
    genMHT = genMHTV.pt();
    genHT  = htGenJets;

    for (unsigned int ii = 0; ii < ptVals->size(); ++ii){
      if (genJets.size()>ii){
	  j = &((genJets)[ii]);
	  *(ptVals->at(ii))  = j->pt();
	  *(etaVals->at(ii)) = j->eta();
	  *(dphiValsMHT->at(ii))  = fabs(reco::deltaPhi(j->phi(),genMHTV.phi()));
	  
	  if (dphiValsMET->size())
	    *(dphiValsMET->at(ii))  = fabs(reco::deltaPhi(j->phi(),metVec.phi()));
      }
    }

    gjet = genJets.begin();
    for (; gjet != genJets.end(); ++gjet) {
      double tmpDPhi = fabs(reco::deltaPhi(gjet->phi(),
					   genMHTV.phi()));
      if (tmpDPhi < dphiMHTMin)
	dphiMHTMin = tmpDPhi;
      
      if (dphiValsMET->size()) {
	tmpDPhi = fabs(reco::deltaPhi(gjet->phi(),metVec.phi()));
	if (tmpDPhi < dphiMETMin)
	  dphiMETMin = tmpDPhi;
      }
    }
  }
  

  inline void computeJetVarsNoGen(const edm::View<pat::Jet> ra2Jets,
				  const std::vector<reco::GenParticle> gens,const int nGens,const double maxDR,
				  const reco::GenMET metVec,
				  double& genHT, double& genMHT,
				  std::vector<double*>* dphiValsMHT,
				  std::vector<double*>* dphiValsMET,
				  double& dphiMHTMin, double& dphiMETMin,
				  double& boson1MinDR, double& boson2MinDR,
				  int& nJetsPt30Eta50, int& nJetsPt50Eta25
				  ) {
    
    const pat::Jet  *j;
    std::vector<const pat::Jet*> genRemoved_mhtJets;
    std::vector<int> jetIndex;
    std::vector<double> bestDR;
    jetIndex.reserve(nGens);
    bestDR.reserve(nGens);
    for (int i = 0; i < nGens; i++){
      jetIndex.push_back(-1);
      bestDR.push_back(1000.);
    }
    int iJet(0);
    edm::View<pat::Jet>::const_iterator jet = ra2Jets.begin();
    for (; jet!= ra2Jets.end(); ++jet) {
      for (unsigned int i = 0; i < jetIndex.size(); i++){
	double dR = reco::deltaR(gens.at(i).eta(),gens.at(i).phi(),jet->eta(),jet->phi());
	if (dR < bestDR.at(i)) {
	  bestDR.at(i) = dR;
	  jetIndex.at(i) = iJet;
	}
      }
      ++iJet;
    }//presumably found the matched jets

    //now recompute event without those jets
    double htNoGen(0.);
    reco::MET::LorentzVector mhtNoGen(0,0,0,0);
    iJet = 0;  jet = ra2Jets.begin();
    for (; jet!= ra2Jets.end(); ++jet) {
      if (((iJet == jetIndex.at(0) && bestDR.at(0) > maxDR) ||iJet != jetIndex.at(0)) ||
	  (jetIndex.size()>1&&((iJet == jetIndex.at(1) && bestDR.at(1) > maxDR) || iJet != jetIndex.at(1)))) {
	if (jet->pt() > 30 && fabs(jet->eta()) < 5.0) {
	  mhtNoGen -= jet->p4();
	  genRemoved_mhtJets.push_back(&(*jet));
	  ++nJetsPt30Eta50;
	  if (jet->pt() > 50 && fabs(jet->eta()) < 2.5) {
	    htNoGen += jet->pt();
	    ++nJetsPt50Eta25;
	  }
	}
      }
      ++iJet;
    }
    reco::MET ra2MHTNoGen = reco::MET(mhtNoGen, reco::MET::Point());
    genMHT = ra2MHTNoGen.pt();
    genHT  = htNoGen;
    
    for (unsigned int ii = 0; ii < dphiValsMHT->size(); ++ii){
      if (genRemoved_mhtJets.size()>ii){
	j = genRemoved_mhtJets.at(ii);
	*(dphiValsMHT->at(ii))  = fabs(reco::deltaPhi(j->phi(),ra2MHTNoGen.phi()));
	
	if (dphiValsMET->size())
	  *(dphiValsMET->at(ii))  = fabs(reco::deltaPhi(j->phi(),metVec.phi()));
      }
    }
    unsigned int pjet = 0;
    for (;pjet < genRemoved_mhtJets.size();++pjet) {
      double dR = reco::deltaR(gens.at(0).eta(),gens.at(0).phi(),
			       genRemoved_mhtJets.at(pjet)->eta(),
			       genRemoved_mhtJets.at(pjet)->phi());
      if (boson1MinDR > dR)
	boson1MinDR = dR;
      if (gens.size()>1) {
	dR = reco::deltaR(gens.at(1).eta(),gens.at(1).phi(),
			  genRemoved_mhtJets.at(pjet)->eta(),
			  genRemoved_mhtJets.at(pjet)->phi());
	if (boson2MinDR > dR)
	  boson2MinDR = dR;
      }
    }//end loop over cleaned jet collection
  }
  
}
  
#endif
