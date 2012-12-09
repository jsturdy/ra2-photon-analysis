#ifndef ZInvisibleBkgds_Photons_plugins_SpecialCleanedMETProducer_h
#define ZInvisibleBkgds_Photons_plugins_SpecialCleanedMETProducer_h

#include <memory>
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/METReco/interface/MET.h"

namespace zinvtools {

  template<class PATObjType, int>
    class SpecialCleanedMETProducer : public edm::EDProducer {
  public:
    /// collection of Object objects
    typedef std::vector<PATObjType> ObjectCollection;
    /// presistent reference to a Object
    typedef edm::Ref<ObjectCollection> ObjectRef;
    /// references to Object collection
    typedef edm::RefProd<ObjectCollection> ObjectRefProd;
    /// vector of references to Object objects all in the same collection
    typedef edm::RefVector<ObjectCollection> ObjectRefVector;
    /// iterator over a vector of references to Object objects all in the same collection
    typedef typename ObjectRefVector::iterator object_iterator;
    
    
    explicit SpecialCleanedMETProducer(const edm::ParameterSet & pset);
    ~SpecialCleanedMETProducer() {};
    
    virtual void produce(edm::Event & ev, const edm::EventSetup & es);
    
  private:
    
    edm::InputTag inputMETLabel_;
    edm::InputTag inputObjectLabel_;
    edm::InputTag inputJetLabel_;
    edm::InputTag inputPFCandidateLabel_;
    double matchDR_, matchDPtRel_;
    int partID_;
    bool debug_;
    //int           objectsToRemove_;
  };
}

//template <class PATObjType>
//zinvtools::SpecialCleanedMETProducer<PATObjType>::~SpecialCleanedMETProducer() {
//}

template <class PATObjType, int objectsToRemove_>
  zinvtools::SpecialCleanedMETProducer<PATObjType, objectsToRemove_>::SpecialCleanedMETProducer(const edm::ParameterSet & pset) :
  inputMETLabel_        ( pset.getParameter<edm::InputTag>("inputMET") ),
  inputObjectLabel_     ( pset.getParameter<edm::InputTag>("inputObjects") ),
  inputJetLabel_        ( pset.getParameter<edm::InputTag>("inputJets") ),
  inputPFCandidateLabel_( pset.getParameter<edm::InputTag>("inputPFCands") ),
  matchDR_              ( pset.getParameter<double>("matchDR") ),
  matchDPtRel_          ( pset.getParameter<double>("matchDPtRel") ),
  partID_               ( pset.getParameter<int>("particleId") ),
  debug_                ( pset.getParameter<bool>("debug") )
  //objectsToRemove_  ( pset.getParameter<edm::InputTag>("objectsToRemove") )
{
  using namespace pat;
  produces<std::vector<reco::MET> >("pfcandjec");
  produces<std::vector<reco::MET> >("pfcand");
  produces<std::vector<reco::MET> >("selected");
}

template <class PATObjType, int objectsToRemove_>
  void 
  zinvtools::SpecialCleanedMETProducer<PATObjType, objectsToRemove_>::produce(edm::Event & ev, const edm::EventSetup & es) 
{
  using namespace pat;

  // read in the met
  edm::Handle<edm::View<reco::MET> > inputMET;
  ev.getByLabel(inputMETLabel_, inputMET);

  // read in the Jets
  edm::Handle<edm::View<pat::Jet> > inputJets;
  ev.getByLabel(inputJetLabel_, inputJets);

  // read in the PFCands
  edm::Handle<edm::View<reco::PFCandidate> > inputPFCands;
  ev.getByLabel(inputPFCandidateLabel_, inputPFCands);

  // read in the objects
  edm::Handle<ObjectCollection > inputObjects;
  ev.getByLabel(inputObjectLabel_, inputObjects);

  // calculate MHT
  std::auto_ptr<std::vector<reco::MET> > metpsel(new std::vector<reco::MET>());
  reco::MET::LorentzVector metsel(0,0,0,0);
  std::auto_ptr<std::vector<reco::MET> > metppfc(new std::vector<reco::MET>());
  reco::MET::LorentzVector metpfc(0,0,0,0);
  std::auto_ptr<std::vector<reco::MET> > metpjec(new std::vector<reco::MET>());
  reco::MET::LorentzVector metjec(0,0,0,0);
  //met = inputMET;
  metsel = (*inputMET)[0].p4();
  metpfc = (*inputMET)[0].p4();
  metjec = (*inputMET)[0].p4();
  //reco::MET::LorentzVector met(*inputMET);
  //reco::MET::LorentzVector met = (*inputMET)[0];
  double jecCorrFactor = 1;
  double corrFactor    = 1;
  double rawCorrFactor = 1;
  //loop over objects to remove from MET
  if (objectsToRemove_ < 0) {
    for (typename ObjectCollection::const_iterator itObj = inputObjects->begin(); itObj != inputObjects->end(); ++itObj) {
      
      if (debug_) {
	std::cout<<"object::pt    eta   phi"<<std::endl;
	std::cout<<"      ::"
	  //<<itObj->particleId()
		 <<"  "<<itObj->pt()
		 <<"  "<<itObj->eta()
		 <<"  "<<itObj->phi()
		 <<std::endl;
      }
      metsel += itObj->p4();
      //loop over jets match deltar to itObj
      for (edm::View<pat::Jet>::const_iterator jet = inputJets->begin();
	   jet != inputJets->end();
	   ++jet) {
	double deltaR = reco::deltaR(jet->eta(),jet->phi(),itObj->eta(), itObj->phi());
	double relPt = fabs(itObj->pt()-jet->pt())/itObj->pt();
	if (deltaR < matchDR_) {
	  pat::Jet uncorrJet = jet->correctedJet("Uncorrected","none",jet->currentJECSet());
	  jecCorrFactor = jet->pt()/uncorrJet.pt();
	  corrFactor    = jet->jecFactor(jet->currentJECLevel(),"none",jet->currentJECSet());
	  rawCorrFactor = jet->jecFactor("Uncorrected","none",jet->currentJECSet());
	  if (debug_) {
	    std::cout<<"found matching jet::pt     corrFacts    dR  #cons   "<<std::endl;
	    std::cout<<"                  ::"
		     <<"  ("<<uncorrJet.pt()<<"," <<jet->pt()<<")"
		     <<"  "<<jet->jecFactor(jet->currentJECLevel(),"none",jet->currentJECSet())
		     <<"  "<<deltaR
		     <<"  "<<jet->nConstituents()
		     <<std::endl;
	  }
	  std::vector<reco::PFCandidatePtr> pfConsts = jet->getPFConstituents();
	  for (unsigned int pfConst = 0;
	       pfConst < pfConsts.size();
	       ++pfConst) {
	    deltaR = reco::deltaR(pfConsts[pfConst]->eta(),pfConsts[pfConst]->phi(),itObj->eta(), itObj->phi());
	    relPt = fabs(itObj->pt()-pfConsts[pfConst]->pt())/itObj->pt();
	    if (deltaR < matchDR_ && relPt < matchDPtRel_) {
	      if (debug_) {
		metpfc += pfConsts[pfConst]->p4();
		metjec += jet->jecFactor(jet->currentJECLevel(),"none",jet->currentJECSet())*(pfConsts[pfConst]->p4());
		std::cout<<"   constituent::"
			 <<pfConsts[pfConst]->particleId()<<"::"
			 <<"  "<<pfConsts[pfConst]->pt()
			 <<"  "<<pfConsts[pfConst]->eta()
			 <<"  "<<pfConsts[pfConst]->phi()
			 <<"  "<<deltaR
			 <<std::endl;
	      }
	    }
	  }
	}
      }
    }
  }
  else if (!(inputObjects->size() < objectsToRemove_))
    for (int obj = 0; obj < objectsToRemove_; ++obj) {

      metsel += (*inputObjects)[obj].p4();
      // match objects to pf candidates in matched jet
      for (edm::View<pat::Jet>::const_iterator jet = inputJets->begin();
	   jet != inputJets->end();
	   ++jet) {
	double deltaR = reco::deltaR(jet->eta(),jet->phi(),(*inputObjects)[obj].eta(), (*inputObjects)[obj].phi());
	double relPt = fabs((*inputObjects)[obj].pt()-jet->pt())/(*inputObjects)[obj].pt();
	if (deltaR < matchDR_) {
	  pat::Jet uncorrJet = jet->correctedJet("Uncorrected","none",jet->currentJECSet());
	  jecCorrFactor = jet->pt()/uncorrJet.pt();
	  corrFactor    = 1/jet->jecFactor(jet->currentJECLevel(),"none",jet->currentJECSet());
	  rawCorrFactor = 1/jet->jecFactor("Uncorrected","none",jet->currentJECSet());
	  std::vector<reco::PFCandidatePtr> pfConsts = jet->getPFConstituents();
	  for (unsigned int pfConst = 0;
	       pfConst < pfConsts.size();
	       ++pfConst) {
	    deltaR = reco::deltaR(pfConsts[pfConst]->eta(),pfConsts[pfConst]->phi(),(*inputObjects)[obj].eta(), (*inputObjects)[obj].phi());
	    relPt = fabs((*inputObjects)[obj].pt()-pfConsts[pfConst]->pt())/(*inputObjects)[obj].pt();
	    if (deltaR < matchDR_ && relPt < matchDPtRel_) {
	      math::XYZTLorentzVector candP4 = pfConsts[pfConst]->p4();
	      if (debug_) {
		std::cout<<"pT uncorr:"<<candP4.pt()<<std::endl;
		std::cout<<"rawCorrFactor:"<<rawCorrFactor<<std::endl;
		std::cout<<"jecCorrFactor:"<<jecCorrFactor<<std::endl;
		std::cout<<"pT corr:"<<(corrFactor*candP4).pt()<<std::endl;
	      }
	      metpfc += candP4;
	      metjec += (rawCorrFactor*candP4);
	    }
	  }
	}
      }
    }
  
  metpsel->push_back(reco::MET(metsel, reco::MET::Point()));
  metppfc->push_back(reco::MET(metpfc, reco::MET::Point()));
  metpjec->push_back(reco::MET(metjec, reco::MET::Point()));
  
  ev.put(metpsel, "selected");
  ev.put(metppfc, "pfcand");
  ev.put(metpjec, "pfcandjec");
  
}

#endif
