
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include <DataFormats/Candidate/interface/Candidate.h>


class ZCandFilter : public edm::EDFilter {

  public:

    explicit ZCandFilter(const edm::ParameterSet & iConfig);
    ~ZCandFilter() {}

  private:

    virtual bool filter(edm::Event & iEvent, const edm::EventSetup & iSetup);
    
    edm::InputTag candSrc_;
    double minPt_;
};


ZCandFilter::ZCandFilter(const edm::ParameterSet & iConfig) {
  candSrc_ = iConfig.getParameter<edm::InputTag>("ZCandSource");
  minPt_   = iConfig.getParameter<double>("minPt");
}


bool ZCandFilter::filter(edm::Event & iEvent, const edm::EventSetup & iSetup) {
  edm::Handle<std::vector<reco::CompositeCandidate> > cands;
  //edm::Handle<edm::View<reco::CompositeCandidate> > cands;
  iEvent.getByLabel(candSrc_, cands);
  bool result = false;
  if (cands->size() > 0)
    if ((*cands)[0].pt() > minPt_)
      result = true;
  return result;
}


#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(ZCandFilter);
