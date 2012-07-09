// system include files
#include <memory>
#include <iostream>
#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/Photon.h"

class ValueMapProducer : public edm::EDProducer {
public:
  explicit ValueMapProducer(const edm::ParameterSet&);
  ~ValueMapProducer();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data --------------------------
  edm::InputTag photonSrc_;
  edm::InputTag src_;
};

ValueMapProducer::ValueMapProducer(const edm::ParameterSet& pset)
{
  photonSrc_ = pset.getParameter<edm::InputTag>("photonSrc");
  src_       = pset.getParameter<edm::InputTag>("src");
  produces<edm::ValueMap<float> >().setBranchAlias("PhotonToRhoMap");
}


ValueMapProducer::~ValueMapProducer()
{

}

// ------------ method called to produce the data  ------------
void
ValueMapProducer::produce(edm::Event& ev, const edm::EventSetup& es)
{
  edm::Handle<double> eventRho;
  ev.getByLabel( src_, eventRho );

  edm::Handle<edm::View<pat::Photon> > photons;
  ev.getByLabel(photonSrc_,photons);
  std::vector<float> values;
  values.reserve(photons->size());

  unsigned int photonIdx = 0;
  for(edm::View<pat::Photon>::const_iterator photon = photons->begin();
      photon != photons->end(); ++photon) {
    
    values.push_back(*eventRho);
    ++photonIdx;
  }

  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  filler.insert(photons, values.begin(), values.end());
  filler.fill();

  // put value map into event
  ev.put(out);
}

// ------------ method called once each job just before starting event loop  ------------
void
ValueMapProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
ValueMapProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(ValueMapProducer);
