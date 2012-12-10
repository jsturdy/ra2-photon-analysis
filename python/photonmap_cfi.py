import FWCore.ParameterSet.Config as cms

photonmap = cms.EDProducer("ValueMapProducer",
    photonSrc = cms.InputTag("patPhotons"),
    src       = cms.InputTag("kt6PFJets","rho"),
)
