import FWCore.ParameterSet.Config as cms

photonmap = cms.EDProducer("ValueMapProducer",
    photonSrc = cms.InputTag("patPhotonsAlt"),
    src       = cms.InputTag("kt6PFJetsForIsolation","rho"),
)
