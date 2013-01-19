import FWCore.ParameterSet.Config as cms

photonmap = cms.EDProducer("ValueMapProducer",
    photonSrc = cms.InputTag("patPhotonsRA2"),
    src       = cms.InputTag("kt6PFJets","rho"),
)
