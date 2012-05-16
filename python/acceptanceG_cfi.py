import FWCore.ParameterSet.Config as cms

acceptanceG = cms.EDProducer("AcceptanceG",
    debug    = cms.bool(True),
    debugString    = cms.string("acceptanceG"),
    genLabel = cms.InputTag("genParticles"),
    genStatus = cms.int32(3),

    jetLabel = cms.InputTag("patJetsAK5PFPt50Eta25"),

    photonLabel  = cms.InputTag("patPhotons"),
    photonMinPt  = cms.double(20.),
    photonMaxEta = cms.double(2.4)
)
