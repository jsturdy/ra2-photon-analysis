import FWCore.ParameterSet.Config as cms

genbosoncompG = cms.EDProducer("GenBosonComparison",
    debug          = cms.bool(False),
    debugString    = cms.string("genbosoncompG"),
    genLabel       = cms.InputTag("genParticles"),
    genStatus      = cms.int32(3),
    genPDGId       = cms.int32(22),

    doPUReweight = cms.bool(False),
    puWeight     = cms.InputTag("puWeight"),

    jetHTLabel  = cms.InputTag("patJetsAK5PFPt50Eta25NoPhotonIDIso"),
    jetMHTLabel = cms.InputTag("patJetsAK5PFPt30NoPhotonIDIso"),
    htLabel     = cms.InputTag("htPFNoPhotIDIso"),
    mhtLabel    = cms.InputTag("mhtPFNoPhotIDIso"),

    photonEBMaxEta = cms.double(1.4442),
    photonEEMinEta = cms.double(1.566),
    photonEEMaxEta = cms.double(2.5)
)
