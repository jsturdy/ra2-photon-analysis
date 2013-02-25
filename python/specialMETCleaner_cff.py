import FWCore.ParameterSet.Config as cms

specialMuonCleanedMET = cms.EDProducer("SpecialPATMuonCleanedMETProducer",
    inputObjects    = cms.InputTag("patMuonsPFIDIso"),
    inputMET        = cms.InputTag("newMETwPhiCorr"),
    inputJets       = cms.InputTag("newJetsMET"),
    inputPFCands    = cms.InputTag("particleFlow"),
    objectsToRemove = cms.int32(2),
    particleId      = cms.int32(3),
    matchDR         = cms.double(0.3),
    matchDPtRel     = cms.double(0.3),
    debug           = cms.bool(False),
)

specialElectronCleanedMET = cms.EDProducer("SpecialPATElectronCleanedMETProducer",
    inputObjects    = cms.InputTag("patElectronsIDIso"),
    inputMET        = cms.InputTag("newMETwPhiCorr"),
    inputJets       = cms.InputTag("newJetsMET"),
    inputPFCands    = cms.InputTag("particleFlow"),
    objectsToRemove = cms.int32(2),
    particleId      = cms.int32(2),
    matchDR         = cms.double(0.3),
    matchDPtRel     = cms.double(0.3),
    debug           = cms.bool(False),
)

specialPhotonCleanedMET = cms.EDProducer("SpecialPATPhotonCleanedMETProducer",
    inputObjects    = cms.InputTag("patPhotonsRA2"),
    inputMET        = cms.InputTag("newMETwPhiCorr"),
    inputJets       = cms.InputTag("newJetsMET"),
    inputPFCands    = cms.InputTag("particleFlow"),
    objectsToRemove = cms.int32(1),
    particleId      = cms.int32(4),
    matchDR         = cms.double(0.3),
    matchDPtRel     = cms.double(0.3),
    debug           = cms.bool(False),
)
