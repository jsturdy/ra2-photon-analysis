import FWCore.ParameterSet.Config as cms

specialMuonCollection = cms.EDProducer("SpecialPATMuonCollection",
    debug          = cms.bool(False),
    debugString    = cms.string("specialMuonCollection"),
    candidateLabel = cms.InputTag("zMuMuCands"),
    mZ             = cms.double(91.19)
)
