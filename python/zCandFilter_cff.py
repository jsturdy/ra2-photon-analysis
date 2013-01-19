
import FWCore.ParameterSet.Config as cms

zCandFilter = cms.EDFilter(
  "ZCandFilter",
  ZCandSource = cms.InputTag("zToMuMu"),
  minPt = cms.double(20.0)
)
