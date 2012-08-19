
from PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.muonCountFilter_cfi import *

muonIDCutTight  = cms.string('isGlobalMuon && isPFMuon && globalTrack.normalizedChi2 < 10 &&'
                             'globalTrack.hitPattern.numberOfValidMuonHits > 0 && '
                             'numberOfMatchedStations > 1 && dB < 0.2 && '
                             'fabs(innerTrack.dz(vertex.position)) < 0.5 && '
                             'innerTrack.hitPattern.numberOfValidPixelHits > 0 && '
                             'track.hitPattern.trackerLayersWithMeasurement > 5'
                            )

muonISOCutTight  = cms.string('(pfIsolationR04.sumChargedHadronPt+max(0,pfIsolationR04.sumNeutralHadronEt+pfIsolationR04.sumPhotonEt)-0.5*pfIsolationR04.sumPUPt)/pt'
                             )



patMuonsID = cms.EDFilter(
  "PATMuonSelector",
   src = cms.InputTag('patMuons'),
   cut = muonIDCutTight,
   filter = cms.bool(False),
)
patMuonsIDPFIso = patMuonsID.clone(src = cms.InputTag('patMuonsID'))
patMuonsIDPFIso.cut = muonISOCutTight

patMuonsIDIso = patMuonsID.clone(
   src = cms.InputTag('patMuonsID'),
   cut = muoncombiso03cut,
)

patMuonRefsID = cms.EDFilter(
  "PATMuonRefSelector",
   src = cms.InputTag('patMuons'),
   cut = muonIDCutTight,
   filter = cms.bool(False)
)

zinvMuons = cms.Sequence(
    patMuonsID
  * patMuonsIDPFIso
)

countMuonsID = countPatMuons.clone()
countMuonsID.src = cms.InputTag("patMuonsID")
countMuonsID.minNumber = cms.uint32(1)

countMuonsIDPFIso = countPatMuons.clone()
countMuonsIDPFIso.src = cms.InputTag("patMuonsIDPFIso")
countMuonsIDPFIso.minNumber = cms.uint32(1)

countMuonRefsID = countPatMuons.clone()
countMuonRefsID.src = cms.InputTag("patMuonsRefsID")
countMuonRefsID.minNumber = cms.uint32(1)
