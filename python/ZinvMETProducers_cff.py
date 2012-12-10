import FWCore.ParameterSet.Config as cms
from ZInvisibleBkgds.Photons.specialMETCleaner_cff import *

pfType1MetNoPhotonIDPFIso = specialPhotonCleanedMET.clone()
pfType1MetNoPhotonIDPFIso.inputObjects = cms.InputTag("patPhotonsIDPFIso")
pfType1MetNoPhotonID = specialPhotonCleanedMET.clone()
pfType1MetNoPhotonID.inputObjects = cms.InputTag("patPhotonsID")

pfType1MetNoMuon = specialMuonCleanedMET.clone()
pfType1MetNoMuon.inputObjects = cms.InputTag("patMuonsPF")

pfType1MetNoElectron = specialElectronCleanedMET.clone()
pfType1MetNoElectron.inputObjects = cms.InputTag("patElectronsPF")

photonMETCollections = cms.Sequence(
    pfType1MetNoPhotonID
    *pfType1MetNoPhotonIDPFIso
)

zmumuMETCollections = cms.Sequence(
    pfType1MetNoMuon
)

zelelMETCollections = cms.Sequence(
    pfType1MetNoElectron
)

