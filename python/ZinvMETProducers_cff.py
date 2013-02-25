import FWCore.ParameterSet.Config as cms
from ZInvisibleBkgds.Photons.specialMETCleaner_cff import *

pfType1MetNoPhotonIDPFIso = specialPhotonCleanedMET.clone()
pfType1MetNoPhotonIDPFIso.inputObjects = cms.InputTag("patPhotonsIDPFIso")
pfType1MetNoPhotonFitTemplate = specialPhotonCleanedMET.clone()
pfType1MetNoPhotonFitTemplate.inputObjects = cms.InputTag("patFitTemplatePhotons")
pfType1MetNoPhotonJetFake = specialPhotonCleanedMET.clone()
pfType1MetNoPhotonJetFake.inputObjects = cms.InputTag("patJetFakePhotons")
pfType1MetNoPhotonID = specialPhotonCleanedMET.clone()
pfType1MetNoPhotonID.inputObjects = cms.InputTag("patPhotonsID")

pfType1MetNoMuon = specialMuonCleanedMET.clone()
pfType1MetNoMuon.inputObjects = cms.InputTag("patMuonsPFIDIso")

pfType1MetNoElectron = specialElectronCleanedMET.clone()
pfType1MetNoElectron.inputObjects = cms.InputTag("patElectronsIDIso")

photonMETCollections = cms.Sequence(
    pfType1MetNoPhotonID
  * pfType1MetNoPhotonIDPFIso
)
photonTemplateMETCollections = cms.Sequence(
    pfType1MetNoPhotonFitTemplate
  * pfType1MetNoPhotonJetFake
)

zmumuMETCollections = cms.Sequence(
    pfType1MetNoMuon
)

zelelMETCollections = cms.Sequence(
    pfType1MetNoElectron
)

