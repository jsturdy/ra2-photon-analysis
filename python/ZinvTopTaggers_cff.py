import FWCore.ParameterSet.Config as cms
from UserCode.TopTagger.topTagger_cfi import *

photonTopTaggerID5Loose = topTagger5Loose.clone()
photonTopTaggerID5Loose.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
photonTopTaggerID5Loose.metSrc = cms.InputTag("pfType1MetNoPhotonID","pfcand")
photonTopTaggerIDPFIso5Loose = topTagger5Loose.clone()
photonTopTaggerIDPFIso5Loose.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
photonTopTaggerIDPFIso5Loose.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","pfcand")
photonTopTaggerID5M = topTagger5M.clone()
photonTopTaggerID5M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
photonTopTaggerID5M.metSrc = cms.InputTag("pfType1MetNoPhotonID","pfcand")
photonTopTaggerIDPFIso5M = topTagger5M.clone()
photonTopTaggerIDPFIso5M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
photonTopTaggerIDPFIso5M.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","pfcand")

diMuonTopTagger5Loose = topTagger5Loose.clone()
diMuonTopTagger5Loose.jetSrc = cms.InputTag("patJetsPFNoMuonPt30")
diMuonTopTagger5Loose.metSrc = cms.InputTag("pfType1MetNoMuon","pfcand")
diMuonTopTagger5M = topTagger5M.clone()
diMuonTopTagger5M.jetSrc = cms.InputTag("patJetsPFNoMuonPt30")
diMuonTopTagger5M.metSrc = cms.InputTag("pfType1MetNoMuon","pfcand")

diElectronTopTagger5Loose = topTagger5Loose.clone()
diElectronTopTagger5Loose.jetSrc = cms.InputTag("patJetsPFNoElectronPt30")
diElectronTopTagger5Loose.metSrc = cms.InputTag("pfType1MetNoElectron","pfcand")
diElectronTopTagger5M = topTagger5M.clone()
diElectronTopTagger5M.jetSrc = cms.InputTag("patJetsPFNoElectronPt30")
diElectronTopTagger5M.metSrc = cms.InputTag("pfType1MetNoElectron","pfcand")

zInvTopTagger5Loose = topTagger5Loose.clone()
zInvTopTagger5Loose.jetSrc = cms.InputTag("patJetsPFchsPt30")
zInvTopTagger5Loose.metSrc = cms.InputTag("patMETsPF")
zInvTopTagger5M = topTagger5M.clone()
zInvTopTagger5M.jetSrc = cms.InputTag("patJetsPFchsPt30")
zInvTopTagger5M.metSrc = cms.InputTag("patMETsPF")


photonTopTaggers = cms.Sequence(
    photonTopTaggerID5Loose
    *photonTopTaggerID5M
    #*photonTopTaggerIDPFIso5Loose
    #*photonTopTaggerIDPFIso5M
)

zinvTopTaggers = cms.Sequence(
    zInvTopTagger5Loose
    *zInvTopTagger5M
)

zmumuTopTaggers = cms.Sequence(
    diMuonTopTagger5Loose
    *diMuonTopTagger5M
)

zelelTopTaggers = cms.Sequence(
    diElectronTopTagger5Loose
    *diElectronTopTagger5M
)

