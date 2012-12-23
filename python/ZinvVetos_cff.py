import FWCore.ParameterSet.Config as cms
from SandBox.Skims.basicJetSelector_cfi import selectedPatJets
patJetsForIndirectTauVetoPhotonsID = selectedPatJets.clone()
patJetsForIndirectTauVetoPhotonsID.src = cms.InputTag("patJetsPFNoPhotonIDSpecial")
patJetsForIndirectTauVetoPhotonsID.cut = cms.string('pt > 15 && abs(eta) < 2.4 && bDiscriminator("combinedSecondaryVertexBJetTags") <= 0.898')

patJetsForIndirectTauVetoPhotonsIDPFIso = patJetsForIndirectTauVetoPhotonsID.clone()
patJetsForIndirectTauVetoPhotonsIDPFIso.src = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecial")

#######
patJetsForIndirectTauVetoDiMuons = patJetsForIndirectTauVetoPhotonsID.clone()
patJetsForIndirectTauVetoDiMuons.src = cms.InputTag("patJetsPFNoMuon")

#######
patJetsForIndirectTauVetoDiElectrons = patJetsForIndirectTauVetoPhotonsID.clone()
patJetsForIndirectTauVetoDiElectrons.src = cms.InputTag("patJetsPFNoElectron")

#######
patJetsForIndirectTauVetoZInv = patJetsForIndirectTauVetoPhotonsID.clone()
patJetsForIndirectTauVetoZInv.src = cms.InputTag("patJetsPF")


#================ configure filters =======================#

from SandBox.HadSTopSkims.muonVetoFilter_cfi import muonVetoSTop
sTopPFMuonVeto       = muonVetoSTop.clone()
sTopPFMuonVetoDiMuon = muonVetoSTop.clone(MaxMuons = cms.int32(2))

from SandBox.HadSTopSkims.electronVetoFilter_cfi import electronVetoSTop
sTopPFElectronVeto           = electronVetoSTop.clone()
sTopPFElectronVetoDiElectron = electronVetoSTop.clone(MaxElectrons = cms.int32(2))

from SandBox.HadSTopSkims.trackIsolationMaker_cfi import trackIsolationMaker
sTopTrkIsolationMaker = trackIsolationMaker.clone()

from SandBox.HadSTopSkims.isolatedTrackVeto_cfi import isolatedTrackVeto
sTopIsoTrkVeto          = isolatedTrackVeto.clone(
    CandidatesCharge     = cms.InputTag("sTopTrkIsolationMaker","pfcandschg"),
    CandidatesTrkIso     = cms.InputTag("sTopTrkIsolationMaker","pfcandstrkiso"),
    CandidatesPt         = cms.InputTag("sTopTrkIsolationMaker","pfcandspt"),
    CandidatesDZ         = cms.InputTag("sTopTrkIsolationMaker","pfcandsdzpv"),
)
sTopIsoTrkVetoDiLeptons = sTopIsoTrkVeto.clone(MaxChargedCandidates = cms.int32(2))

from SandBox.HadSTopSkims.indirectTauVeto_cfi import indirectTauVeto
sTopTauVetoDiMuon = indirectTauVeto.clone(
    JetSource     = cms.InputTag("patJetsForIndirectTauVetoDiMuons"),
    METSource     = cms.InputTag("pfType1MetNoMuon","pfcand"),
)
sTopTauVetoDiElectron = indirectTauVeto.clone(
    JetSource     = cms.InputTag("patJetsForIndirectTauVetoDiElectrons"),
    METSource     = cms.InputTag("pfType1MetNoElectron","pfcand"),
)
sTopTauVetoPhotonID = indirectTauVeto.clone(
    JetSource     = cms.InputTag("patJetsForIndirectTauVetoPhotonsID"),
    METSource     = cms.InputTag("pfType1MetNoPhotonID","pfcand"),
)
sTopTauVetoPhotonIDPFIso = indirectTauVeto.clone(
    JetSource     = cms.InputTag("patJetsForIndirectTauVetoPhotonsIDPFIso"),
    METSource     = cms.InputTag("pfType1MetNoPhotonIDPFIso","pfcand"),
)
sTopTauVetoZInv = indirectTauVeto.clone(
    JetSource   = cms.InputTag("patJetsForIndirectTauVetoZInv"),
    METSource   = cms.InputTag("patMETsPF"),
)




photonVetos = cms.Sequence(
    patJetsForIndirectTauVetoPhotonsID
#    *patJetsForIndirectTauVetoPhotonsIDPFIso
    *sTopPFMuonVeto
    *sTopPFElectronVeto
    *sTopTrkIsolationMaker
    *sTopIsoTrkVeto
    *sTopTauVetoPhotonID
#    *sTopTauVetoPhotonIDPFIso
)

zinvVetos = cms.Sequence(
    patJetsForIndirectTauVetoZInv
    *sTopPFMuonVeto
    *sTopPFElectronVeto
    *sTopTrkIsolationMaker
    *sTopIsoTrkVeto
    *sTopTauVetoZInv
)

zmumuVetos = cms.Sequence(
    patJetsForIndirectTauVetoDiMuons
    *sTopPFMuonVetoDiMuon
    *sTopPFElectronVeto
    *sTopTrkIsolationMaker
    *sTopIsoTrkVetoDiLeptons
    *sTopTauVetoDiMuon
)

zelelVetos = cms.Sequence(
    patJetsForIndirectTauVetoDiElectrons
    *sTopPFMuonVeto
    *sTopPFElectronVetoDiElectron
    *sTopTrkIsolationMaker
    *sTopIsoTrkVetoDiLeptons
    *sTopTauVetoDiElectron
)

