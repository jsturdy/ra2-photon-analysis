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
patJetsForIndirectTauVetoZInv.src = cms.InputTag("newJetsMET")


#================ configure filters =======================#

#from SandBox.HadSTopSkims.muonVetoFilter_cfi import muonVetoSTop
#sTopPFMuonVeto       = muonVetoSTop.clone()
#sTopPFMuonVetoDiMuon = muonVetoSTop.clone(MaxMuons = cms.int32(2))
from SandBox.Stop.sakLooseMuonSelector_cfi import sakLooseMuonSelector
sTopPFMuonVeto       = sakLooseMuonSelector.clone()
sTopPFMuonVetoDiMuon = sakLooseMuonSelector.clone(MaxMuons = cms.int32(2))
from SandBox.Skims.RA2Leptons_cff import countMuonsIDIso,countPFMuonsIDIso
countRA2MuonsIDIso          = countMuonsIDIso.clone(minNumber = cms.uint32(0),
                                                    maxNumber = cms.uint32(0))
countRA2MuonsPFIDIso        = countPFMuonsIDIso.clone(minNumber = cms.uint32(0),
                                                      maxNumber = cms.uint32(0))
countRA2MuonsIDIsoDiMuons   = countMuonsIDIso.clone(minNumber = cms.uint32(2),maxNumber = cms.uint32(2))
countRA2MuonsPFIDIsoDiMuons = countPFMuonsIDIso.clone(minNumber = cms.uint32(2),maxNumber = cms.uint32(2))

##from SandBox.HadSTopSkims.electronVetoFilter_cfi import electronVetoSTop
##sTopPFElectronVeto           = electronVetoSTop.clone()
##sTopPFElectronVetoDiElectron = electronVetoSTop.clone(MaxElectrons = cms.int32(2))
from SandBox.Stop.sakLooseElectronSelector_cfi import sakLooseElectronSelector
sTopPFElectronVeto           = sakLooseElectronSelector.clone()
sTopPFElectronVetoDiElectron = sakLooseElectronSelector.clone(MaxElectrons = cms.int32(2))
from SandBox.Skims.RA2Leptons_cff import countElectronsIDIso
countRA2ElectronsIDIso            = countElectronsIDIso.clone(minNumber = cms.uint32(0),
                                                              maxNumber = cms.uint32(0))
countRA2ElectronsIDIsoDiElectrons = countElectronsIDIso.clone(minNumber = cms.uint32(2),maxNumber = cms.uint32(2))

#from SandBox.HadSTopSkims.trackIsolationMaker_cfi import trackIsolationMaker
#sTopTrkIsolationMaker = trackIsolationMaker.clone()
from SandBox.Stop.trackIsolationMaker_cfi import trackIsolationFilter
sTopTrkIsolationMaker = trackIsolationFilter.clone(doTrkIsoVeto = cms.bool(False))
sTopTrkIsolationMakerDiLeptons = sTopTrkIsolationMaker.clone(MaxChargedCandidates = cms.uint32(2))

##from SandBox.HadSTopSkims.isolatedTrackVeto_cfi import isolatedTrackVeto
##sTopIsoTrkVeto          = isolatedTrackVeto.clone(
##    CandidatesCharge     = cms.InputTag("sTopTrkIsolationMaker","pfcandschg"),
##    CandidatesTrkIso     = cms.InputTag("sTopTrkIsolationMaker","pfcandstrkiso"),
##    CandidatesPt         = cms.InputTag("sTopTrkIsolationMaker","pfcandspt"),
##    CandidatesDZ         = cms.InputTag("sTopTrkIsolationMaker","pfcandsdzpv"),
##)
##sTopIsoTrkVetoDiLeptons = sTopIsoTrkVeto.clone(MaxChargedCandidates = cms.int32(2))

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
    METSource   = cms.InputTag("newMETwPhiCorr"),
)




photonVetos = cms.Sequence(
    #patJetsForIndirectTauVetoPhotonsID
    #*patJetsForIndirectTauVetoPhotonsIDPFIso
    sTopPFMuonVeto
    *sTopPFElectronVeto
    *sTopTrkIsolationMaker
    #*sTopIsoTrkVeto
    #*sTopTauVetoPhotonID
    #*sTopTauVetoPhotonIDPFIso
)

zinvVetos = cms.Sequence(
    #patJetsForIndirectTauVetoZInv
    sTopPFMuonVeto
    *sTopPFElectronVeto
    *sTopTrkIsolationMaker
    #*sTopIsoTrkVeto
    #*sTopTauVetoZInv
)

zmumuVetos = cms.Sequence(
    #patJetsForIndirectTauVetoDiMuons
    sTopPFMuonVetoDiMuon
    *sTopPFElectronVeto
    *sTopTrkIsolationMakerDiLeptons
#    *sTopIsoTrkVetoDiLeptons
    #*sTopTauVetoDiMuon
)

zelelVetos = cms.Sequence(
    #patJetsForIndirectTauVetoDiElectrons
    sTopPFMuonVeto
    *sTopPFElectronVetoDiElectron
    *sTopTrkIsolationMakerDiLeptons
#    *sTopIsoTrkVetoDiLeptons
    #*sTopTauVetoDiElectron
)

