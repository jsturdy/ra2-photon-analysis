import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

### electron cleaned jets (same basic cleaning as against muons)
patJetsPFNoElectron = cleanPatJets.clone()
patJetsPFNoElectron.src                                       = cms.InputTag('newJetsMET')
patJetsPFNoElectron.checkOverlaps.photons.src                 = cms.InputTag('patPhotonsRA2')
patJetsPFNoElectron.checkOverlaps.taus.src                    = cms.InputTag('selectedPatTausPF')
patJetsPFNoElectron.checkOverlaps.electrons.src               = cms.InputTag('patElectrons')
patJetsPFNoElectron.checkOverlaps.electrons.algorithm         = cms.string('byDeltaR')
patJetsPFNoElectron.checkOverlaps.electrons.preselection      = cms.string('')
patJetsPFNoElectron.checkOverlaps.electrons.deltaR            = cms.double(0.2) #changed from 0.1
patJetsPFNoElectron.checkOverlaps.electrons.pairCut           = cms.string('')
patJetsPFNoElectron.checkOverlaps.electrons.requireNoOverlaps = cms.bool(True)
patJetsPFNoElectron.checkOverlaps.tkIsoElectrons.src          = cms.InputTag('patElectrons')
patJetsPFNoElectron.checkOverlaps.muons.src                   = cms.InputTag('patMuonsPF')

patJetsPFNoElectronID = patJetsPFNoElectron.clone()
patJetsPFNoElectronID.checkOverlaps.electrons.src               = cms.InputTag('patElectronsID')
patJetsPFNoElectronID.checkOverlaps.tkIsoElectrons.src          = cms.InputTag('patElectronsID')

patJetsPFNoElectronIDIso = patJetsPFNoElectron.clone()
patJetsPFNoElectronIDIso.checkOverlaps.electrons.src            = cms.InputTag('patElectronsIDIso')
patJetsPFNoElectronIDIso.checkOverlaps.tkIsoElectrons.src       = cms.InputTag('patElectronsIDIso')

####
patJetsPFNoElectronPt30 = patJetsPFNoElectronPt30NoElectron.clone()
patJetsPFNoElectronPt30.checkOverlaps.electrons.src                 = cms.InputTag('patElectrons')
patJetsPFNoElectronPt30.checkOverlaps.tkIsoElectrons.src            = cms.InputTag('patElectrons')

patJetsPFNoElectronPt30ID = patJetsPFNoElectronPt30.clone()
patJetsPFNoElectronPt30ID.checkOverlaps.electrons.src               = cms.InputTag('patElectronsID')
patJetsPFNoElectronPt30.checkOverlaps.tkIsoElectrons.src            = cms.InputTag('patElectronsID')

patJetsPFNoElectronPt30IDIso = patJetsPFNoElectronPt30.clone()
patJetsPFNoElectronPt30IDIso.checkOverlaps.electrons.src            = cms.InputTag('patElectronsIDIso')
patJetsPFNoElectronPt30.checkOverlaps.tkIsoElectrons.src            = cms.InputTag('patElectronsIDIso')
####
patJetsPFNoElectronPt50Eta25     = patJetsPFNoElectronPt30.clone()
patJetsPFNoElectronPt50Eta25.src = cms.InputTag('patJetsPFchsPt50Eta25')

patJetsPFNoElectronPt50Eta25ID = patJetsPFNoElectronPt50Eta25.clone()
patJetsPFNoElectronPt50Eta25ID.checkOverlaps.electrons.src      = cms.InputTag('patElectronsID')
patJetsPFNoElectronPt50Eta25ID.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectronsID')

patJetsPFNoElectronPt50Eta25IDIso = patJetsPFNoElectronPt50Eta25.clone()
patJetsPFNoElectronPt50Eta25IDIso.checkOverlaps.electrons.src      = cms.InputTag('patElectronsIDIso')
patJetsPFNoElectronPt50Eta25IDIso.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectronsIDIso')


electronCleanedPFJetsPF = cms.Sequence(
    patJetsPFNoElectron               *
    patJetsPFNoElectronID             *
    patJetsPFNoElectronIDIso          *
    patJetsPFNoElectronPt30           *
    patJetsPFNoElectronPt30ID         *
    patJetsPFNoElectronPt30IDIso      *
    patJetsPFNoElectronPt50Eta25      *
    patJetsPFNoElectronPt50Eta25ID    *
    patJetsPFNoElectronPt50Eta25IDIso
)

# PFJets - filters

countJetsPFNoElectronPt50Eta25           = countPatJets.clone()
countJetsPFNoElectronPt50Eta25.src       = cms.InputTag('patJetsPFNoElectronPt50Eta25')
countJetsPFNoElectronPt50Eta25.minNumber = cms.uint32(3)

countJetsPFNoElectronPt50Eta25ID     = countJetsPFNoElectronPt50Eta25.clone()
countJetsPFNoElectronPt50Eta25ID.src = cms.InputTag('patJetsPFNoElectronPt50Eta25ID')

countJetsPFNoElectronPt50Eta25IDIso     = countJetsPFNoElectronPt50Eta25.clone()
countJetsPFNoElectronPt50Eta25IDIso.src = cms.InputTag('patJetsPFNoElectronPt50Eta25IDIso')

######
# b-tagged jets
#from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from ZInvisibleBkgds.Photons.specialJetSelector_cff import selectedRA2PatJets

CSVL = 0.244
CSVM = 0.679
CSVT = 0.898 
csvPoint = CSVM

###patJetsPFNoElectron
patCSVJetsPFNoElectron = selectedRA2PatJets.clone()
patCSVJetsPFNoElectron.src = cms.InputTag('patJetsPFNoElectron')
patCSVJetsPFNoElectron.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(csvPoint))

patCSVJetsPFNoElectronPt30Eta24 = selectedRA2PatJets.clone()
patCSVJetsPFNoElectronPt30Eta24.src = cms.InputTag('patCSVJetsPFNoElectron')
patCSVJetsPFNoElectronPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVJetsPFNoElectronPt30Eta50 = selectedRA2PatJets.clone()
patCSVJetsPFNoElectronPt30Eta50.src = cms.InputTag('patCSVJetsPFNoElectron')
patCSVJetsPFNoElectronPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVJetsPFNoElectronPt50Eta25 = selectedRA2PatJets.clone()
patCSVJetsPFNoElectronPt50Eta25.src = cms.InputTag('patCSVJetsPFNoElectron')
patCSVJetsPFNoElectronPt50Eta25.cut = cms.string('pt > 30 && abs(eta) < 2.5')

patCSVMVAJetsPFNoElectron = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoElectron.src = cms.InputTag('patJetsPFNoElectron')
patCSVMVAJetsPFNoElectron.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > %f'%(csvPoint))

patCSVMVAJetsPFNoElectronPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoElectronPt30Eta24.src = cms.InputTag('patCSVMVAJetsPFNoElectron')
patCSVMVAJetsPFNoElectronPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsPFNoElectronPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoElectronPt30Eta50.src = cms.InputTag('patCSVMVAJetsPFNoElectron')
patCSVMVAJetsPFNoElectronPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsPFNoElectronPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoElectronPt50Eta25.src = cms.InputTag('patCSVMVAJetsPFNoElectron')
patCSVMVAJetsPFNoElectronPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

### count the b-jets
###patJetsPFNoElectron
### create the jet collections

zinvBJetsPFNoElectron = cms.Sequence(
      patCSVJetsPFNoElectron
    * patCSVJetsPFNoElectronPt30Eta24
    #* patCSVJetsPFNoElectronPt30Eta50
    #* patCSVJetsPFNoElectronPt50Eta25
    #* patCSVMVAJetsPFNoElectron
    #* patCSVMVAJetsPFNoElectronPt30Eta24
    #* patCSVMVAJetsPFNoElectronPt30Eta50
    #* patCSVMVAJetsPFNoElectronPt50Eta25
)
