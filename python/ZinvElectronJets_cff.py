from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

### electron cleaned jets (same basic cleaning as against muons)
patJetsAK5PFNoElectron = cleanPatJets.clone()
patJetsAK5PFNoElectron.src                                       = cms.InputTag('patJetsAK5PF')
patJetsAK5PFNoElectron.checkOverlaps.photons.src                 = cms.InputTag('patPhotonsAlt')
patJetsAK5PFNoElectron.checkOverlaps.taus.src                    = cms.InputTag('selectedPatTausPF')
patJetsAK5PFNoElectron.checkOverlaps.electrons.src               = cms.InputTag('patElectrons')
patJetsAK5PFNoElectron.checkOverlaps.electrons.algorithm         = cms.string('byDeltaR')
patJetsAK5PFNoElectron.checkOverlaps.electrons.preselection      = cms.string('')
patJetsAK5PFNoElectron.checkOverlaps.electrons.deltaR            = cms.double(0.1)
patJetsAK5PFNoElectron.checkOverlaps.electrons.pairCut           = cms.string('')
patJetsAK5PFNoElectron.checkOverlaps.electrons.requireNoOverlaps = cms.bool(True)
patJetsAK5PFNoElectron.checkOverlaps.tkIsoElectrons.src          = cms.InputTag('patElectrons')
patJetsAK5PFNoElectron.checkOverlaps.muons.src                   = cms.InputTag('patMuons')

patJetsAK5PFNoElectronID = patJetsAK5PFNoElectron.clone()
patJetsAK5PFNoElectronID.checkOverlaps.electrons.src               = cms.InputTag('patElectronsID')
patJetsAK5PFNoElectronID.checkOverlaps.tkIsoElectrons.src          = cms.InputTag('patElectronsID')

patJetsAK5PFNoElectronIDIso = patJetsAK5PFNoElectron.clone()
patJetsAK5PFNoElectronIDIso.checkOverlaps.electrons.src            = cms.InputTag('patElectronsIDIso')
patJetsAK5PFNoElectronIDIso.checkOverlaps.tkIsoElectrons.src       = cms.InputTag('patElectronsIDIso')

patJetsPFNoElectron = patJetsPFNoPFElectron.clone()
patJetsPFNoElectron.src                 = cms.InputTag('patJetsPF')

patJetsPFNoElectronID = patJetsPFNoElectron.clone()
patJetsPFNoElectronID.checkOverlaps.electrons.src               = cms.InputTag('patElectronsID')
patJetsPFNoElectronID.checkOverlaps.tkIsoElectrons.src          = cms.InputTag('patElectronsID')

patJetsPFNoElectronIDIso = patJetsPFNoElectron.clone()
patJetsPFNoElectronIDIso.checkOverlaps.electrons.src            = cms.InputTag('patElectronsIDIso')
patJetsPFNoElectronIDIso.checkOverlaps.tkIsoElectrons.src       = cms.InputTag('patElectronsIDIso')

####
patJetsAK5PFNoElectronPt30ID     = patJetsAK5PFNoElectron.clone()
patJetsAK5PFNoElectronPt30ID.src = cms.InputTag('patJetsAK5PFPt30')

patJetsAK5PFNoElectronPt30IDIso = patJetsAK5PFNoElectronPt30.clone()
patJetsAK5PFNoElectronPt30IDIso.checkOverlaps.electrons.src            = cms.InputTag('patElectronsIDIso')
patJetsAK5PFNoElectronPt30.checkOverlaps.tkIsoElectrons.src            = cms.InputTag('patElectronsIDIso')

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
patJetsAK5PFNoElectronPt50Eta25     = patJetsAK5PFNoElectronPt30.clone()
patJetsAK5PFNoElectronPt50Eta25.src = cms.InputTag('patJetsAK5PFPt50Eta25')

patJetsAK5PFNoElectronPt50Eta25ID = patJetsAK5PFNoElectronPt50Eta25.clone()
patJetsAK5PFNoElectronPt50Eta25ID.checkOverlaps.electrons.src      = cms.InputTag('patElectronsID')
patJetsAK5PFNoElectronPt50Eta25ID.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectronsID')

patJetsAK5PFNoElectronPt50Eta25IDIso = patJetsAK5PFNoElectronPt50Eta25.clone()
patJetsAK5PFNoElectronPt50Eta25IDIso.checkOverlaps.electrons.src      = cms.InputTag('patElectronsIDIso')
patJetsAK5PFNoElectronPt50Eta25IDIso.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectronsIDIso')

patJetsPFNoElectronPt50Eta25     = patJetsPFNoElectronPt30.clone()
patJetsPFNoElectronPt50Eta25.src = cms.InputTag('patJetsPFPt50Eta25')

patJetsPFNoElectronPt50Eta25ID = patJetsPFNoElectronPt50Eta25.clone()
patJetsPFNoElectronPt50Eta25ID.checkOverlaps.electrons.src      = cms.InputTag('patElectronsID')
patJetsPFNoElectronPt50Eta25ID.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectronsID')

patJetsPFNoElectronPt50Eta25IDIso = patJetsPFNoElectronPt50Eta25.clone()
patJetsPFNoElectronPt50Eta25IDIso.checkOverlaps.electrons.src      = cms.InputTag('patElectronsIDIso')
patJetsPFNoElectronPt50Eta25IDIso.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectronsIDIso')


electronCleanedPFJetsAK5PF = cms.Sequence(
    patJetsAK5PFNoElectron               *
    patJetsAK5PFNoElectronID             *
    patJetsAK5PFNoElectronIDIso          *
    patJetsAK5PFNoElectronPt30           *
    patJetsAK5PFNoElectronPt30ID         *
    patJetsAK5PFNoElectronPt30IDIso      *
    patJetsAK5PFNoElectronPt50Eta25      *
    patJetsAK5PFNoElectronPt50Eta25ID    *
    patJetsAK5PFNoElectronPt50Eta25IDIso
)

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

countJetsAK5PFNoElectronPt50Eta25           = countPatJets.clone()
countJetsAK5PFNoElectronPt50Eta25.src       = cms.InputTag('patJetsAK5PFNoElectronPt50Eta25')
countJetsAK5PFNoElectronPt50Eta25.minNumber = cms.uint32(3)

countJetsAK5PFNoElectronPt50Eta25ID     = countJetsAK5PFNoElectronPt50Eta25.clone()
countJetsAK5PFNoElectronPt50Eta25ID.src = cms.InputTag('patJetsAK5PFNoElectronPt50Eta25ID')

countJetsAK5PFNoElectronPt50Eta25IDIso     = countJetsAK5PFNoElectronPt50Eta25.clone()
countJetsAK5PFNoElectronPt50Eta25IDIso.src = cms.InputTag('patJetsAK5PFNoElectronPt50Eta25IDIso')

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

patSSVHEMBJetsAK5PFNoElectron = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PFNoElectron.src = cms.InputTag('patJetsAK5PFNoElectron')
patSSVHEMBJetsAK5PFNoElectron.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFNoElectron = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PFNoElectron.src = cms.InputTag('patJetsAK5PFNoElectron')
patSSVHPTBJetsAK5PFNoElectron.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsAK5PFNoElectronPt30 = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PFNoElectronPt30.src = cms.InputTag('patJetsAK5PFNoElectronPt30')
patSSVHEMBJetsAK5PFNoElectronPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFNoElectronPt30 = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PFNoElectronPt30.src = cms.InputTag('patJetsAK5PFNoElectronPt30')
patSSVHPTBJetsAK5PFNoElectronPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsAK5PFNoElectronPt50Eta25 = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PFNoElectronPt50Eta25.src = cms.InputTag('patJetsAK5PFNoElectronPt50Eta25')
patSSVHEMBJetsAK5PFNoElectronPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFNoElectronPt50Eta25 = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PFNoElectronPt50Eta25.src = cms.InputTag('patJetsAK5PFNoElectronPt50Eta25')
patSSVHPTBJetsAK5PFNoElectronPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patCSVJetsAK5PFNoElectron = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoElectron.src = cms.InputTag('patJetsAK5PFNoElectron')
patCSVJetsAK5PFNoElectron.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898')

patCSVJetsAK5PFNoElectronPt30Eta24 = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoElectronPt30Eta24.src = cms.InputTag('patCSVJetsAK5PFNoElectron')
patCSVJetsAK5PFNoElectronPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVJetsAK5PFNoElectronPt30Eta50 = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoElectronPt30Eta50.src = cms.InputTag('patCSVJetsAK5PFNoElectron')
patCSVJetsAK5PFNoElectronPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVJetsAK5PFNoElectronPt50Eta25 = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoElectronPt50Eta25.src = cms.InputTag('patCSVJetsAK5PFNoElectron')
patCSVJetsAK5PFNoElectronPt50Eta25.cut = cms.string('pt > 30 && abs(eta) < 2.5')

patCSVMVAJetsAK5PFNoElectron = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoElectron.src = cms.InputTag('patJetsAK5PFNoElectron')
patCSVMVAJetsAK5PFNoElectron.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > 0.898')

patCSVMVAJetsAK5PFNoElectronPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoElectronPt30Eta24.src = cms.InputTag('patCSVMVAJetsAK5PFNoElectron')
patCSVMVAJetsAK5PFNoElectronPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsAK5PFNoElectronPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoElectronPt30Eta50.src = cms.InputTag('patCSVMVAJetsAK5PFNoElectron')
patCSVMVAJetsAK5PFNoElectronPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsAK5PFNoElectronPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoElectronPt50Eta25.src = cms.InputTag('patCSVMVAJetsAK5PFNoElectron')
patCSVMVAJetsAK5PFNoElectronPt50Eta25.cut = cms.string('pt > 30 && abs(eta) < 2.5')

###patJetsPFNoElectron
patSSVHEMBJetsPFNoElectron = selectedRA2PatJets.clone()
patSSVHEMBJetsPFNoElectron.src = cms.InputTag('patJetsPFNoElectron')
patSSVHEMBJetsPFNoElectron.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPFNoElectron = selectedRA2PatJets.clone()
patSSVHPTBJetsPFNoElectron.src = cms.InputTag('patJetsPFNoElectron')
patSSVHPTBJetsPFNoElectron.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsPFNoElectronPt30 = selectedRA2PatJets.clone()
patSSVHEMBJetsPFNoElectronPt30.src = cms.InputTag('patJetsPFNoElectronPt30')
patSSVHEMBJetsPFNoElectronPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPFNoElectronPt30 = selectedRA2PatJets.clone()
patSSVHPTBJetsPFNoElectronPt30.src = cms.InputTag('patJetsPFNoElectronPt30')
patSSVHPTBJetsPFNoElectronPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsPFNoElectronPt50Eta25 = selectedRA2PatJets.clone()
patSSVHEMBJetsPFNoElectronPt50Eta25.src = cms.InputTag('patJetsPFNoElectronPt50Eta25')
patSSVHEMBJetsPFNoElectronPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPFNoElectronPt50Eta25 = selectedRA2PatJets.clone()
patSSVHPTBJetsPFNoElectronPt50Eta25.src = cms.InputTag('patJetsPFNoElectronPt50Eta25')
patSSVHPTBJetsPFNoElectronPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patCSVJetsPFNoElectron = selectedRA2PatJets.clone()
patCSVJetsPFNoElectron.src = cms.InputTag('patJetsPFNoElectron')
patCSVJetsPFNoElectron.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898')

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
patCSVMVAJetsPFNoElectron.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > 0.898')

patCSVMVAJetsPFNoElectronPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoElectronPt30Eta24.src = cms.InputTag('patCSVMVAJetsPFNoElectron')
patCSVMVAJetsPFNoElectronPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsPFNoElectronPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoElectronPt30Eta50.src = cms.InputTag('patCSVMVAJetsPFNoElectron')
patCSVMVAJetsPFNoElectronPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsPFNoElectronPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoElectronPt50Eta25.src = cms.InputTag('patCSVMVAJetsPFNoElectron')
patCSVMVAJetsPFNoElectronPt50Eta25.cut = cms.string('pt > 30 && abs(eta) < 2.5')

### count the b-jets
countSSVHEMBJetsAK5PFNoElectron = countPatJets.clone()
countSSVHEMBJetsAK5PFNoElectron.src = cms.InputTag('patSSVHEMBJetsAK5PFNoElectron')
countSSVHEMBJetsAK5PFNoElectron.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PFNoElectron = countPatJets.clone()
countSSVHPTBJetsAK5PFNoElectron.src = cms.InputTag('patSSVHPTBJetsAK5PFNoElectron')
countSSVHPTBJetsAK5PFNoElectron.minNumber = cms.uint32(1)

countSSVHEMBJetsAK5PFNoElectronPt30 = countPatJets.clone()
countSSVHEMBJetsAK5PFNoElectronPt30.src = cms.InputTag('patSSVHEMBJetsAK5PFNoElectronPt30')
countSSVHEMBJetsAK5PFNoElectronPt30.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PFNoElectronPt30 = countPatJets.clone()
countSSVHPTBJetsAK5PFNoElectronPt30.src = cms.InputTag('patSSVHPTBJetsAK5PFNoElectronPt30')
countSSVHPTBJetsAK5PFNoElectronPt30.minNumber = cms.uint32(1)

countSSVHEMBJetsAK5PFNoElectronPt50Eta25 = countPatJets.clone()
countSSVHEMBJetsAK5PFNoElectronPt50Eta25.src = cms.InputTag('patSSVHEMBJetsAK5PFNoElectronPt50Eta25')
countSSVHEMBJetsAK5PFNoElectronPt50Eta25.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PFNoElectronPt50Eta25 = countPatJets.clone()
countSSVHPTBJetsAK5PFNoElectronPt50Eta25.src = cms.InputTag('patSSVHPTBJetsAK5PFNoElectronPt50Eta25')
countSSVHPTBJetsAK5PFNoElectronPt50Eta25.minNumber = cms.uint32(1)

###patJetsPFNoElectron
countSSVHEMBJetsPFNoElectron = countPatJets.clone()
countSSVHEMBJetsPFNoElectron.src = cms.InputTag('patSSVHEMBJetsPFNoElectron')
countSSVHEMBJetsPFNoElectron.minNumber = cms.uint32(1)

countSSVHPTBJetsPFNoElectron = countPatJets.clone()
countSSVHPTBJetsPFNoElectron.src = cms.InputTag('patSSVHPTBJetsPFNoElectron')
countSSVHPTBJetsPFNoElectron.minNumber = cms.uint32(1)

countSSVHEMBJetsPFNoElectronPt30 = countPatJets.clone()
countSSVHEMBJetsPFNoElectronPt30.src = cms.InputTag('patSSVHEMBJetsPFNoElectronPt30')
countSSVHEMBJetsPFNoElectronPt30.minNumber = cms.uint32(1)

countSSVHPTBJetsPFNoElectronPt30 = countPatJets.clone()
countSSVHPTBJetsPFNoElectronPt30.src = cms.InputTag('patSSVHPTBJetsPFNoElectronPt30')
countSSVHPTBJetsPFNoElectronPt30.minNumber = cms.uint32(1)

countSSVHEMBJetsPFNoElectronPt50Eta25 = countPatJets.clone()
countSSVHEMBJetsPFNoElectronPt50Eta25.src = cms.InputTag('patSSVHEMBJetsPFNoElectronPt50Eta25')
countSSVHEMBJetsPFNoElectronPt50Eta25.minNumber = cms.uint32(1)

countSSVHPTBJetsPFNoElectronPt50Eta25 = countPatJets.clone()
countSSVHPTBJetsPFNoElectronPt50Eta25.src = cms.InputTag('patSSVHPTBJetsPFNoElectronPt50Eta25')
countSSVHPTBJetsPFNoElectronPt50Eta25.minNumber = cms.uint32(1)


zinvBVeto = cms.Sequence(
    ~countSSVHEMBJetsAK5PFNoElectron
)
zinvBVetoPt30 = cms.Sequence(
    ~countSSVHEMBJetsAK5PFNoElectronPt30
)
zinvBVetoPt50Eta25 = cms.Sequence(
    ~countSSVHEMBJetsAK5PFNoElectronPt50Eta25
)
### create the jet collections

zinvBJetsAK5PFNoElectron = cms.Sequence(
#      patSSVHEMBJetsAK5PFNoElectron
#    * patSSVHPTBJetsAK5PFNoElectron
#    * patSSVHEMBJetsAK5PFNoElectronPt30
#    * patSSVHPTBJetsAK5PFNoElectronPt30
#    * patSSVHEMBJetsAK5PFNoElectronPt50Eta25
#    * patSSVHPTBJetsAK5PFNoElectronPt50Eta25
      patCSVJetsAK5PFNoElectron
    * patCSVJetsAK5PFNoElectronPt30Eta24
    #* patCSVJetsAK5PFNoElectronPt30Eta50
    #* patCSVJetsAK5PFNoElectronPt50Eta25
    #* patCSVMVAJetsAK5PFNoElectron
    #* patCSVMVAJetsAK5PFNoElectronPt30Eta24
    #* patCSVMVAJetsAK5PFNoElectronPt30Eta50
    #* patCSVMVAJetsAK5PFNoElectronPt50Eta25
)

zinvBJetsPFNoElectron = cms.Sequence(
#      patSSVHEMBJetsPFNoElectron
#    * patSSVHPTBJetsPFNoElectron
#    * patSSVHEMBJetsPFNoElectronPt30
#    * patSSVHPTBJetsPFNoElectronPt30
#    * patSSVHEMBJetsPFNoElectronPt50Eta25
#    * patSSVHPTBJetsPFNoElectronPt50Eta25
      patCSVJetsPFNoElectron
    * patCSVJetsPFNoElectronPt30Eta24
    #* patCSVJetsPFNoElectronPt30Eta50
    #* patCSVJetsPFNoElectronPt50Eta25
    #* patCSVMVAJetsPFNoElectron
    #* patCSVMVAJetsPFNoElectronPt30Eta24
    #* patCSVMVAJetsPFNoElectronPt30Eta50
    #* patCSVMVAJetsPFNoElectronPt50Eta25
)
