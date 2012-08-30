
from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

patJetsAK5PFNoPhotonPt30 = cleanPatJets.clone()
patJetsAK5PFNoPhotonPt30.src = cms.InputTag('patJetsAK5PFPt30')
patJetsAK5PFNoPhotonPt30.checkOverlaps.photons.src             = cms.InputTag('patPhotonsAlt')
patJetsAK5PFNoPhotonPt30.checkOverlaps.photons.algorithm         = cms.string('byDeltaR')
patJetsAK5PFNoPhotonPt30.checkOverlaps.photons.preselection      = cms.string('')
patJetsAK5PFNoPhotonPt30.checkOverlaps.photons.deltaR            = cms.double(0.1)
patJetsAK5PFNoPhotonPt30.checkOverlaps.photons.pairCut           = cms.string('')
patJetsAK5PFNoPhotonPt30.checkOverlaps.photons.requireNoOverlaps = cms.bool(True)
patJetsAK5PFNoPhotonPt30.checkOverlaps.taus.src                = cms.InputTag('selectedPatTausPF')
patJetsAK5PFNoPhotonPt30.checkOverlaps.electrons.src           = cms.InputTag('patElectrons')
patJetsAK5PFNoPhotonPt30.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectrons')
patJetsAK5PFNoPhotonPt30.checkOverlaps.muons.src               = cms.InputTag('patMuons')

patJetsAK5PFNoPhotonPt30ID                         = patJetsAK5PFNoPhotonPt30.clone()
patJetsAK5PFNoPhotonPt30ID.checkOverlaps.photons.src = cms.InputTag('patPhotonsID')

patJetsAK5PFNoPhotonPt30IDIso                         = patJetsAK5PFNoPhotonPt30.clone()
patJetsAK5PFNoPhotonPt30IDIso.checkOverlaps.photons.src = cms.InputTag('patPhotonsIDPFIso')

patJetsPFNoPhotonPt30                         = patJetsAK5PFPt30.clone()
patJetsPFNoPhotonPt30.src = cms.InputTag('patJetsPF')

patJetsPFNoPhotonPt30ID                         = patJetsPFNoPhotonPt30.clone()
patJetsPFNoPhotonPt30ID.checkOverlaps.photons.src = cms.InputTag('patPhotonsID')

patJetsPFNoPhotonPt30IDIso                         = patJetsPFNoPhotonPt30.clone()
patJetsPFNoPhotonPt30IDIso.checkOverlaps.photons.src = cms.InputTag('patPhotonsIDPFIso')

#####
patJetsAK5PFNoPhotonPt50Eta25     = patJetsAK5PFNoPhotonPt30.clone()
patJetsAK5PFNoPhotonPt50Eta25.src = cms.InputTag('patJetsAK5PFPt50Eta25')

patJetsAK5PFNoPhotonPt50Eta25ID                         = patJetsAK5PFNoPhotonPt50Eta25.clone()
patJetsAK5PFNoPhotonPt50Eta25ID.checkOverlaps.photons.src = cms.InputTag('patPhotonsID')

patJetsAK5PFNoPhotonPt50Eta25IDIso                           = patJetsAK5PFNoPhotonPt50Eta25.clone()
patJetsAK5PFNoPhotonPt50Eta25IDIso.checkOverlaps.photons.src = cms.InputTag('patPhotonsIDPFIso')

patJetsPFNoPhotonPt50Eta25     = patJetsPFNoPhotonPt30.clone()
patJetsPFNoPhotonPt50Eta25.src = cms.InputTag('patJetsPFPt50Eta25')

patJetsPFNoPhotonPt50Eta25ID                         = patJetsPFNoPhotonPt50Eta25.clone()
patJetsPFNoPhotonPt50Eta25ID.checkOverlaps.photons.src = cms.InputTag('patPhotonsID')

patJetsPFNoPhotonPt50Eta25IDIso                           = patJetsPFNoPhotonPt50Eta25.clone()
patJetsPFNoPhotonPt50Eta25IDIso.checkOverlaps.photons.src = cms.InputTag('patPhotonsIDPFIso')

photonCleanedPFJetsAK5PF = cms.Sequence(
    patJetsAK5PFNoPhotonPt30           *
    patJetsAK5PFNoPhotonPt30ID         *
    patJetsAK5PFNoPhotonPt30IDIso      *
    patJetsAK5PFNoPhotonPt50Eta25      *
    patJetsAK5PFNoPhotonPt50Eta25ID    *
    patJetsAK5PFNoPhotonPt50Eta25IDIso
)
photonCleanedPFJetsPF = cms.Sequence(
    patJetsPFNoPhotonPt30           *
    patJetsPFNoPhotonPt30ID         *
    patJetsPFNoPhotonPt30IDIso      *
    patJetsPFNoPhotonPt50Eta25      *
    patJetsPFNoPhotonPt50Eta25ID    *
    patJetsPFNoPhotonPt50Eta25IDIso
)


#############

countJetsAK5PFNoPhotonPt50Eta25DiJets           = countPatJets.clone()
countJetsAK5PFNoPhotonPt50Eta25DiJets.src       = cms.InputTag('patJetsAK5PFNoPhotonPt50Eta25')
countJetsAK5PFNoPhotonPt50Eta25DiJets.minNumber = cms.uint32(2)

countJetsAK5PFNoPhotonPt50Eta25           = countPatJets.clone()
countJetsAK5PFNoPhotonPt50Eta25.src       = cms.InputTag('patJetsAK5PFNoPhotonPt50Eta25')
countJetsAK5PFNoPhotonPt50Eta25.minNumber = cms.uint32(3)

countJetsAK5PFNoPhotonPt50Eta25ID     = countJetsAK5PFNoPhotonPt50Eta25.clone()
countJetsAK5PFNoPhotonPt50Eta25ID.src = cms.InputTag('patJetsAK5PFNoPhotonPt50Eta25ID')

countJetsAK5PFNoPhotonPt50Eta25IDIso     = countJetsAK5PFNoPhotonPt50Eta25.clone()
countJetsAK5PFNoPhotonPt50Eta25IDIso.src = cms.InputTag('patJetsAK5PFNoPhotonPt50Eta25IDIso')

countJetsPFNoPhotonPt50Eta25DiJets           = countPatJets.clone()
countJetsPFNoPhotonPt50Eta25DiJets.src       = cms.InputTag('patJetsPFNoPhotonPt50Eta25')
countJetsPFNoPhotonPt50Eta25DiJets.minNumber = cms.uint32(2)

countJetsPFNoPhotonPt50Eta25           = countPatJets.clone()
countJetsPFNoPhotonPt50Eta25.src       = cms.InputTag('patJetsPFNoPhotonPt50Eta25')
countJetsPFNoPhotonPt50Eta25.minNumber = cms.uint32(3)

countJetsPFNoPhotonPt50Eta25ID     = countJetsPFNoPhotonPt50Eta25.clone()
countJetsPFNoPhotonPt50Eta25ID.src = cms.InputTag('patJetsPFNoPhotonPt50Eta25ID')

countJetsPFNoPhotonPt50Eta25IDIso     = countJetsPFNoPhotonPt50Eta25.clone()
countJetsPFNoPhotonPt50Eta25IDIso.src = cms.InputTag('patJetsPFNoPhotonPt50Eta25IDIso')


######
# b-tagged jets
#from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from ZInvisibleBkgds.Photons.specialJetSelector_cff import selectedRA2PatJets

patSSVHEMBJetsAK5PFNoPhoton = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PFNoPhoton.src = cms.InputTag('patJetsAK5PFNoPhoton')
patSSVHEMBJetsAK5PFNoPhoton.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFNoPhoton = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PFNoPhoton.src = cms.InputTag('patJetsAK5PFNoPhoton')
patSSVHPTBJetsAK5PFNoPhoton.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsAK5PFNoPhotonPt30 = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PFNoPhotonPt30.src = cms.InputTag('patJetsAK5PFNoPhotonPt30')
patSSVHEMBJetsAK5PFNoPhotonPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFNoPhotonPt30 = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PFNoPhotonPt30.src = cms.InputTag('patJetsAK5PFNoPhotonPt30')
patSSVHPTBJetsAK5PFNoPhotonPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsAK5PFNoPhotonPt50Eta25 = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PFNoPhotonPt50Eta25.src = cms.InputTag('patJetsAK5PFNoPhotonPt50Eta25')
patSSVHEMBJetsAK5PFNoPhotonPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFNoPhotonPt50Eta25 = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PFNoPhotonPt50Eta25.src = cms.InputTag('patJetsAK5PFNoPhotonPt50Eta25')
patSSVHPTBJetsAK5PFNoPhotonPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patCSVJetsAK5PFNoPhoton = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoPhoton.src = cms.InputTag('patJetsAK5PFNoPhoton')
patCSVJetsAK5PFNoPhoton.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898')

patCSVJetsAK5PFNoPhotonPt30Eta24 = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoPhotonPt30Eta24.src = cms.InputTag('patCSVJetsAK5PFNoPhoton')
patCSVJetsAK5PFNoPhotonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVJetsAK5PFNoPhotonPt30Eta50 = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoPhotonPt30Eta50.src = cms.InputTag('patCSVJetsAK5PFNoPhoton')
patCSVJetsAK5PFNoPhotonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVJetsAK5PFNoPhotonPt50Eta25 = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoPhotonPt50Eta25.src = cms.InputTag('patCSVJetsAK5PFNoPhoton')
patCSVJetsAK5PFNoPhotonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMVAJetsAK5PFNoPhoton = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoPhoton.src = cms.InputTag('patJetsAK5PFNoPhoton')
patCSVMVAJetsAK5PFNoPhoton.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > 0.898')

patCSVMVAJetsAK5PFNoPhotonPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoPhotonPt30Eta24.src = cms.InputTag('patCSVMVAJetsAK5PFNoPhoton')
patCSVMVAJetsAK5PFNoPhotonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsAK5PFNoPhotonPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoPhotonPt30Eta50.src = cms.InputTag('patCSVMVAJetsAK5PFNoPhoton')
patCSVMVAJetsAK5PFNoPhotonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsAK5PFNoPhotonPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoPhotonPt50Eta25.src = cms.InputTag('patCSVMVAJetsAK5PFNoPhoton')
patCSVMVAJetsAK5PFNoPhotonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

###patJetsPFNoPhoton
patSSVHEMBJetsPFNoPhoton = selectedRA2PatJets.clone()
patSSVHEMBJetsPFNoPhoton.src = cms.InputTag('patJetsPFNoPhoton')
patSSVHEMBJetsPFNoPhoton.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPFNoPhoton = selectedRA2PatJets.clone()
patSSVHPTBJetsPFNoPhoton.src = cms.InputTag('patJetsPFNoPhoton')
patSSVHPTBJetsPFNoPhoton.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsPFNoPhotonPt30 = selectedRA2PatJets.clone()
patSSVHEMBJetsPFNoPhotonPt30.src = cms.InputTag('patJetsPFNoPhotonPt30')
patSSVHEMBJetsPFNoPhotonPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPFNoPhotonPt30 = selectedRA2PatJets.clone()
patSSVHPTBJetsPFNoPhotonPt30.src = cms.InputTag('patJetsPFNoPhotonPt30')
patSSVHPTBJetsPFNoPhotonPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsPFNoPhotonPt50Eta25 = selectedRA2PatJets.clone()
patSSVHEMBJetsPFNoPhotonPt50Eta25.src = cms.InputTag('patJetsPFNoPhotonPt50Eta25')
patSSVHEMBJetsPFNoPhotonPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPFNoPhotonPt50Eta25 = selectedRA2PatJets.clone()
patSSVHPTBJetsPFNoPhotonPt50Eta25.src = cms.InputTag('patJetsPFNoPhotonPt50Eta25')
patSSVHPTBJetsPFNoPhotonPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patCSVJetsPFNoPhoton = selectedRA2PatJets.clone()
patCSVJetsPFNoPhoton.src = cms.InputTag('patJetsPFNoPhoton')
patCSVJetsPFNoPhoton.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898')

patCSVJetsPFNoPhotonPt30Eta24 = selectedRA2PatJets.clone()
patCSVJetsPFNoPhotonPt30Eta24.src = cms.InputTag('patCSVJetsPFNoPhoton')
patCSVJetsPFNoPhotonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVJetsPFNoPhotonPt30Eta50 = selectedRA2PatJets.clone()
patCSVJetsPFNoPhotonPt30Eta50.src = cms.InputTag('patCSVJetsPFNoPhoton')
patCSVJetsPFNoPhotonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVJetsPFNoPhotonPt50Eta25 = selectedRA2PatJets.clone()
patCSVJetsPFNoPhotonPt50Eta25.src = cms.InputTag('patCSVJetsPFNoPhoton')
patCSVJetsPFNoPhotonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMVAJetsPFNoPhoton = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoPhoton.src = cms.InputTag('patJetsPFNoPhoton')
patCSVMVAJetsPFNoPhoton.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > 0.898')

patCSVMVAJetsPFNoPhotonPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoPhotonPt30Eta24.src = cms.InputTag('patCSVMVAJetsPFNoPhoton')
patCSVMVAJetsPFNoPhotonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsPFNoPhotonPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoPhotonPt30Eta50.src = cms.InputTag('patCSVMVAJetsPFNoPhoton')
patCSVMVAJetsPFNoPhotonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsPFNoPhotonPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoPhotonPt50Eta25.src = cms.InputTag('patCSVMVAJetsPFNoPhoton')
patCSVMVAJetsPFNoPhotonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

### count the b-jets
countSSVHEMBJetsAK5PFNoPhoton = countPatJets.clone()
countSSVHEMBJetsAK5PFNoPhoton.src = cms.InputTag('patSSVHEMBJetsAK5PFNoPhoton')
countSSVHEMBJetsAK5PFNoPhoton.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PFNoPhoton = countPatJets.clone()
countSSVHPTBJetsAK5PFNoPhoton.src = cms.InputTag('patSSVHPTBJetsAK5PFNoPhoton')
countSSVHPTBJetsAK5PFNoPhoton.minNumber = cms.uint32(1)

countSSVHEMBJetsAK5PFNoPhotonPt30 = countPatJets.clone()
countSSVHEMBJetsAK5PFNoPhotonPt30.src = cms.InputTag('patSSVHEMBJetsAK5PFNoPhotonPt30')
countSSVHEMBJetsAK5PFNoPhotonPt30.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PFNoPhotonPt30 = countPatJets.clone()
countSSVHPTBJetsAK5PFNoPhotonPt30.src = cms.InputTag('patSSVHPTBJetsAK5PFNoPhotonPt30')
countSSVHPTBJetsAK5PFNoPhotonPt30.minNumber = cms.uint32(1)

countSSVHEMBJetsAK5PFNoPhotonPt50Eta25 = countPatJets.clone()
countSSVHEMBJetsAK5PFNoPhotonPt50Eta25.src = cms.InputTag('patSSVHEMBJetsAK5PFNoPhotonPt50Eta25')
countSSVHEMBJetsAK5PFNoPhotonPt50Eta25.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PFNoPhotonPt50Eta25 = countPatJets.clone()
countSSVHPTBJetsAK5PFNoPhotonPt50Eta25.src = cms.InputTag('patSSVHPTBJetsAK5PFNoPhotonPt50Eta25')
countSSVHPTBJetsAK5PFNoPhotonPt50Eta25.minNumber = cms.uint32(1)

###patJetsPFNoPhoton
countSSVHEMBJetsPFNoPhoton = countPatJets.clone()
countSSVHEMBJetsPFNoPhoton.src = cms.InputTag('patSSVHEMBJetsPFNoPhoton')
countSSVHEMBJetsPFNoPhoton.minNumber = cms.uint32(1)

countSSVHPTBJetsPFNoPhoton = countPatJets.clone()
countSSVHPTBJetsPFNoPhoton.src = cms.InputTag('patSSVHPTBJetsPFNoPhoton')
countSSVHPTBJetsPFNoPhoton.minNumber = cms.uint32(1)

countSSVHEMBJetsPFNoPhotonPt30 = countPatJets.clone()
countSSVHEMBJetsPFNoPhotonPt30.src = cms.InputTag('patSSVHEMBJetsPFNoPhotonPt30')
countSSVHEMBJetsPFNoPhotonPt30.minNumber = cms.uint32(1)

countSSVHPTBJetsPFNoPhotonPt30 = countPatJets.clone()
countSSVHPTBJetsPFNoPhotonPt30.src = cms.InputTag('patSSVHPTBJetsPFNoPhotonPt30')
countSSVHPTBJetsPFNoPhotonPt30.minNumber = cms.uint32(1)

countSSVHEMBJetsPFNoPhotonPt50Eta25 = countPatJets.clone()
countSSVHEMBJetsPFNoPhotonPt50Eta25.src = cms.InputTag('patSSVHEMBJetsPFNoPhotonPt50Eta25')
countSSVHEMBJetsPFNoPhotonPt50Eta25.minNumber = cms.uint32(1)

countSSVHPTBJetsPFNoPhotonPt50Eta25 = countPatJets.clone()
countSSVHPTBJetsPFNoPhotonPt50Eta25.src = cms.InputTag('patSSVHPTBJetsPFNoPhotonPt50Eta25')
countSSVHPTBJetsPFNoPhotonPt50Eta25.minNumber = cms.uint32(1)


zinvBVetoNoPhoton = cms.Sequence(
    ~countSSVHEMBJetsAK5PFNoPhoton
)
zinvBVetoNoPhotonPt30 = cms.Sequence(
    ~countSSVHEMBJetsAK5PFNoPhotonPt30
)
zinvBVetoPtNoPhoton50Eta25 = cms.Sequence(
    ~countSSVHEMBJetsAK5PFNoPhotonPt50Eta25
)

### create the jet collections
zinvBJetsAK5PFNoPhoton = cms.Sequence(
#      patSSVHEMBJetsAK5PFNoPhoton
#    * patSSVHPTBJetsAK5PFNoPhoton
#    * patSSVHEMBJetsAK5PFNoPhotonPt30
#    * patSSVHPTBJetsAK5PFNoPhotonPt30
#    * patSSVHEMBJetsAK5PFNoPhotonPt50Eta25
#    * patSSVHPTBJetsAK5PFNoPhotonPt50Eta25
      patCSVJetsAK5PFNoPhoton
    * patCSVJetsAK5PFNoPhotonPt30Eta24
    #* patCSVJetsAK5PFNoPhotonPt30Eta50
    #* patCSVJetsAK5PFNoPhotonPt50Eta25
    #* patCSVMVAJetsAK5PFNoPhoton
    #* patCSVMVAJetsAK5PFNoPhotonPt30Eta24
    #* patCSVMVAJetsAK5PFNoPhotonPt30Eta50
    #* patCSVMVAJetsAK5PFNoPhotonPt50Eta25
)

zinvBJetsPFNoPhoton = cms.Sequence(
#      patSSVHEMBJetsPFNoPhoton
#    * patSSVHPTBJetsPFNoPhoton
#    * patSSVHEMBJetsPFNoPhotonPt30
#    * patSSVHPTBJetsPFNoPhotonPt30
#    * patSSVHEMBJetsPFNoPhotonPt50Eta25
#    * patSSVHPTBJetsPFNoPhotonPt50Eta25
      patCSVJetsPFNoPhoton
    * patCSVJetsPFNoPhotonPt30Eta24
    #* patCSVJetsPFNoPhotonPt30Eta50
    #* patCSVJetsPFNoPhotonPt50Eta25
    #* patCSVMVAJetsPFNoPhoton
    #* patCSVMVAJetsPFNoPhotonPt30Eta24
    #* patCSVMVAJetsPFNoPhotonPt30Eta50
    #* patCSVMVAJetsPFNoPhotonPt50Eta25
)

#############
from ZInvisibleBkgds.Photons.specialObjectCleaner_cff import specialPhotonCleanedJets

patJetsAK5PFNoPhotonSpecial              = specialPhotonCleanedJets.clone()
patJetsAK5PFNoPhotonSpecial.jetLabel     = cms.InputTag('patJetsAK5PF')
patJetsAK5PFNoPhotonSpecial.objectLabel  = cms.InputTag('patPhotonsAlt')

patJetsAK5PFNoPhotonSpecialID             = patJetsAK5PFNoPhotonSpecial.clone()
patJetsAK5PFNoPhotonSpecialID.objectLabel = cms.InputTag('patPhotonsID')

patJetsAK5PFNoPhotonSpecialIDIso             = patJetsAK5PFNoPhotonSpecial.clone()
patJetsAK5PFNoPhotonSpecialIDIso.objectLabel = cms.InputTag('patPhotonsIDPFIso')

patJetsPFNoPhotonSpecial     = patJetsAK5PF.clone()
patJetsPFNoPhotonSpecial.src = cms.InputTag('patJetsPF')

patJetsPFNoPhotonSpecialID             = patJetsPFNoPhotonSpecial.clone()
patJetsPFNoPhotonSpecialID.objectLabel = cms.InputTag('patPhotonsID')

patJetsPFNoPhotonSpecialIDIso             = patJetsPFNoPhotonSpecial.clone()
patJetsPFNoPhotonSpecialIDIso.objectLabel = cms.InputTag('patPhotonsIDPFIso')

#MHT Jets
patJetsAK5PFNoPhotonSpecialPt30              = patJetsAK5PFNoPhotonSpecial.clone()
patJetsAK5PFNoPhotonSpecialPt30.jetLabel     = cms.InputTag('patJetsAK5PFPt30')

patJetsAK5PFNoPhotonSpecialPt30ID             = patJetsAK5PFNoPhotonSpecialPt30.clone()
patJetsAK5PFNoPhotonSpecialPt30ID.objectLabel = cms.InputTag('patPhotonsID')

patJetsAK5PFNoPhotonSpecialPt30IDIso             = patJetsAK5PFNoPhotonSpecialPt30.clone()
patJetsAK5PFNoPhotonSpecialPt30IDIso.objectLabel = cms.InputTag('patPhotonsIDPFIso')

patJetsPFNoPhotonSpecialPt30          = patJetsAK5PFPt30.clone()
patJetsPFNoPhotonSpecialPt30.jetLabel = cms.InputTag('patJetsPFPt30')

patJetsPFNoPhotonSpecialPt30ID             = patJetsPFNoPhotonSpecialPt30.clone()
patJetsPFNoPhotonSpecialPt30ID.objectLabel = cms.InputTag('patPhotonsID')

patJetsPFNoPhotonSpecialPt30IDIso             = patJetsPFNoPhotonSpecialPt30.clone()
patJetsPFNoPhotonSpecialPt30IDIso.objectLabel = cms.InputTag('patPhotonsIDPFIso')

#HT Jets
patJetsAK5PFNoPhotonSpecialPt50Eta25          = patJetsAK5PFNoPhotonSpecialPt30.clone()
patJetsAK5PFNoPhotonSpecialPt50Eta25.jetLabel = cms.InputTag('patJetsAK5PFPt50Eta25')

patJetsAK5PFNoPhotonSpecialPt50Eta25ID             = patJetsAK5PFNoPhotonSpecialPt50Eta25.clone()
patJetsAK5PFNoPhotonSpecialPt50Eta25ID.objectLabel = cms.InputTag('patPhotonsID')

patJetsAK5PFNoPhotonSpecialPt50Eta25IDIso             = patJetsAK5PFNoPhotonSpecialPt50Eta25.clone()
patJetsAK5PFNoPhotonSpecialPt50Eta25IDIso.objectLabel = cms.InputTag('patPhotonsIDPFIso')

patJetsPFNoPhotonSpecialPt50Eta25          = patJetsPFNoPhotonSpecialPt30.clone()
patJetsPFNoPhotonSpecialPt50Eta25.jetLabel = cms.InputTag('patJetsPFPt50Eta25')

patJetsPFNoPhotonSpecialPt50Eta25ID             = patJetsPFNoPhotonSpecialPt50Eta25.clone()
patJetsPFNoPhotonSpecialPt50Eta25ID.objectLabel = cms.InputTag('patPhotonsID')

patJetsPFNoPhotonSpecialPt50Eta25IDIso             = patJetsPFNoPhotonSpecialPt50Eta25.clone()
patJetsPFNoPhotonSpecialPt50Eta25IDIso.objectLabel = cms.InputTag('patPhotonsIDPFIso')

specialPhotonCleanedPFJetsAK5PF = cms.Sequence(
    patJetsAK5PFNoPhotonSpecial           *
    patJetsAK5PFNoPhotonSpecialID         *
    patJetsAK5PFNoPhotonSpecialIDIso      *
    patJetsAK5PFNoPhotonSpecialPt30           *
    patJetsAK5PFNoPhotonSpecialPt30ID         *
    patJetsAK5PFNoPhotonSpecialPt30IDIso      *
    patJetsAK5PFNoPhotonSpecialPt50Eta25      *
    patJetsAK5PFNoPhotonSpecialPt50Eta25ID    *
    patJetsAK5PFNoPhotonSpecialPt50Eta25IDIso
)
specialPhotonCleanedPFJetsPF = cms.Sequence(
    patJetsPFNoPhotonSpecial           *
    patJetsPFNoPhotonSpecialID         *
    patJetsPFNoPhotonSpecialIDIso      *
    patJetsPFNoPhotonSpecialPt30           *
    patJetsPFNoPhotonSpecialPt30ID         *
    patJetsPFNoPhotonSpecialPt30IDIso      *
    patJetsPFNoPhotonSpecialPt50Eta25      *
    patJetsPFNoPhotonSpecialPt50Eta25ID    *
    patJetsPFNoPhotonSpecialPt50Eta25IDIso
)

#############

countJetsAK5PFNoPhotonSpecialPt50Eta25DiJets           = countPatJets.clone()
countJetsAK5PFNoPhotonSpecialPt50Eta25DiJets.src       = cms.InputTag('patJetsAK5PFNoPhotonSpecialPt50Eta25')
countJetsAK5PFNoPhotonSpecialPt50Eta25DiJets.minNumber = cms.uint32(2)

countJetsAK5PFNoPhotonSpecialPt50Eta25           = countPatJets.clone()
countJetsAK5PFNoPhotonSpecialPt50Eta25.src       = cms.InputTag('patJetsAK5PFNoPhotonSpecialPt50Eta25')
countJetsAK5PFNoPhotonSpecialPt50Eta25.minNumber = cms.uint32(3)

countJetsAK5PFNoPhotonSpecialPt50Eta25ID     = countJetsAK5PFNoPhotonSpecialPt50Eta25.clone()
countJetsAK5PFNoPhotonSpecialPt50Eta25ID.src = cms.InputTag('patJetsAK5PFNoPhotonSpecialPt50Eta25ID')

countJetsAK5PFNoPhotonSpecialPt50Eta25IDIso     = countJetsAK5PFNoPhotonSpecialPt50Eta25.clone()
countJetsAK5PFNoPhotonSpecialPt50Eta25IDIso.src = cms.InputTag('patJetsAK5PFNoPhotonSpecialPt50Eta25IDIso')

countJetsPFNoPhotonSpecialPt50Eta25DiJets           = countPatJets.clone()
countJetsPFNoPhotonSpecialPt50Eta25DiJets.src       = cms.InputTag('patJetsPFNoPhotonSpecialPt50Eta25')
countJetsPFNoPhotonSpecialPt50Eta25DiJets.minNumber = cms.uint32(2)

countJetsPFNoPhotonSpecialPt50Eta25           = countPatJets.clone()
countJetsPFNoPhotonSpecialPt50Eta25.src       = cms.InputTag('patJetsPFNoPhotonSpecialPt50Eta25')
countJetsPFNoPhotonSpecialPt50Eta25.minNumber = cms.uint32(3)

countJetsPFNoPhotonSpecialPt50Eta25ID     = countJetsPFNoPhotonSpecialPt50Eta25.clone()
countJetsPFNoPhotonSpecialPt50Eta25ID.src = cms.InputTag('patJetsPFNoPhotonSpecialPt50Eta25ID')

countJetsPFNoPhotonSpecialPt50Eta25IDIso     = countJetsPFNoPhotonSpecialPt50Eta25.clone()
countJetsPFNoPhotonSpecialPt50Eta25IDIso.src = cms.InputTag('patJetsPFNoPhotonSpecialPt50Eta25IDIso')


######
# b-tagged jets
#from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from ZInvisibleBkgds.Photons.specialJetSelector_cff import selectedRA2PatJets

patSSVHEMBJetsAK5PFNoPhotonSpecial = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PFNoPhotonSpecial.src = cms.InputTag('patJetsAK5PFNoPhotonSpecial')
patSSVHEMBJetsAK5PFNoPhotonSpecial.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFNoPhotonSpecial = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PFNoPhotonSpecial.src = cms.InputTag('patJetsAK5PFNoPhotonSpecial')
patSSVHPTBJetsAK5PFNoPhotonSpecial.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsAK5PFNoPhotonSpecialPt30 = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PFNoPhotonSpecialPt30.src = cms.InputTag('patJetsAK5PFNoPhotonSpecialPt30')
patSSVHEMBJetsAK5PFNoPhotonSpecialPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFNoPhotonSpecialPt30 = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PFNoPhotonSpecialPt30.src = cms.InputTag('patJetsAK5PFNoPhotonSpecialPt30')
patSSVHPTBJetsAK5PFNoPhotonSpecialPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsAK5PFNoPhotonSpecialPt50Eta25 = selectedRA2PatJets.clone()
patSSVHEMBJetsAK5PFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patJetsAK5PFNoPhotonSpecialPt50Eta25')
patSSVHEMBJetsAK5PFNoPhotonSpecialPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsAK5PFNoPhotonSpecialPt50Eta25 = selectedRA2PatJets.clone()
patSSVHPTBJetsAK5PFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patJetsAK5PFNoPhotonSpecialPt50Eta25')
patSSVHPTBJetsAK5PFNoPhotonSpecialPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patCSVJetsAK5PFNoPhotonSpecial = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoPhotonSpecial.src = cms.InputTag('patJetsAK5PFNoPhotonSpecial')
patCSVJetsAK5PFNoPhotonSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898')

patCSVJetsAK5PFNoPhotonSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoPhotonSpecialPt30Eta24.src = cms.InputTag('patCSVJetsAK5PFNoPhotonSpecial')
patCSVJetsAK5PFNoPhotonSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVJetsAK5PFNoPhotonSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoPhotonSpecialPt30Eta50.src = cms.InputTag('patCSVJetsAK5PFNoPhotonSpecial')
patCSVJetsAK5PFNoPhotonSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVJetsAK5PFNoPhotonSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVJetsAK5PFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patCSVJetsAK5PFNoPhotonSpecial')
patCSVJetsAK5PFNoPhotonSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMVAJetsAK5PFNoPhotonSpecial = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoPhotonSpecial.src = cms.InputTag('patJetsAK5PFNoPhotonSpecial')
patCSVMVAJetsAK5PFNoPhotonSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > 0.898')

patCSVMVAJetsAK5PFNoPhotonSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoPhotonSpecialPt30Eta24.src = cms.InputTag('patCSVMVAJetsAK5PFNoPhotonSpecial')
patCSVMVAJetsAK5PFNoPhotonSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsAK5PFNoPhotonSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoPhotonSpecialPt30Eta50.src = cms.InputTag('patCSVMVAJetsAK5PFNoPhotonSpecial')
patCSVMVAJetsAK5PFNoPhotonSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsAK5PFNoPhotonSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patCSVMVAJetsAK5PFNoPhotonSpecial')
patCSVMVAJetsAK5PFNoPhotonSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

###patJetsPFNoPhotonSpecial
patSSVHEMBJetsPFNoPhotonSpecial = selectedRA2PatJets.clone()
patSSVHEMBJetsPFNoPhotonSpecial.src = cms.InputTag('patJetsPFNoPhotonSpecial')
patSSVHEMBJetsPFNoPhotonSpecial.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPFNoPhotonSpecial = selectedRA2PatJets.clone()
patSSVHPTBJetsPFNoPhotonSpecial.src = cms.InputTag('patJetsPFNoPhotonSpecial')
patSSVHPTBJetsPFNoPhotonSpecial.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsPFNoPhotonSpecialPt30 = selectedRA2PatJets.clone()
patSSVHEMBJetsPFNoPhotonSpecialPt30.src = cms.InputTag('patJetsPFNoPhotonSpecialPt30')
patSSVHEMBJetsPFNoPhotonSpecialPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPFNoPhotonSpecialPt30 = selectedRA2PatJets.clone()
patSSVHPTBJetsPFNoPhotonSpecialPt30.src = cms.InputTag('patJetsPFNoPhotonSpecialPt30')
patSSVHPTBJetsPFNoPhotonSpecialPt30.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patSSVHEMBJetsPFNoPhotonSpecialPt50Eta25 = selectedRA2PatJets.clone()
patSSVHEMBJetsPFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patJetsPFNoPhotonSpecialPt50Eta25')
patSSVHEMBJetsPFNoPhotonSpecialPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighEffBJetTags") > 1.74')

patSSVHPTBJetsPFNoPhotonSpecialPt50Eta25 = selectedRA2PatJets.clone()
patSSVHPTBJetsPFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patJetsPFNoPhotonSpecialPt50Eta25')
patSSVHPTBJetsPFNoPhotonSpecialPt50Eta25.cut = cms.string('bDiscriminator("simpleSecondaryVertexHighPurBJetTags") > 2.00')

patCSVJetsPFNoPhotonSpecial = selectedRA2PatJets.clone()
patCSVJetsPFNoPhotonSpecial.src = cms.InputTag('patJetsPFNoPhotonSpecial')
patCSVJetsPFNoPhotonSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898')

patCSVJetsPFNoPhotonSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVJetsPFNoPhotonSpecialPt30Eta24.src = cms.InputTag('patCSVJetsPFNoPhotonSpecial')
patCSVJetsPFNoPhotonSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVJetsPFNoPhotonSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVJetsPFNoPhotonSpecialPt30Eta50.src = cms.InputTag('patCSVJetsPFNoPhotonSpecial')
patCSVJetsPFNoPhotonSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVJetsPFNoPhotonSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVJetsPFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patCSVJetsPFNoPhotonSpecial')
patCSVJetsPFNoPhotonSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMVAJetsPFNoPhotonSpecial = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoPhotonSpecial.src = cms.InputTag('patJetsPFNoPhotonSpecial')
patCSVMVAJetsPFNoPhotonSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > 0.898')

patCSVMVAJetsPFNoPhotonSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoPhotonSpecialPt30Eta24.src = cms.InputTag('patCSVMVAJetsPFNoPhotonSpecial')
patCSVMVAJetsPFNoPhotonSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsPFNoPhotonSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoPhotonSpecialPt30Eta50.src = cms.InputTag('patCSVMVAJetsPFNoPhotonSpecial')
patCSVMVAJetsPFNoPhotonSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsPFNoPhotonSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patCSVMVAJetsPFNoPhotonSpecial')
patCSVMVAJetsPFNoPhotonSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

### count the b-jets
countSSVHEMBJetsAK5PFNoPhotonSpecial = countPatJets.clone()
countSSVHEMBJetsAK5PFNoPhotonSpecial.src = cms.InputTag('patSSVHEMBJetsAK5PFNoPhotonSpecial')
countSSVHEMBJetsAK5PFNoPhotonSpecial.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PFNoPhotonSpecial = countPatJets.clone()
countSSVHPTBJetsAK5PFNoPhotonSpecial.src = cms.InputTag('patSSVHPTBJetsAK5PFNoPhotonSpecial')
countSSVHPTBJetsAK5PFNoPhotonSpecial.minNumber = cms.uint32(1)

countSSVHEMBJetsAK5PFNoPhotonSpecialPt30 = countPatJets.clone()
countSSVHEMBJetsAK5PFNoPhotonSpecialPt30.src = cms.InputTag('patSSVHEMBJetsAK5PFNoPhotonSpecialPt30')
countSSVHEMBJetsAK5PFNoPhotonSpecialPt30.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PFNoPhotonSpecialPt30 = countPatJets.clone()
countSSVHPTBJetsAK5PFNoPhotonSpecialPt30.src = cms.InputTag('patSSVHPTBJetsAK5PFNoPhotonSpecialPt30')
countSSVHPTBJetsAK5PFNoPhotonSpecialPt30.minNumber = cms.uint32(1)

countSSVHEMBJetsAK5PFNoPhotonSpecialPt50Eta25 = countPatJets.clone()
countSSVHEMBJetsAK5PFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patSSVHEMBJetsAK5PFNoPhotonSpecialPt50Eta25')
countSSVHEMBJetsAK5PFNoPhotonSpecialPt50Eta25.minNumber = cms.uint32(1)

countSSVHPTBJetsAK5PFNoPhotonSpecialPt50Eta25 = countPatJets.clone()
countSSVHPTBJetsAK5PFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patSSVHPTBJetsAK5PFNoPhotonSpecialPt50Eta25')
countSSVHPTBJetsAK5PFNoPhotonSpecialPt50Eta25.minNumber = cms.uint32(1)

###patJetsPFNoPhotonSpecial
countSSVHEMBJetsPFNoPhotonSpecial = countPatJets.clone()
countSSVHEMBJetsPFNoPhotonSpecial.src = cms.InputTag('patSSVHEMBJetsPFNoPhotonSpecial')
countSSVHEMBJetsPFNoPhotonSpecial.minNumber = cms.uint32(1)

countSSVHPTBJetsPFNoPhotonSpecial = countPatJets.clone()
countSSVHPTBJetsPFNoPhotonSpecial.src = cms.InputTag('patSSVHPTBJetsPFNoPhotonSpecial')
countSSVHPTBJetsPFNoPhotonSpecial.minNumber = cms.uint32(1)

countSSVHEMBJetsPFNoPhotonSpecialPt30 = countPatJets.clone()
countSSVHEMBJetsPFNoPhotonSpecialPt30.src = cms.InputTag('patSSVHEMBJetsPFNoPhotonSpecialPt30')
countSSVHEMBJetsPFNoPhotonSpecialPt30.minNumber = cms.uint32(1)

countSSVHPTBJetsPFNoPhotonSpecialPt30 = countPatJets.clone()
countSSVHPTBJetsPFNoPhotonSpecialPt30.src = cms.InputTag('patSSVHPTBJetsPFNoPhotonSpecialPt30')
countSSVHPTBJetsPFNoPhotonSpecialPt30.minNumber = cms.uint32(1)

countSSVHEMBJetsPFNoPhotonSpecialPt50Eta25 = countPatJets.clone()
countSSVHEMBJetsPFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patSSVHEMBJetsPFNoPhotonSpecialPt50Eta25')
countSSVHEMBJetsPFNoPhotonSpecialPt50Eta25.minNumber = cms.uint32(1)

countSSVHPTBJetsPFNoPhotonSpecialPt50Eta25 = countPatJets.clone()
countSSVHPTBJetsPFNoPhotonSpecialPt50Eta25.src = cms.InputTag('patSSVHPTBJetsPFNoPhotonSpecialPt50Eta25')
countSSVHPTBJetsPFNoPhotonSpecialPt50Eta25.minNumber = cms.uint32(1)


zinvBVetoNoPhotonSpecial = cms.Sequence(
    ~countSSVHEMBJetsAK5PFNoPhotonSpecial
)
zinvBVetoNoPhotonSpecialPt30 = cms.Sequence(
    ~countSSVHEMBJetsAK5PFNoPhotonSpecialPt30
)
zinvBVetoNoPhotonSpecialPt50Eta25 = cms.Sequence(
    ~countSSVHEMBJetsAK5PFNoPhotonSpecialPt50Eta25
)

### create the jet collections
zinvBJetsAK5PFNoPhotonSpecial = cms.Sequence(
    #  patSSVHEMBJetsAK5PFNoPhotonSpecial
    #* patSSVHPTBJetsAK5PFNoPhotonSpecial
    #* patSSVHEMBJetsAK5PFNoPhotonSpecialPt30
    #* patSSVHPTBJetsAK5PFNoPhotonSpecialPt30
    #* patSSVHEMBJetsAK5PFNoPhotonSpecialPt50Eta25
    #* patSSVHPTBJetsAK5PFNoPhotonSpecialPt50Eta25
      patCSVJetsAK5PFNoPhotonSpecial
    * patCSVJetsAK5PFNoPhotonSpecialPt30Eta24
    #* patCSVJetsAK5PFNoPhotonSpecialPt30Eta50
    #* patCSVJetsAK5PFNoPhotonSpecialPt50Eta25
    #* patCSVMVAJetsAK5PFNoPhotonSpecial
    #* patCSVMVAJetsAK5PFNoPhotonSpecialPt30Eta24
    #* patCSVMVAJetsAK5PFNoPhotonSpecialPt30Eta50
    #* patCSVMVAJetsAK5PFNoPhotonSpecialPt50Eta25
)

zinvBJetsPFNoPhotonSpecial = cms.Sequence(
    #  patSSVHEMBJetsPFNoPhotonSpecial
    #* patSSVHPTBJetsPFNoPhotonSpecial
    #* patSSVHEMBJetsPFNoPhotonSpecialPt30
    #* patSSVHPTBJetsPFNoPhotonSpecialPt30
    #* patSSVHEMBJetsPFNoPhotonSpecialPt50Eta25
    #* patSSVHPTBJetsPFNoPhotonSpecialPt50Eta25
      patCSVJetsPFNoPhotonSpecial
    * patCSVJetsPFNoPhotonSpecialPt30Eta24
    #* patCSVJetsPFNoPhotonSpecialPt30Eta50
    #* patCSVJetsPFNoPhotonSpecialPt50Eta25
    #* patCSVMVAJetsPFNoPhotonSpecial
    #* patCSVMVAJetsPFNoPhotonSpecialPt30Eta24
    #* patCSVMVAJetsPFNoPhotonSpecialPt30Eta50
    #* patCSVMVAJetsPFNoPhotonSpecialPt50Eta25
)
