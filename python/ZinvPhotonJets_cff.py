import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

patJetsAK5PFNoPhotonID = cleanPatJets.clone()
patJetsAK5PFNoPhotonID.src = cms.InputTag('patJetsAK5PF')
patJetsAK5PFNoPhotonID.checkOverlaps.photons.src               = cms.InputTag('patPhotonsID')
patJetsAK5PFNoPhotonID.checkOverlaps.photons.algorithm         = cms.string('byDeltaR')
patJetsAK5PFNoPhotonID.checkOverlaps.photons.preselection      = cms.string('')
patJetsAK5PFNoPhotonID.checkOverlaps.photons.deltaR            = cms.double(0.1)
patJetsAK5PFNoPhotonID.checkOverlaps.photons.pairCut           = cms.string('')
patJetsAK5PFNoPhotonID.checkOverlaps.photons.requireNoOverlaps = cms.bool(True)
patJetsAK5PFNoPhotonID.checkOverlaps.taus.src                = cms.InputTag('selectedPatTausPF')
patJetsAK5PFNoPhotonID.checkOverlaps.electrons.src           = cms.InputTag('patElectronsPF')
patJetsAK5PFNoPhotonID.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectronsPF')
patJetsAK5PFNoPhotonID.checkOverlaps.muons.src               = cms.InputTag('patMuonsPF')

patJetsAK5PFNoPhotonIDPt30 = patJetsAK5PFNoPhotonID.clone()
patJetsAK5PFNoPhotonIDPt30.src = cms.InputTag('patJetsAK5PFPt30')

patJetsPFNoPhotonID     = patJetsAK5PFNoPhotonID.clone()
patJetsPFNoPhotonID.src = cms.InputTag('patJetsPF')

patJetsPFNoPhotonIDPt30                         = patJetsAK5PFNoPhotonIDPt30.clone()
patJetsPFNoPhotonIDPt30.src = cms.InputTag('patJetsPFchsPt30')

#####
patJetsAK5PFNoPhotonIDPt50Eta25     = patJetsAK5PFNoPhotonIDPt30.clone()
patJetsAK5PFNoPhotonIDPt50Eta25.src = cms.InputTag('patJetsAK5PFPt50Eta25')

patJetsPFNoPhotonIDPt50Eta25     = patJetsPFNoPhotonIDPt30.clone()
patJetsPFNoPhotonIDPt50Eta25.src = cms.InputTag('patJetsPFchsPt50Eta25')

####ID/PF Iso
patJetsAK5PFNoPhotonIDPFIso = patJetsAK5PFNoPhotonID.clone()
patJetsAK5PFNoPhotonIDPFIso.src = cms.InputTag('patJetsAK5PF')
patJetsAK5PFNoPhotonIDPFIso.checkOverlaps.photons.src             = cms.InputTag('patPhotonsIDPFIso')

patJetsAK5PFNoPhotonIDPFIsoPt30 = patJetsAK5PFNoPhotonIDPFIso.clone()
patJetsAK5PFNoPhotonIDPFIsoPt30.src = cms.InputTag('patJetsAK5PFPt30')

patJetsPFNoPhotonIDPFIso     = patJetsAK5PFNoPhotonIDPFIso.clone()
patJetsPFNoPhotonIDPFIso.src = cms.InputTag('patJetsPF')

patJetsPFNoPhotonIDPFIsoPt30                         = patJetsAK5PFNoPhotonIDPFIsoPt30.clone()
patJetsPFNoPhotonIDPFIsoPt30.src = cms.InputTag('patJetsPFchsPt30')

#####
patJetsAK5PFNoPhotonIDPFIsoPt50Eta25     = patJetsAK5PFNoPhotonIDPFIsoPt30.clone()
patJetsAK5PFNoPhotonIDPFIsoPt50Eta25.src = cms.InputTag('patJetsAK5PFPt50Eta25')

patJetsPFNoPhotonIDPFIsoPt50Eta25     = patJetsPFNoPhotonIDPFIsoPt30.clone()
patJetsPFNoPhotonIDPFIsoPt50Eta25.src = cms.InputTag('patJetsPFchsPt50Eta25')


#####
photonCleanedPFJetsAK5PF = cms.Sequence(
    patJetsAK5PFNoPhotonID               *
    patJetsAK5PFNoPhotonIDPt30           *
    patJetsAK5PFNoPhotonIDPt50Eta25      *
    patJetsAK5PFNoPhotonIDPFIso          *
    patJetsAK5PFNoPhotonIDPFIsoPt30      *
    patJetsAK5PFNoPhotonIDPFIsoPt50Eta25
)
photonCleanedPFJetsPF = cms.Sequence(
    patJetsPFNoPhotonID               *
    patJetsPFNoPhotonIDPt30           *
    patJetsPFNoPhotonIDPt50Eta25      *
    patJetsPFNoPhotonIDPFIso          *
    patJetsPFNoPhotonIDPFIsoPt30      *
    patJetsPFNoPhotonIDPFIsoPt50Eta25
)


#############

countJetsAK5PFNoPhotonIDPFIsoPt50Eta25DiJets           = countPatJets.clone()
countJetsAK5PFNoPhotonIDPFIsoPt50Eta25DiJets.src       = cms.InputTag('patJetsAK5PFNoPhotonIDPFIsoPt50Eta25')
countJetsAK5PFNoPhotonIDPFIsoPt50Eta25DiJets.minNumber = cms.uint32(2)

countJetsAK5PFNoPhotonIDPFIsoPt50Eta25           = countPatJets.clone()
countJetsAK5PFNoPhotonIDPFIsoPt50Eta25.src       = cms.InputTag('patJetsAK5PFNoPhotonIDPFIsoPt50Eta25')
countJetsAK5PFNoPhotonIDPFIsoPt50Eta25.minNumber = cms.uint32(3)

countJetsPFNoPhotonIDPFIsoPt50Eta25DiJets           = countPatJets.clone()
countJetsPFNoPhotonIDPFIsoPt50Eta25DiJets.src       = cms.InputTag('patJetsPFNoPhotonIDPFIsoPt50Eta25')
countJetsPFNoPhotonIDPFIsoPt50Eta25DiJets.minNumber = cms.uint32(2)

countJetsPFNoPhotonIDPFIsoPt50Eta25           = countPatJets.clone()
countJetsPFNoPhotonIDPFIsoPt50Eta25.src       = cms.InputTag('patJetsPFNoPhotonIDPFIsoPt50Eta25')
countJetsPFNoPhotonIDPFIsoPt50Eta25.minNumber = cms.uint32(3)


######
# b-tagged jets
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from ZInvisibleBkgds.Photons.specialJetSelector_cff import selectedRA2PatJets
CSVL = 0.244
CSVM = 0.679
CSVT = 0.898 
csvPoint = CSVT

patCSVMJetsAK5PFNoPhotonID = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonID.src = cms.InputTag('patJetsAK5PFNoPhotonID')
patCSVMJetsAK5PFNoPhotonID.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVM))

patCSVMJetsAK5PFNoPhotonIDPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDPt30Eta24.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonID')
patCSVMJetsAK5PFNoPhotonIDPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsAK5PFNoPhotonIDPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDPt30Eta50.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonID')
patCSVMJetsAK5PFNoPhotonIDPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsAK5PFNoPhotonIDPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonID')
patCSVMJetsAK5PFNoPhotonIDPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

######
patCSVMJetsAK5PFNoPhotonIDPFIso = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDPFIso.src = cms.InputTag('patJetsAK5PFNoPhotonIDPFIso')
patCSVMJetsAK5PFNoPhotonIDPFIso.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVM))

patCSVMJetsAK5PFNoPhotonIDPFIsoPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDPFIsoPt30Eta24.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIso')
patCSVMJetsAK5PFNoPhotonIDPFIsoPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsAK5PFNoPhotonIDPFIsoPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDPFIsoPt30Eta50.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIso')
patCSVMJetsAK5PFNoPhotonIDPFIsoPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsAK5PFNoPhotonIDPFIsoPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDPFIsoPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIso')
patCSVMJetsAK5PFNoPhotonIDPFIsoPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

###patJetsPFNoPhoton
patCSVMJetsPFNoPhotonID = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonID.src = cms.InputTag('patJetsPFNoPhotonID')
patCSVMJetsPFNoPhotonID.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVM))

patCSVMJetsPFNoPhotonIDPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDPt30Eta24.src = cms.InputTag('patCSVMJetsPFNoPhotonID')
patCSVMJetsPFNoPhotonIDPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsPFNoPhotonIDPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDPt30Eta50.src = cms.InputTag('patCSVMJetsPFNoPhotonID')
patCSVMJetsPFNoPhotonIDPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsPFNoPhotonIDPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDPt50Eta25.src = cms.InputTag('patCSVMJetsPFNoPhotonID')
patCSVMJetsPFNoPhotonIDPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

######
patCSVMJetsPFNoPhotonIDPFIso = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDPFIso.src = cms.InputTag('patJetsPFNoPhotonIDPFIso')
patCSVMJetsPFNoPhotonIDPFIso.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVM))

patCSVMJetsPFNoPhotonIDPFIsoPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDPFIsoPt30Eta24.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIso')
patCSVMJetsPFNoPhotonIDPFIsoPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsPFNoPhotonIDPFIsoPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDPFIsoPt30Eta50.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIso')
patCSVMJetsPFNoPhotonIDPFIsoPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsPFNoPhotonIDPFIsoPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDPFIsoPt50Eta25.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIso')
patCSVMJetsPFNoPhotonIDPFIsoPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')


############Tight
patCSVTJetsAK5PFNoPhotonID = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonID.src = cms.InputTag('patJetsAK5PFNoPhotonID')
patCSVTJetsAK5PFNoPhotonID.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVT))

patCSVTJetsAK5PFNoPhotonIDPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDPt30Eta24.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonID')
patCSVTJetsAK5PFNoPhotonIDPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsAK5PFNoPhotonIDPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDPt30Eta50.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonID')
patCSVTJetsAK5PFNoPhotonIDPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsAK5PFNoPhotonIDPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDPt50Eta25.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonID')
patCSVTJetsAK5PFNoPhotonIDPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

######
patCSVTJetsAK5PFNoPhotonIDPFIso = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDPFIso.src = cms.InputTag('patJetsAK5PFNoPhotonIDPFIso')
patCSVTJetsAK5PFNoPhotonIDPFIso.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVT))

patCSVTJetsAK5PFNoPhotonIDPFIsoPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDPFIsoPt30Eta24.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonIDPFIso')
patCSVTJetsAK5PFNoPhotonIDPFIsoPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsAK5PFNoPhotonIDPFIsoPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDPFIsoPt30Eta50.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonIDPFIso')
patCSVTJetsAK5PFNoPhotonIDPFIsoPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsAK5PFNoPhotonIDPFIsoPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDPFIsoPt50Eta25.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonIDPFIso')
patCSVTJetsAK5PFNoPhotonIDPFIsoPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

###patJetsPFNoPhoton
patCSVTJetsPFNoPhotonID = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonID.src = cms.InputTag('patJetsPFNoPhotonID')
patCSVTJetsPFNoPhotonID.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVT))

patCSVTJetsPFNoPhotonIDPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDPt30Eta24.src = cms.InputTag('patCSVTJetsPFNoPhotonID')
patCSVTJetsPFNoPhotonIDPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsPFNoPhotonIDPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDPt30Eta50.src = cms.InputTag('patCSVTJetsPFNoPhotonID')
patCSVTJetsPFNoPhotonIDPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsPFNoPhotonIDPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDPt50Eta25.src = cms.InputTag('patCSVTJetsPFNoPhotonID')
patCSVTJetsPFNoPhotonIDPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

######
patCSVTJetsPFNoPhotonIDPFIso = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDPFIso.src = cms.InputTag('patJetsPFNoPhotonIDPFIso')
patCSVTJetsPFNoPhotonIDPFIso.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVT))

patCSVTJetsPFNoPhotonIDPFIsoPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDPFIsoPt30Eta24.src = cms.InputTag('patCSVTJetsPFNoPhotonIDPFIso')
patCSVTJetsPFNoPhotonIDPFIsoPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsPFNoPhotonIDPFIsoPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDPFIsoPt30Eta50.src = cms.InputTag('patCSVTJetsPFNoPhotonIDPFIso')
patCSVTJetsPFNoPhotonIDPFIsoPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsPFNoPhotonIDPFIsoPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDPFIsoPt50Eta25.src = cms.InputTag('patCSVTJetsPFNoPhotonIDPFIso')
patCSVTJetsPFNoPhotonIDPFIsoPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')


### count the b-jets
countCSVMJetsAK5PFNoPhotonID = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonID.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonID')
countCSVMJetsAK5PFNoPhotonID.minNumber = cms.uint32(1)

countCSVMJetsAK5PFNoPhotonIDPt30 = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonIDPt30.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPt30')
countCSVMJetsAK5PFNoPhotonIDPt30.minNumber = cms.uint32(1)

countCSVMJetsAK5PFNoPhotonIDPt50Eta25 = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonIDPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPt50Eta25')
countCSVMJetsAK5PFNoPhotonIDPt50Eta25.minNumber = cms.uint32(1)

countCSVMJetsAK5PFNoPhotonIDPFIso = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonIDPFIso.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIso')
countCSVMJetsAK5PFNoPhotonIDPFIso.minNumber = cms.uint32(1)

countCSVMJetsAK5PFNoPhotonIDPFIsoPt30 = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonIDPFIsoPt30.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIsoPt30')
countCSVMJetsAK5PFNoPhotonIDPFIsoPt30.minNumber = cms.uint32(1)

countCSVMJetsAK5PFNoPhotonIDPFIsoPt50Eta25 = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonIDPFIsoPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIsoPt50Eta25')
countCSVMJetsAK5PFNoPhotonIDPFIsoPt50Eta25.minNumber = cms.uint32(1)

###patJetsPFNoPhoton
countCSVMJetsPFNoPhotonID = countPatJets.clone()
countCSVMJetsPFNoPhotonID.src = cms.InputTag('patCSVMJetsPFNoPhotonID')
countCSVMJetsPFNoPhotonID.minNumber = cms.uint32(1)

countCSVMJetsPFNoPhotonIDPt30 = countPatJets.clone()
countCSVMJetsPFNoPhotonIDPt30.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPt30')
countCSVMJetsPFNoPhotonIDPt30.minNumber = cms.uint32(1)

countCSVMJetsPFNoPhotonIDPt50Eta25 = countPatJets.clone()
countCSVMJetsPFNoPhotonIDPt50Eta25.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPt50Eta25')
countCSVMJetsPFNoPhotonIDPt50Eta25.minNumber = cms.uint32(1)

countCSVMJetsPFNoPhotonIDPFIso = countPatJets.clone()
countCSVMJetsPFNoPhotonIDPFIso.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIso')
countCSVMJetsPFNoPhotonIDPFIso.minNumber = cms.uint32(1)

countCSVMJetsPFNoPhotonIDPFIsoPt30 = countPatJets.clone()
countCSVMJetsPFNoPhotonIDPFIsoPt30.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIsoPt30')
countCSVMJetsPFNoPhotonIDPFIsoPt30.minNumber = cms.uint32(1)

countCSVMJetsPFNoPhotonIDPFIsoPt50Eta25 = countPatJets.clone()
countCSVMJetsPFNoPhotonIDPFIsoPt50Eta25.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIsoPt50Eta25')
countCSVMJetsPFNoPhotonIDPFIsoPt50Eta25.minNumber = cms.uint32(1)


zinvBVetoNoPhotonID = cms.Sequence(
    ~countCSVMJetsAK5PFNoPhotonID
)
zinvBVetoNoPhotonIDPt30 = cms.Sequence(
    ~countCSVMJetsAK5PFNoPhotonIDPt30
)
zinvBVetoPtNoPhotonID50Eta25 = cms.Sequence(
    ~countCSVMJetsAK5PFNoPhotonIDPt50Eta25
)

### create the jet collections
zinvBJetsAK5PFNoPhotonID = cms.Sequence(
      patCSVMJetsAK5PFNoPhotonID
    * patCSVMJetsAK5PFNoPhotonIDPt30Eta24
    * patCSVTJetsAK5PFNoPhotonID
    * patCSVTJetsAK5PFNoPhotonIDPt30Eta24
    #* patCSVJetsAK5PFNoPhotonIDPt30Eta50
    #* patCSVJetsAK5PFNoPhotonIDPt50Eta25
)

zinvBJetsAK5PFNoPhotonIDPFIso = cms.Sequence(
      patCSVMJetsAK5PFNoPhotonIDPFIso
    * patCSVMJetsAK5PFNoPhotonIDPFIsoPt30Eta24
    * patCSVTJetsAK5PFNoPhotonIDPFIso
    * patCSVTJetsAK5PFNoPhotonIDPFIsoPt30Eta24
    #* patCSVJetsAK5PFNoPhotonIDPFIsoPt30Eta50
    #* patCSVJetsAK5PFNoPhotonIDPFIsoPt50Eta25
)

zinvBJetsPFNoPhotonID = cms.Sequence(
      patCSVMJetsPFNoPhotonID
    * patCSVMJetsPFNoPhotonIDPt30Eta24
    * patCSVTJetsPFNoPhotonID
    * patCSVTJetsPFNoPhotonIDPt30Eta24
    #* patCSVJetsPFNoPhotonIDPt30Eta50
    #* patCSVJetsPFNoPhotonIDPt50Eta25
)

zinvBJetsPFNoPhotonIDPFIso = cms.Sequence(
      patCSVMJetsPFNoPhotonIDPFIso
    * patCSVMJetsPFNoPhotonIDPFIsoPt30Eta24
    * patCSVTJetsPFNoPhotonIDPFIso
    * patCSVTJetsPFNoPhotonIDPFIsoPt30Eta24
    #* patCSVJetsPFNoPhotonIDPFIsoPt30Eta50
    #* patCSVJetsPFNoPhotonIDPFIsoPt50Eta25
)

#############
from ZInvisibleBkgds.Photons.specialObjectCleaner_cff import specialPhotonCleanedJets

patJetsAK5PFNoPhotonIDSpecial              = specialPhotonCleanedJets.clone()
patJetsAK5PFNoPhotonIDSpecial.jetLabel     = cms.InputTag('patJetsAK5PF')
patJetsAK5PFNoPhotonIDSpecial.objectLabel  = cms.InputTag('patPhotonsID')

patJetsPFNoPhotonIDSpecial     = patJetsAK5PFNoPhotonIDSpecial.clone()
patJetsPFNoPhotonIDSpecial.src = cms.InputTag('patJetsPF')

patJetsAK5PFNoPhotonIDPFIsoSpecial     = patJetsAK5PFNoPhotonIDSpecial.clone()
patJetsAK5PFNoPhotonIDPFIsoSpecial.objectLabel  = cms.InputTag('patPhotonsIDPFIso')
patJetsAK5PFNoPhotonIDPFIsoSpecial.src = cms.InputTag('patJetsAK5PF')

patJetsPFNoPhotonIDPFIsoSpecial     = patJetsAK5PFNoPhotonIDSpecial.clone()
patJetsPFNoPhotonIDPFIsoSpecial.src = cms.InputTag('patJetsPF')

#MHT Jets
patJetsAK5PFNoPhotonIDSpecialPt30              = patJetsAK5PFNoPhotonIDSpecial.clone()
patJetsAK5PFNoPhotonIDSpecialPt30.jetLabel     = cms.InputTag('patJetsAK5PFPt30')

patJetsPFNoPhotonIDSpecialPt30          = patJetsPFNoPhotonIDSpecial.clone()
patJetsPFNoPhotonIDSpecialPt30.jetLabel = cms.InputTag('patJetsPFchsPt30')

patJetsAK5PFNoPhotonIDPFIsoSpecialPt30              = patJetsAK5PFNoPhotonIDPFIsoSpecial.clone()
patJetsAK5PFNoPhotonIDPFIsoSpecialPt30.jetLabel     = cms.InputTag('patJetsAK5PFPt30')

patJetsPFNoPhotonIDPFIsoSpecialPt30          = patJetsPFNoPhotonIDPFIsoSpecial.clone()
patJetsPFNoPhotonIDPFIsoSpecialPt30.jetLabel = cms.InputTag('patJetsPFchsPt30')

#HT Jets
patJetsAK5PFNoPhotonIDSpecialPt50Eta25          = patJetsAK5PFNoPhotonIDSpecialPt30.clone()
patJetsAK5PFNoPhotonIDSpecialPt50Eta25.jetLabel = cms.InputTag('patJetsAK5PFPt50Eta25')

patJetsPFNoPhotonIDSpecialPt50Eta25          = patJetsPFNoPhotonIDSpecialPt30.clone()
patJetsPFNoPhotonIDSpecialPt50Eta25.jetLabel = cms.InputTag('patJetsPFchsPt50Eta25')

patJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25          = patJetsAK5PFNoPhotonIDPFIsoSpecialPt30.clone()
patJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25.jetLabel = cms.InputTag('patJetsAK5PFPt50Eta25')

patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25          = patJetsPFNoPhotonIDPFIsoSpecialPt30.clone()
patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25.jetLabel = cms.InputTag('patJetsPFchsPt50Eta25')

specialPhotonCleanedPFJetsAK5PF = cms.Sequence(
    patJetsAK5PFNoPhotonIDSpecial           *
    patJetsAK5PFNoPhotonIDSpecialPt30       *
    patJetsAK5PFNoPhotonIDSpecialPt50Eta25  *
    patJetsAK5PFNoPhotonIDPFIsoSpecial      *
    patJetsAK5PFNoPhotonIDPFIsoSpecialPt30  *
    patJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25
)
specialPhotonCleanedPFJetsPF = cms.Sequence(
    patJetsPFNoPhotonIDSpecial           *
    patJetsPFNoPhotonIDSpecialPt30       *
    patJetsPFNoPhotonIDSpecialPt50Eta25  *
    patJetsPFNoPhotonIDPFIsoSpecial      *
    patJetsPFNoPhotonIDPFIsoSpecialPt30  *
    patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25
)

#############

countJetsAK5PFNoPhotonSpecialPt50Eta25DiJets           = countPatJets.clone()
countJetsAK5PFNoPhotonSpecialPt50Eta25DiJets.src       = cms.InputTag('patJetsAK5PFNoPhotonSpecialPt50Eta25')
countJetsAK5PFNoPhotonSpecialPt50Eta25DiJets.minNumber = cms.uint32(2)

countJetsAK5PFNoPhotonSpecialPt50Eta25           = countPatJets.clone()
countJetsAK5PFNoPhotonSpecialPt50Eta25.src       = cms.InputTag('patJetsAK5PFNoPhotonSpecialPt50Eta25')
countJetsAK5PFNoPhotonSpecialPt50Eta25.minNumber = cms.uint32(3)

countJetsPFNoPhotonSpecialPt50Eta25DiJets           = countPatJets.clone()
countJetsPFNoPhotonSpecialPt50Eta25DiJets.src       = cms.InputTag('patJetsPFNoPhotonSpecialPt50Eta25')
countJetsPFNoPhotonSpecialPt50Eta25DiJets.minNumber = cms.uint32(2)

countJetsPFNoPhotonSpecialPt50Eta25           = countPatJets.clone()
countJetsPFNoPhotonSpecialPt50Eta25.src       = cms.InputTag('patJetsPFNoPhotonSpecialPt50Eta25')
countJetsPFNoPhotonSpecialPt50Eta25.minNumber = cms.uint32(3)


######
# b-tagged jets
#from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from ZInvisibleBkgds.Photons.specialJetSelector_cff import selectedRA2PatJets

###Medium
patCSVMJetsAK5PFNoPhotonIDSpecial = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDSpecial.src = cms.InputTag('patJetsAK5PFNoPhotonIDSpecial')
patCSVMJetsAK5PFNoPhotonIDSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVM))

patCSVMJetsAK5PFNoPhotonIDSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDSpecialPt30Eta24.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDSpecial')
patCSVMJetsAK5PFNoPhotonIDSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsAK5PFNoPhotonIDSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDSpecialPt30Eta50.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDSpecial')
patCSVMJetsAK5PFNoPhotonIDSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsAK5PFNoPhotonIDSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDSpecialPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDSpecial')
patCSVMJetsAK5PFNoPhotonIDSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMJetsAK5PFNoPhotonIDPFIsoSpecial = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDPFIsoSpecial.src = cms.InputTag('patJetsAK5PFNoPhotonIDPFIsoSpecial')
patCSVMJetsAK5PFNoPhotonIDPFIsoSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVM))

patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta24.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIsoSpecial')
patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta50.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIsoSpecial')
patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIsoSpecial')
patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

###patJetsPFNoPhotonIDSpecial
patCSVMJetsPFNoPhotonIDSpecial = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDSpecial.src = cms.InputTag('patJetsPFNoPhotonIDSpecial')
patCSVMJetsPFNoPhotonIDSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVM))

patCSVMJetsPFNoPhotonIDSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDSpecialPt30Eta24.src = cms.InputTag('patCSVMJetsPFNoPhotonIDSpecial')
patCSVMJetsPFNoPhotonIDSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsPFNoPhotonIDSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDSpecialPt30Eta50.src = cms.InputTag('patCSVMJetsPFNoPhotonIDSpecial')
patCSVMJetsPFNoPhotonIDSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsPFNoPhotonIDSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDSpecialPt50Eta25.src = cms.InputTag('patCSVMJetsPFNoPhotonIDSpecial')
patCSVMJetsPFNoPhotonIDSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMJetsPFNoPhotonIDPFIsoSpecial = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDPFIsoSpecial.src = cms.InputTag('patJetsPFNoPhotonIDPFIsoSpecial')
patCSVMJetsPFNoPhotonIDPFIsoSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVM))

patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta24.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIsoSpecial')
patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta50.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIsoSpecial')
patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsPFNoPhotonIDPFIsoSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsPFNoPhotonIDPFIsoSpecialPt50Eta25.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIsoSpecial')
patCSVMJetsPFNoPhotonIDPFIsoSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

###Tight
patCSVTJetsAK5PFNoPhotonIDSpecial = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDSpecial.src = cms.InputTag('patJetsAK5PFNoPhotonIDSpecial')
patCSVTJetsAK5PFNoPhotonIDSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVT))

patCSVTJetsAK5PFNoPhotonIDSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDSpecialPt30Eta24.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonIDSpecial')
patCSVTJetsAK5PFNoPhotonIDSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsAK5PFNoPhotonIDSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDSpecialPt30Eta50.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonIDSpecial')
patCSVTJetsAK5PFNoPhotonIDSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsAK5PFNoPhotonIDSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDSpecialPt50Eta25.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonIDSpecial')
patCSVTJetsAK5PFNoPhotonIDSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVTJetsAK5PFNoPhotonIDPFIsoSpecial = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDPFIsoSpecial.src = cms.InputTag('patJetsAK5PFNoPhotonIDPFIsoSpecial')
patCSVTJetsAK5PFNoPhotonIDPFIsoSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVT))

patCSVTJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta24.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonIDPFIsoSpecial')
patCSVTJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta50.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonIDPFIsoSpecial')
patCSVTJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25.src = cms.InputTag('patCSVTJetsAK5PFNoPhotonIDPFIsoSpecial')
patCSVTJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

###patJetsPFNoPhotonIDSpecial
patCSVTJetsPFNoPhotonIDSpecial = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDSpecial.src = cms.InputTag('patJetsPFNoPhotonIDSpecial')
patCSVTJetsPFNoPhotonIDSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVT))

patCSVTJetsPFNoPhotonIDSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDSpecialPt30Eta24.src = cms.InputTag('patCSVTJetsPFNoPhotonIDSpecial')
patCSVTJetsPFNoPhotonIDSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsPFNoPhotonIDSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDSpecialPt30Eta50.src = cms.InputTag('patCSVTJetsPFNoPhotonIDSpecial')
patCSVTJetsPFNoPhotonIDSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsPFNoPhotonIDSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDSpecialPt50Eta25.src = cms.InputTag('patCSVTJetsPFNoPhotonIDSpecial')
patCSVTJetsPFNoPhotonIDSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVTJetsPFNoPhotonIDPFIsoSpecial = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDPFIsoSpecial.src = cms.InputTag('patJetsPFNoPhotonIDPFIsoSpecial')
patCSVTJetsPFNoPhotonIDPFIsoSpecial.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVT))

patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta24.src = cms.InputTag('patCSVTJetsPFNoPhotonIDPFIsoSpecial')
patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta50.src = cms.InputTag('patCSVTJetsPFNoPhotonIDPFIsoSpecial')
patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsPFNoPhotonIDPFIsoSpecialPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsPFNoPhotonIDPFIsoSpecialPt50Eta25.src = cms.InputTag('patCSVTJetsPFNoPhotonIDPFIsoSpecial')
patCSVTJetsPFNoPhotonIDPFIsoSpecialPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

### count the b-jets
countCSVMJetsAK5PFNoPhotonIDSpecial = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonIDSpecial.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDSpecial')
countCSVMJetsAK5PFNoPhotonIDSpecial.minNumber = cms.uint32(1)

countCSVMJetsAK5PFNoPhotonIDSpecialPt30 = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonIDSpecialPt30.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDSpecialPt30')
countCSVMJetsAK5PFNoPhotonIDSpecialPt30.minNumber = cms.uint32(1)

countCSVMJetsAK5PFNoPhotonIDSpecialPt50Eta25 = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonIDSpecialPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDSpecialPt50Eta25')
countCSVMJetsAK5PFNoPhotonIDSpecialPt50Eta25.minNumber = cms.uint32(1)

countCSVMJetsAK5PFNoPhotonIDPFIso = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonIDPFIso.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIso')
countCSVMJetsAK5PFNoPhotonIDPFIso.minNumber = cms.uint32(1)

countCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt30 = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt30.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt30')
countCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt30.minNumber = cms.uint32(1)

countCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25 = countPatJets.clone()
countCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25')
countCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25.minNumber = cms.uint32(1)

###patJetsPFNoPhoton
countCSVMJetsPFNoPhotonIDSpecial = countPatJets.clone()
countCSVMJetsPFNoPhotonIDSpecial.src = cms.InputTag('patCSVMJetsPFNoPhotonIDSpecial')
countCSVMJetsPFNoPhotonIDSpecial.minNumber = cms.uint32(1)

countCSVMJetsPFNoPhotonIDSpecialPt30 = countPatJets.clone()
countCSVMJetsPFNoPhotonIDSpecialPt30.src = cms.InputTag('patCSVMJetsPFNoPhotonIDSpecialPt30')
countCSVMJetsPFNoPhotonIDSpecialPt30.minNumber = cms.uint32(1)

countCSVMJetsPFNoPhotonIDSpecialPt50Eta25 = countPatJets.clone()
countCSVMJetsPFNoPhotonIDSpecialPt50Eta25.src = cms.InputTag('patCSVMJetsPFNoPhotonIDSpecialPt50Eta25')
countCSVMJetsPFNoPhotonIDSpecialPt50Eta25.minNumber = cms.uint32(1)

countCSVMJetsPFNoPhotonIDPFIsoSpecial = countPatJets.clone()
countCSVMJetsPFNoPhotonIDPFIsoSpecial.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIsoSpecial')
countCSVMJetsPFNoPhotonIDPFIsoSpecial.minNumber = cms.uint32(1)

countCSVMJetsPFNoPhotonIDPFIsoSpecialPt30 = countPatJets.clone()
countCSVMJetsPFNoPhotonIDPFIsoSpecialPt30.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30')
countCSVMJetsPFNoPhotonIDPFIsoSpecialPt30.minNumber = cms.uint32(1)

countCSVMJetsPFNoPhotonIDPFIsoSpecialPt50Eta25 = countPatJets.clone()
countCSVMJetsPFNoPhotonIDPFIsoSpecialPt50Eta25.src = cms.InputTag('patCSVMJetsPFNoPhotonIDPFIsoSpecialPt50Eta25')
countCSVMJetsPFNoPhotonIDPFIsoSpecialPt50Eta25.minNumber = cms.uint32(1)


zinvBVetoNoPhotonIDSpecial = cms.Sequence(
    ~countCSVMJetsAK5PFNoPhotonIDSpecial
)
zinvBVetoNoPhotonIDSpecialPt30 = cms.Sequence(
    ~countCSVMJetsAK5PFNoPhotonIDSpecialPt30
)
zinvBVetoNoPhotonIDSpecialPt50Eta25 = cms.Sequence(
    ~countCSVMJetsAK5PFNoPhotonIDSpecialPt50Eta25
)

### create the jet collections
zinvBJetsAK5PFNoPhotonIDSpecial = cms.Sequence(
      patCSVMJetsAK5PFNoPhotonIDSpecial
    * patCSVMJetsAK5PFNoPhotonIDSpecialPt30Eta24
    * patCSVTJetsAK5PFNoPhotonIDSpecial
    * patCSVTJetsAK5PFNoPhotonIDSpecialPt30Eta24
    #* patCSVJetsAK5PFNoPhotonIDSpecialPt30Eta50
    #* patCSVJetsAK5PFNoPhotonIDSpecialPt50Eta25
)

zinvBJetsAK5PFNoPhotonIDPFIsoSpecial = cms.Sequence(
      patCSVMJetsAK5PFNoPhotonIDPFIsoSpecial
    * patCSVMJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta24
    * patCSVTJetsAK5PFNoPhotonIDPFIsoSpecial
    * patCSVTJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta24
    #* patCSVJetsAK5PFNoPhotonIDPFIsoSpecialPt30Eta50
    #* patCSVJetsAK5PFNoPhotonIDPFIsoSpecialPt50Eta25
)

zinvBJetsPFNoPhotonIDSpecial = cms.Sequence(
      patCSVMJetsPFNoPhotonIDSpecial
    * patCSVMJetsPFNoPhotonIDSpecialPt30Eta24
    * patCSVTJetsPFNoPhotonIDSpecial
    * patCSVTJetsPFNoPhotonIDSpecialPt30Eta24
    #* patCSVJetsPFNoPhotonIDSpecialPt30Eta50
    #* patCSVJetsPFNoPhotonIDSpecialPt50Eta25
)

zinvBJetsPFNoPhotonIDPFIsoSpecial = cms.Sequence(
      patCSVMJetsPFNoPhotonIDPFIsoSpecial
    * patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta24
    * patCSVTJetsPFNoPhotonIDPFIsoSpecial
    * patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta24
    #* patCSVJetsPFNoPhotonIDPFIsoSpecialPt30Eta50
    #* patCSVJetsPFNoPhotonIDPFIsoSpecialPt50Eta25
)
