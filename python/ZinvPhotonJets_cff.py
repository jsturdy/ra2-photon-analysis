import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

patJetsPFNoPhotonID = cleanPatJets.clone()
patJetsPFNoPhotonID.src = cms.InputTag('newJetsMET')
patJetsPFNoPhotonID.checkOverlaps.photons.src               = cms.InputTag('patPhotonsID')
patJetsPFNoPhotonID.checkOverlaps.photons.algorithm         = cms.string('byDeltaR')
patJetsPFNoPhotonID.checkOverlaps.photons.preselection      = cms.string('')
patJetsPFNoPhotonID.checkOverlaps.photons.deltaR            = cms.double(0.2) #changed from 0.1
patJetsPFNoPhotonID.checkOverlaps.photons.pairCut           = cms.string('')
patJetsPFNoPhotonID.checkOverlaps.photons.requireNoOverlaps = cms.bool(True)
patJetsPFNoPhotonID.checkOverlaps.taus.src                = cms.InputTag('selectedPatTausPF')
patJetsPFNoPhotonID.checkOverlaps.electrons.src           = cms.InputTag('patElectrons')
patJetsPFNoPhotonID.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectrons')
patJetsPFNoPhotonID.checkOverlaps.muons.src               = cms.InputTag('patMuonsPF')

patJetsPFNoPhotonIDPt30 = patJetsPFNoPhotonID.clone()
patJetsPFNoPhotonIDPt30.src = cms.InputTag('patJetsPFchsPt30')

patJetsPFNoPhotonIDPt50Eta25     = patJetsPFNoPhotonIDPt30.clone()
patJetsPFNoPhotonIDPt50Eta25.src = cms.InputTag('patJetsPFchsPt50Eta25')

####ID/PF Iso
patJetsPFNoPhotonIDPFIso                           = patJetsPFNoPhotonID.clone()
patJetsPFNoPhotonIDPFIso.checkOverlaps.photons.src = cms.InputTag('patPhotonsIDPFIso')
patJetsPFNoPhotonIDPFIsoPt30                       = patJetsPFNoPhotonIDPFIso.clone()
patJetsPFNoPhotonIDPFIsoPt30.src                   = cms.InputTag('patJetsPFchsPt30')
patJetsPFNoPhotonIDPFIsoPt50Eta25                  = patJetsPFNoPhotonIDPFIsoPt30.clone()
patJetsPFNoPhotonIDPFIsoPt50Eta25.src              = cms.InputTag('patJetsPFchsPt50Eta25')

###Photons which are Jet-fake candidates
patJetsPFNoJetFakePhoton                           = patJetsPFNoPhotonIDPFIso.clone()
patJetsPFNoJetFakePhoton.checkOverlaps.photons.src = cms.InputTag('patJetFakePhotons')
patJetsPFNoJetFakePhotonPt30                       = patJetsPFNoJetFakePhoton.clone()
patJetsPFNoJetFakePhotonPt30.src                   = cms.InputTag('patJetsPFchsPt30')
patJetsPFNoJetFakePhotonPt50Eta25                  = patJetsPFNoJetFakePhoton.clone()
patJetsPFNoJetFakePhotonPt50Eta25.src              = cms.InputTag('patJetsPFchsPt50Eta25')

###Photons which are Jet-fake candidates
patJetsPFNoFitTemplatePhoton                           = patJetsPFNoPhotonIDPFIso.clone()
patJetsPFNoFitTemplatePhoton.checkOverlaps.photons.src = cms.InputTag('patFitTemplatePhotons')
patJetsPFNoFitTemplatePhotonPt30                       = patJetsPFNoFitTemplatePhoton.clone()
patJetsPFNoFitTemplatePhotonPt30.src                   = cms.InputTag('patJetsPFchsPt30')
patJetsPFNoFitTemplatePhotonPt50Eta25                  = patJetsPFNoFitTemplatePhoton.clone()
patJetsPFNoFitTemplatePhotonPt50Eta25.src              = cms.InputTag('patJetsPFchsPt50Eta25')

#####

#####
photonCleanedPFJetsPF = cms.Sequence(
    #  patJetsPFNoPhotonID               
    #* patJetsPFNoPhotonIDPt30           
    #* patJetsPFNoPhotonIDPt50Eta25      
      patJetsPFNoPhotonIDPFIso          
    * patJetsPFNoPhotonIDPFIsoPt30      
    * patJetsPFNoPhotonIDPFIsoPt50Eta25
)
photonTemplateCleanedPFJetsPF = cms.Sequence(
      patJetsPFNoFitTemplatePhoton          
    * patJetsPFNoFitTemplatePhotonPt30      
    * patJetsPFNoFitTemplatePhotonPt50Eta25
)
photonJetFakeCleanedPFJetsPF = cms.Sequence(
      patJetsPFNoJetFakePhoton          
    * patJetsPFNoJetFakePhotonPt30      
    * patJetsPFNoJetFakePhotonPt50Eta25
)


#############

countJetsPFNoPhotonIDPFIsoPt50Eta25DiJets           = countPatJets.clone()
countJetsPFNoPhotonIDPFIsoPt50Eta25DiJets.src       = cms.InputTag('patJetsPFNoPhotonIDPFIsoPt50Eta25')
countJetsPFNoPhotonIDPFIsoPt50Eta25DiJets.minNumber = cms.uint32(2)

countJetsPFNoPhotonIDPFIsoPt50Eta25           = countPatJets.clone()
countJetsPFNoPhotonIDPFIsoPt50Eta25.src       = cms.InputTag('patJetsPFNoPhotonIDPFIsoPt50Eta25')
countJetsPFNoPhotonIDPFIsoPt50Eta25.minNumber = cms.uint32(3)


######
# b-tagged jets
from SandBox.Skims.basicJetSelector_cfi import *
#from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from ZInvisibleBkgds.Photons.specialJetSelector_cff import selectedRA2PatJets
CSVL = 0.244
CSVM = 0.679
CSVT = 0.898 
csvPoint = CSVT

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
    ~countCSVMJetsPFNoPhotonID
)
zinvBVetoNoPhotonIDPt30 = cms.Sequence(
    ~countCSVMJetsPFNoPhotonIDPt30
)
zinvBVetoPtNoPhotonID50Eta25 = cms.Sequence(
    ~countCSVMJetsPFNoPhotonIDPt50Eta25
)

### create the jet collections
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

patJetsPFNoPhotonIDSpecial              = specialPhotonCleanedJets.clone()
patJetsPFNoPhotonIDSpecial.jetLabel     = cms.InputTag('newJetsMET')
patJetsPFNoPhotonIDSpecial.objectLabel  = cms.InputTag('patPhotonsID')

patJetsPFNoPhotonIDPFIsoSpecial              = patJetsPFNoPhotonIDSpecial.clone()
patJetsPFNoPhotonIDPFIsoSpecial.objectLabel  = cms.InputTag('patPhotonsIDPFIso')

patJetsPFNoPhotonFitTemplateSpecial              = patJetsPFNoPhotonIDSpecial.clone()
patJetsPFNoPhotonFitTemplateSpecial.objectLabel  = cms.InputTag('patFitTemplatePhotons')

patJetsPFNoPhotonJetFakeSpecial              = patJetsPFNoPhotonIDSpecial.clone()
patJetsPFNoPhotonJetFakeSpecial.objectLabel  = cms.InputTag('patJetFakePhotons')


patJetsPFNoPhotonIDSpecialPt30          = patJetsPFNoPhotonIDSpecial.clone()
patJetsPFNoPhotonIDSpecialPt30.jetLabel = cms.InputTag('patJetsPFchsPt30')

patJetsPFNoPhotonIDPFIsoSpecialPt30          = patJetsPFNoPhotonIDPFIsoSpecial.clone()
patJetsPFNoPhotonIDPFIsoSpecialPt30.jetLabel = cms.InputTag('patJetsPFchsPt30')

patJetsPFNoPhotonFitTemplateSpecialPt30          = patJetsPFNoPhotonFitTemplateSpecial.clone()
patJetsPFNoPhotonFitTemplateSpecialPt30.jetLabel = cms.InputTag('patJetsPFchsPt30')

patJetsPFNoPhotonJetFakeSpecialPt30          = patJetsPFNoPhotonJetFakeSpecial.clone()
patJetsPFNoPhotonJetFakeSpecialPt30.jetLabel = cms.InputTag('patJetsPFchsPt30')

patJetsPFNoPhotonIDSpecialPt50Eta25          = patJetsPFNoPhotonIDSpecialPt30.clone()
patJetsPFNoPhotonIDSpecialPt50Eta25.jetLabel = cms.InputTag('patJetsPFchsPt50Eta25')

patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25          = patJetsPFNoPhotonIDPFIsoSpecialPt30.clone()
patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25.jetLabel = cms.InputTag('patJetsPFchsPt50Eta25')

patJetsPFNoPhotonFitTemplateSpecialPt50Eta25          = patJetsPFNoPhotonFitTemplateSpecialPt30.clone()
patJetsPFNoPhotonFitTemplateSpecialPt50Eta25.jetLabel = cms.InputTag('patJetsPFchsPt50Eta25')

patJetsPFNoPhotonJetFakeSpecialPt50Eta25          = patJetsPFNoPhotonJetFakeSpecialPt30.clone()
patJetsPFNoPhotonJetFakeSpecialPt50Eta25.jetLabel = cms.InputTag('patJetsPFchsPt50Eta25')

specialPhotonCleanedPFJetsPF = cms.Sequence(
      patJetsPFNoPhotonIDSpecial           
    * patJetsPFNoPhotonIDSpecialPt30       
    * patJetsPFNoPhotonIDSpecialPt50Eta25  
    * patJetsPFNoPhotonIDPFIsoSpecial      
    * patJetsPFNoPhotonIDPFIsoSpecialPt30  
    * patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25
)
specialFitTemplatePhotonCleanedPFJetsPF = cms.Sequence(
      patJetsPFNoPhotonFitTemplateSpecial      
    * patJetsPFNoPhotonFitTemplateSpecialPt30  
    * patJetsPFNoPhotonFitTemplateSpecialPt50Eta25
)
specialJetFakePhotonCleanedPFJetsPF = cms.Sequence(
      patJetsPFNoPhotonJetFakeSpecial      
    * patJetsPFNoPhotonJetFakeSpecialPt30  
    * patJetsPFNoPhotonJetFakeSpecialPt50Eta25
)

#############
countJetsPFNoPhotonSpecialPt50Eta25DiJets           = countPatJets.clone()
countJetsPFNoPhotonSpecialPt50Eta25DiJets.src       = cms.InputTag('patJetsPFNoPhotonSpecialPt50Eta25')
countJetsPFNoPhotonSpecialPt50Eta25DiJets.minNumber = cms.uint32(2)

countJetsPFNoPhotonSpecialPt50Eta25           = countPatJets.clone()
countJetsPFNoPhotonSpecialPt50Eta25.src       = cms.InputTag('patJetsPFNoPhotonSpecialPt50Eta25')
countJetsPFNoPhotonSpecialPt50Eta25.minNumber = cms.uint32(3)


######
# b-tagged jets
from SandBox.Skims.basicJetSelector_cfi import *
#from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from ZInvisibleBkgds.Photons.specialJetSelector_cff import selectedRA2PatJets

###Medium
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
    ~countCSVMJetsPFNoPhotonIDSpecial
)
zinvBVetoNoPhotonIDSpecialPt30 = cms.Sequence(
    ~countCSVMJetsPFNoPhotonIDSpecialPt30
)
zinvBVetoNoPhotonIDSpecialPt50Eta25 = cms.Sequence(
    ~countCSVMJetsPFNoPhotonIDSpecialPt50Eta25
)

### create the jet collections
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
