import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

patJetsAK5PFNoMuon = cleanPatJets.clone()
patJetsAK5PFNoMuon.src = cms.InputTag('patJetsAK5PF')
patJetsAK5PFNoMuon.checkOverlaps.muons.src               = cms.InputTag('specialMuonCollection')
patJetsAK5PFNoMuon.checkOverlaps.muons.algorithm         = cms.string('byDeltaR')
patJetsAK5PFNoMuon.checkOverlaps.muons.preselection      = cms.string('')
patJetsAK5PFNoMuon.checkOverlaps.muons.deltaR            = cms.double(0.1)
patJetsAK5PFNoMuon.checkOverlaps.muons.pairCut           = cms.string('')
patJetsAK5PFNoMuon.checkOverlaps.muons.requireNoOverlaps = cms.bool(True)
patJetsAK5PFNoMuon.checkOverlaps.taus.src                = cms.InputTag('selectedPatTausPF')
patJetsAK5PFNoMuon.checkOverlaps.electrons.src           = cms.InputTag('patElectrons')
patJetsAK5PFNoMuon.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectrons')
patJetsAK5PFNoMuon.checkOverlaps.photons.src             = cms.InputTag('patPhotonsAlt')

patJetsAK5PFNoMuonPt30     = patJetsAK5PFNoMuon.clone()
patJetsAK5PFNoMuonPt30.src = cms.InputTag('patJetsAK5PFPt30')

patJetsPFNoMuon     = patJetsAK5PFNoMuonPt30.clone()
patJetsPFNoMuon.src = cms.InputTag('patJetsPF')

patJetsPFNoMuonPt30     = patJetsAK5PFNoMuonPt30.clone()
patJetsPFNoMuonPt30.src = cms.InputTag('patJetsPFchsPt30')

#####
patJetsAK5PFNoMuonPt50Eta25     = patJetsAK5PFNoMuonPt30.clone()
patJetsAK5PFNoMuonPt50Eta25.src = cms.InputTag('patJetsAK5PFPt50Eta25')

patJetsPFNoMuonPt50Eta25     = patJetsPFNoMuonPt30.clone()
patJetsPFNoMuonPt50Eta25.src = cms.InputTag('patJetsPFchsPt50Eta25')


muonCleanedPFJetsAK5PF = cms.Sequence(
    patJetsAK5PFNoMuon               *
    patJetsAK5PFNoMuonPt30           *
    patJetsAK5PFNoMuonPt50Eta25      
)
muonCleanedPFJetsPF = cms.Sequence(
    patJetsPFNoMuon               *
    patJetsPFNoMuonPt30           *
    patJetsPFNoMuonPt50Eta25      
)


######
# b-tagged jets
#from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from ZInvisibleBkgds.Photons.specialJetSelector_cff import selectedRA2PatJets
#selectedBasicPatJets = cms.EDFilter("RA2BasicJetSelector",
#    src = cms.InputTag("patJets"),
#    cut = cms.string("")
#)
CSVL = 0.244
CSVM = 0.679
CSVT = 0.898 
csvPoint = CSVT

####Medium
patCSVMJetsAK5PFNoMuon = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoMuon.src = cms.InputTag('patJetsAK5PFNoMuon')
patCSVMJetsAK5PFNoMuon.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVM))

patCSVMJetsAK5PFNoMuonPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoMuonPt30Eta24.src = cms.InputTag('patCSVMJetsAK5PFNoMuon')
patCSVMJetsAK5PFNoMuonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsAK5PFNoMuonPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoMuonPt30Eta50.src = cms.InputTag('patCSVMJetsAK5PFNoMuon')
patCSVMJetsAK5PFNoMuonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsAK5PFNoMuonPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFNoMuonPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PFNoMuon')
patCSVMJetsAK5PFNoMuonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

####Tight
patCSVTJetsAK5PFNoMuon = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoMuon.src = cms.InputTag('patJetsAK5PFNoMuon')
patCSVTJetsAK5PFNoMuon.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVT))

patCSVTJetsAK5PFNoMuonPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoMuonPt30Eta24.src = cms.InputTag('patCSVTJetsAK5PFNoMuon')
patCSVTJetsAK5PFNoMuonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsAK5PFNoMuonPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoMuonPt30Eta50.src = cms.InputTag('patCSVTJetsAK5PFNoMuon')
patCSVTJetsAK5PFNoMuonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsAK5PFNoMuonPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFNoMuonPt50Eta25.src = cms.InputTag('patCSVTJetsAK5PFNoMuon')
patCSVTJetsAK5PFNoMuonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

####
patCSVMVAJetsAK5PFNoMuon = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoMuon.src = cms.InputTag('patJetsAK5PFNoMuon')
patCSVMVAJetsAK5PFNoMuon.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > %f'%(csvPoint))

patCSVMVAJetsAK5PFNoMuonPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoMuonPt30Eta24.src = cms.InputTag('patCSVMVAJetsAK5PFNoMuon')
patCSVMVAJetsAK5PFNoMuonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsAK5PFNoMuonPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoMuonPt30Eta50.src = cms.InputTag('patCSVMVAJetsAK5PFNoMuon')
patCSVMVAJetsAK5PFNoMuonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsAK5PFNoMuonPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFNoMuonPt50Eta25.src = cms.InputTag('patCSVMVAJetsAK5PFNoMuon')
patCSVMVAJetsAK5PFNoMuonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

###patJetsPFNoMuon
###Medium
patCSVMJetsPFNoMuon = selectedRA2PatJets.clone()
patCSVMJetsPFNoMuon.src = cms.InputTag('patJetsPFNoMuon')
patCSVMJetsPFNoMuon.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVM))

patCSVMJetsPFNoMuonPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsPFNoMuonPt30Eta24.src = cms.InputTag('patCSVMJetsPFNoMuon')
patCSVMJetsPFNoMuonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsPFNoMuonPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsPFNoMuonPt30Eta50.src = cms.InputTag('patCSVMJetsPFNoMuon')
patCSVMJetsPFNoMuonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsPFNoMuonPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsPFNoMuonPt50Eta25.src = cms.InputTag('patCSVMJetsPFNoMuon')
patCSVMJetsPFNoMuonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

###Tight
patCSVTJetsPFNoMuon = selectedRA2PatJets.clone()
patCSVTJetsPFNoMuon.src = cms.InputTag('patJetsPFNoMuon')
patCSVTJetsPFNoMuon.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVT))

patCSVTJetsPFNoMuonPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsPFNoMuonPt30Eta24.src = cms.InputTag('patCSVTJetsPFNoMuon')
patCSVTJetsPFNoMuonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsPFNoMuonPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsPFNoMuonPt30Eta50.src = cms.InputTag('patCSVTJetsPFNoMuon')
patCSVTJetsPFNoMuonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsPFNoMuonPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsPFNoMuonPt50Eta25.src = cms.InputTag('patCSVTJetsPFNoMuon')
patCSVTJetsPFNoMuonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

##
patCSVMVAJetsPFNoMuon = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoMuon.src = cms.InputTag('patJetsPFNoMuon')
patCSVMVAJetsPFNoMuon.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > %f'%(csvPoint))

patCSVMVAJetsPFNoMuonPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoMuonPt30Eta24.src = cms.InputTag('patCSVMVAJetsPFNoMuon')
patCSVMVAJetsPFNoMuonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsPFNoMuonPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoMuonPt30Eta50.src = cms.InputTag('patCSVMVAJetsPFNoMuon')
patCSVMVAJetsPFNoMuonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsPFNoMuonPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsPFNoMuonPt50Eta25.src = cms.InputTag('patCSVMVAJetsPFNoMuon')
patCSVMVAJetsPFNoMuonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

### count the b-jets
countCSVMJetsAK5PFNoMuonPt50Eta25 = countPatJets.clone()
countCSVMJetsAK5PFNoMuonPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PFNoMuonPt50Eta25')
countCSVMJetsAK5PFNoMuonPt50Eta25.minNumber = cms.uint32(1)

###patJetsPFNoMuon
countCSVMJetsPFNoMuonPt50Eta25 = countPatJets.clone()
countCSVMJetsPFNoMuonPt50Eta25.src = cms.InputTag('patCSVMJetsPFNoMuonPt50Eta25')
countCSVMJetsPFNoMuonPt50Eta25.minNumber = cms.uint32(1)


#zinvBVetoNoMuon = cms.Sequence(
#    ~countCSVMJetsAK5PFNoMuon
#)
#zinvBVetoNoMuonPt30 = cms.Sequence(
#    ~countCSVMJetsAK5PFNoMuonPt30
#)
#zinvBVetoNoMuonPt50Eta25 = cms.Sequence(
#    ~countCSVMJetsAK5PFNoMuonPt50Eta25
#)
### create the jet collections

zinvBJetsAK5PFNoMuon = cms.Sequence(
      patCSVMJetsAK5PFNoMuon
    * patCSVMJetsAK5PFNoMuonPt30Eta24
    * patCSVTJetsAK5PFNoMuon
    * patCSVTJetsAK5PFNoMuonPt30Eta24
    #* patCSVJetsAK5PFNoMuonPt30Eta50
    #* patCSVJetsAK5PFNoMuonPt50Eta25
    #* patCSVMVAJetsAK5PFNoMuon
    #* patCSVMVAJetsAK5PFNoMuonPt30Eta24
    #* patCSVMVAJetsAK5PFNoMuonPt30Eta50
    #* patCSVMVAJetsAK5PFNoMuonPt50Eta25
)

zinvBJetsPFNoMuon = cms.Sequence(
      patCSVMJetsPFNoMuon
    * patCSVMJetsPFNoMuonPt30Eta24
    * patCSVTJetsPFNoMuon
    * patCSVTJetsPFNoMuonPt30Eta24
    #* patCSVJetsPFNoMuonPt30Eta50
    #* patCSVJetsPFNoMuonPt50Eta25
    #* patCSVMVAJetsPFNoMuon
    #* patCSVMVAJetsPFNoMuonPt30Eta24
    #* patCSVMVAJetsPFNoMuonPt30Eta50
    #* patCSVMVAJetsPFNoMuonPt50Eta25
)
