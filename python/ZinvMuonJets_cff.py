import FWCore.ParameterSet.Config as cms

from PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *

patJetsPFNoMuon = cleanPatJets.clone()
patJetsPFNoMuon.src = cms.InputTag('newJetsMET')
patJetsPFNoMuon.checkOverlaps.muons.src               = cms.InputTag('specialMuonCollection')
patJetsPFNoMuon.checkOverlaps.muons.algorithm         = cms.string('byDeltaR')
patJetsPFNoMuon.checkOverlaps.muons.preselection      = cms.string('')
patJetsPFNoMuon.checkOverlaps.muons.deltaR            = cms.double(0.2) #changed from 0.1
patJetsPFNoMuon.checkOverlaps.muons.pairCut           = cms.string('')
patJetsPFNoMuon.checkOverlaps.muons.requireNoOverlaps = cms.bool(True)
patJetsPFNoMuon.checkOverlaps.taus.src                = cms.InputTag('selectedPatTausPF')
patJetsPFNoMuon.checkOverlaps.electrons.src           = cms.InputTag('patElectrons')
patJetsPFNoMuon.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectrons')
patJetsPFNoMuon.checkOverlaps.photons.src             = cms.InputTag('patPhotonsRA2')

patJetsPFNoMuonPt30     = patJetsPFNoMuon.clone()
patJetsPFNoMuonPt30.src = cms.InputTag('patJetsPFchsPt30')

patJetsPFNoMuonPt50Eta25     = patJetsPFNoMuonPt30.clone()
patJetsPFNoMuonPt50Eta25.src = cms.InputTag('patJetsPFchsPt50Eta25')


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
###patJetsPFNoMuon
countCSVMJetsPFNoMuonPt50Eta25 = countPatJets.clone()
countCSVMJetsPFNoMuonPt50Eta25.src = cms.InputTag('patCSVMJetsPFNoMuonPt50Eta25')
countCSVMJetsPFNoMuonPt50Eta25.minNumber = cms.uint32(1)


#zinvBVetoNoMuon = cms.Sequence(
#    ~countCSVMJetsPFNoMuon
#)
#zinvBVetoNoMuonPt30 = cms.Sequence(
#    ~countCSVMJetsPFNoMuonPt30
#)
#zinvBVetoNoMuonPt50Eta25 = cms.Sequence(
#    ~countCSVMJetsPFNoMuonPt50Eta25
#)
### create the jet collections

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
