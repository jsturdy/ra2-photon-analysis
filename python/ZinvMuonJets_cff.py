
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
patJetsAK5PFNoMuon.checkOverlaps.photons.src             = cms.InputTag('patPhotons')

patJetsAK5PFNoMuonPt30 = patJetsAK5PFNoMuon.clone()
patJetsAK5PFNoMuonPt30.src = cms.InputTag('patJetsAK5PFPt30')

patJetsPFNoMuon                         = patJetsAK5PFNoMuonPt30.clone()
patJetsPFNoMuon.src = cms.InputTag('patJetsPF')

patJetsPFNoMuonPt30                         = patJetsAK5PFNoMuonPt30.clone()
patJetsPFNoMuonPt30.src = cms.InputTag('patJetsPFPt30')

#####
patJetsAK5PFNoMuonPt50Eta25     = patJetsAK5PFNoMuonPt30.clone()
patJetsAK5PFNoMuonPt50Eta25.src = cms.InputTag('patJetsAK5PFPt50Eta25')

patJetsPFNoMuonPt50Eta25     = patJetsPFNoMuonPt30.clone()
patJetsPFNoMuonPt50Eta25.src = cms.InputTag('patJetsPFPt50Eta25')


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
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
patCSVJetsAK5PFNoMuon = selectedBasicPatJets.clone()
patCSVJetsAK5PFNoMuon.src = cms.InputTag('patJetsAK5PFNoMuon')
patCSVJetsAK5PFNoMuon.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898')

patCSVJetsAK5PFNoMuonPt30Eta24 = selectedBasicPatJets.clone()
patCSVJetsAK5PFNoMuonPt30Eta24.src = cms.InputTag('patCSVJetsAK5PFNoMuon')
patCSVJetsAK5PFNoMuonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVJetsAK5PFNoMuonPt30Eta50 = selectedBasicPatJets.clone()
patCSVJetsAK5PFNoMuonPt30Eta50.src = cms.InputTag('patCSVJetsAK5PFNoMuon')
patCSVJetsAK5PFNoMuonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVJetsAK5PFNoMuonPt50Eta25 = selectedBasicPatJets.clone()
patCSVJetsAK5PFNoMuonPt50Eta25.src = cms.InputTag('patCSVJetsAK5PFNoMuon')
patCSVJetsAK5PFNoMuonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMVAJetsAK5PFNoMuon = selectedBasicPatJets.clone()
patCSVMVAJetsAK5PFNoMuon.src = cms.InputTag('patJetsAK5PFNoMuon')
patCSVMVAJetsAK5PFNoMuon.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > 0.898')

patCSVMVAJetsAK5PFNoMuonPt30Eta24 = selectedBasicPatJets.clone()
patCSVMVAJetsAK5PFNoMuonPt30Eta24.src = cms.InputTag('patCSVMVAJetsAK5PFNoMuon')
patCSVMVAJetsAK5PFNoMuonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsAK5PFNoMuonPt30Eta50 = selectedBasicPatJets.clone()
patCSVMVAJetsAK5PFNoMuonPt30Eta50.src = cms.InputTag('patCSVMVAJetsAK5PFNoMuon')
patCSVMVAJetsAK5PFNoMuonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsAK5PFNoMuonPt50Eta25 = selectedBasicPatJets.clone()
patCSVMVAJetsAK5PFNoMuonPt50Eta25.src = cms.InputTag('patCSVMVAJetsAK5PFNoMuon')
patCSVMVAJetsAK5PFNoMuonPt50Eta25.cut = cms.string('pt > 30 && abs(eta) < 2.5')

###patJetsPFNoMuon
patCSVJetsPFNoMuon = selectedBasicPatJets.clone()
patCSVJetsPFNoMuon.src = cms.InputTag('patJetsPFNoMuon')
patCSVJetsPFNoMuon.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > 0.898')

patCSVJetsPFNoMuonPt30Eta24 = selectedBasicPatJets.clone()
patCSVJetsPFNoMuonPt30Eta24.src = cms.InputTag('patCSVJetsPFNoMuon')
patCSVJetsPFNoMuonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVJetsPFNoMuonPt30Eta50 = selectedBasicPatJets.clone()
patCSVJetsPFNoMuonPt30Eta50.src = cms.InputTag('patCSVJetsPFNoMuon')
patCSVJetsPFNoMuonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVJetsPFNoMuonPt50Eta25 = selectedBasicPatJets.clone()
patCSVJetsPFNoMuonPt50Eta25.src = cms.InputTag('patCSVJetsPFNoMuon')
patCSVJetsPFNoMuonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMVAJetsPFNoMuon = selectedBasicPatJets.clone()
patCSVMVAJetsPFNoMuon.src = cms.InputTag('patJetsPFNoMuon')
patCSVMVAJetsPFNoMuon.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > 0.898')

patCSVMVAJetsPFNoMuonPt30Eta24 = selectedBasicPatJets.clone()
patCSVMVAJetsPFNoMuonPt30Eta24.src = cms.InputTag('patCSVMVAJetsPFNoMuon')
patCSVMVAJetsPFNoMuonPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsPFNoMuonPt30Eta50 = selectedBasicPatJets.clone()
patCSVMVAJetsPFNoMuonPt30Eta50.src = cms.InputTag('patCSVMVAJetsPFNoMuon')
patCSVMVAJetsPFNoMuonPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsPFNoMuonPt50Eta25 = selectedBasicPatJets.clone()
patCSVMVAJetsPFNoMuonPt50Eta25.src = cms.InputTag('patCSVMVAJetsPFNoMuon')
patCSVMVAJetsPFNoMuonPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

### count the b-jets
countCSVJetsAK5PFNoMuonPt50Eta25 = countPatJets.clone()
countCSVJetsAK5PFNoMuonPt50Eta25.src = cms.InputTag('patCSVJetsAK5PFNoMuonPt50Eta25')
countCSVJetsAK5PFNoMuonPt50Eta25.minNumber = cms.uint32(1)

###patJetsPFNoMuon
countCSVJetsPFNoMuonPt50Eta25 = countPatJets.clone()
countCSVJetsPFNoMuonPt50Eta25.src = cms.InputTag('patCSVJetsPFNoMuonPt50Eta25')
countCSVJetsPFNoMuonPt50Eta25.minNumber = cms.uint32(1)


#zinvBVeto = cms.Sequence(
#    ~countSSVHEMBJetsAK5PFNoMuon
#)
#zinvBVetoPt30 = cms.Sequence(
#    ~countSSVHEMBJetsAK5PFNoMuonPt30
#)
#zinvBVetoPt50Eta25 = cms.Sequence(
#    ~countSSVHEMBJetsAK5PFNoMuonPt50Eta25
#)
### create the jet collections

zinvBJetsAK5PFNoMuon = cms.Sequence(
      patCSVJetsAK5PFNoMuon
    * patCSVJetsAK5PFNoMuonPt30Eta24
    #* patCSVJetsAK5PFNoMuonPt30Eta50
    #* patCSVJetsAK5PFNoMuonPt50Eta25
    #* patCSVMVAJetsAK5PFNoMuon
    #* patCSVMVAJetsAK5PFNoMuonPt30Eta24
    #* patCSVMVAJetsAK5PFNoMuonPt30Eta50
    #* patCSVMVAJetsAK5PFNoMuonPt50Eta25
)

zinvBJetsPFNoMuon = cms.Sequence(
      patCSVJetsPFNoMuon
    * patCSVJetsPFNoMuonPt30Eta24
    #* patCSVJetsPFNoMuonPt30Eta50
    #* patCSVJetsPFNoMuonPt50Eta25
    #* patCSVMVAJetsPFNoMuon
    #* patCSVMVAJetsPFNoMuonPt30Eta24
    #* patCSVMVAJetsPFNoMuonPt30Eta50
    #* patCSVMVAJetsPFNoMuonPt50Eta25
)
