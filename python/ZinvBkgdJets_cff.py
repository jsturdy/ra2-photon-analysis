import FWCore.ParameterSet.Config as cms
CSVL = 0.244
CSVM = 0.679
CSVT = 0.898 
csvPoint = CSVT
# b-tagged jets
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from ZInvisibleBkgds.Photons.specialJetSelector_cff import selectedRA2PatJets

patCSVMJetsAK5PF = selectedRA2PatJets.clone()
patCSVMJetsAK5PF.src = cms.InputTag('patJetsAK5PF')
patCSVMJetsAK5PF.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %d'%(CSVM))

patCSVMJetsAK5PFPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFPt30Eta24.src = cms.InputTag('patCSVMJetsAK5PF')
patCSVMJetsAK5PFPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsAK5PFPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFPt30Eta50.src = cms.InputTag('patCSVMJetsAK5PF')
patCSVMJetsAK5PFPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsAK5PFPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsAK5PFPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PF')
patCSVMJetsAK5PFPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVTJetsAK5PF = selectedRA2PatJets.clone()
patCSVTJetsAK5PF.src = cms.InputTag('patJetsAK5PF')
patCSVTJetsAK5PF.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %d'%(CSVT))

patCSVTJetsAK5PFPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFPt30Eta24.src = cms.InputTag('patCSVTJetsAK5PF')
patCSVTJetsAK5PFPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsAK5PFPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFPt30Eta50.src = cms.InputTag('patCSVTJetsAK5PF')
patCSVTJetsAK5PFPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsAK5PFPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsAK5PFPt50Eta25.src = cms.InputTag('patCSVTJetsAK5PF')
patCSVTJetsAK5PFPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMVAJetsAK5PF = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PF.src = cms.InputTag('patJetsAK5PF')
patCSVMVAJetsAK5PF.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > %d'%(csvPoint))

patCSVMVAJetsAK5PFPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFPt30Eta24.src = cms.InputTag('patCSVMVAJetsAK5PF')
patCSVMVAJetsAK5PFPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsAK5PFPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFPt30Eta50.src = cms.InputTag('patCSVMVAJetsAK5PF')
patCSVMVAJetsAK5PFPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsAK5PFPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsAK5PFPt50Eta25.src = cms.InputTag('patCSVMVAJetsAK5PF')
patCSVMVAJetsAK5PFPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

###patJetsPF
patCSVMJetsPF = selectedRA2PatJets.clone()
patCSVMJetsPF.src = cms.InputTag('patJetsPF')
patCSVMJetsPF.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %d'%(CSVM))

patCSVMJetsPFPt30Eta24 = selectedRA2PatJets.clone()
patCSVMJetsPFPt30Eta24.src = cms.InputTag('patCSVMJetsPF')
patCSVMJetsPFPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMJetsPFPt30Eta50 = selectedRA2PatJets.clone()
patCSVMJetsPFPt30Eta50.src = cms.InputTag('patCSVMJetsPF')
patCSVMJetsPFPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMJetsPFPt50Eta25 = selectedRA2PatJets.clone()
patCSVMJetsPFPt50Eta25.src = cms.InputTag('patCSVMJetsPF')
patCSVMJetsPFPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVTJetsPF = selectedRA2PatJets.clone()
patCSVTJetsPF.src = cms.InputTag('patJetsPF')
patCSVTJetsPF.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %d'%(CSVT))

patCSVTJetsPFPt30Eta24 = selectedRA2PatJets.clone()
patCSVTJetsPFPt30Eta24.src = cms.InputTag('patCSVTJetsPF')
patCSVTJetsPFPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVTJetsPFPt30Eta50 = selectedRA2PatJets.clone()
patCSVTJetsPFPt30Eta50.src = cms.InputTag('patCSVTJetsPF')
patCSVTJetsPFPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVTJetsPFPt50Eta25 = selectedRA2PatJets.clone()
patCSVTJetsPFPt50Eta25.src = cms.InputTag('patCSVTJetsPF')
patCSVTJetsPFPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

patCSVMVAJetsPF = selectedRA2PatJets.clone()
patCSVMVAJetsPF.src = cms.InputTag('patJetsPF')
patCSVMVAJetsPF.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > %d'%(csvPoint))

patCSVMVAJetsPFPt30Eta24 = selectedRA2PatJets.clone()
patCSVMVAJetsPFPt30Eta24.src = cms.InputTag('patCSVMVAJetsPF')
patCSVMVAJetsPFPt30Eta24.cut = cms.string('pt > 30 && abs(eta) < 2.4')

patCSVMVAJetsPFPt30Eta50 = selectedRA2PatJets.clone()
patCSVMVAJetsPFPt30Eta50.src = cms.InputTag('patCSVMVAJetsPF')
patCSVMVAJetsPFPt30Eta50.cut = cms.string('pt > 30 && abs(eta) < 5.0')

patCSVMVAJetsPFPt50Eta25 = selectedRA2PatJets.clone()
patCSVMVAJetsPFPt50Eta25.src = cms.InputTag('patCSVMVAJetsPF')
patCSVMVAJetsPFPt50Eta25.cut = cms.string('pt > 50 && abs(eta) < 2.5')

from PhysicsTools.PatAlgos.selectionLayer1.jetCountFilter_cfi import *
### count the b-jets
countCSVMJetsAK5PF = countPatJets.clone()
countCSVMJetsAK5PF.src = cms.InputTag('patCSVMJetsAK5PF')
countCSVMJetsAK5PF.minNumber = cms.uint32(1)

countCSVTJetsAK5PF = countPatJets.clone()
countCSVTJetsAK5PF.src = cms.InputTag('patCSVTJetsAK5PF')
countCSVTJetsAK5PF.minNumber = cms.uint32(1)

countCSVMJetsAK5PFPt30 = countPatJets.clone()
countCSVMJetsAK5PFPt30.src = cms.InputTag('patCSVMJetsAK5PFPt30')
countCSVMJetsAK5PFPt30.minNumber = cms.uint32(1)

countCSVTJetsAK5PFPt30 = countPatJets.clone()
countCSVTJetsAK5PFPt30.src = cms.InputTag('patCSVTJetsAK5PFPt30')
countCSVTJetsAK5PFPt30.minNumber = cms.uint32(1)

countCSVMJetsAK5PFPt50Eta25 = countPatJets.clone()
countCSVMJetsAK5PFPt50Eta25.src = cms.InputTag('patCSVMJetsAK5PFPt50Eta25')
countCSVMJetsAK5PFPt50Eta25.minNumber = cms.uint32(1)

countCSVTJetsAK5PFPt50Eta25 = countPatJets.clone()
countCSVTJetsAK5PFPt50Eta25.src = cms.InputTag('patCSVTJetsAK5PFPt50Eta25')
countCSVTJetsAK5PFPt50Eta25.minNumber = cms.uint32(1)

###patJetsPF
countCSVMJetsPF = countPatJets.clone()
countCSVMJetsPF.src = cms.InputTag('patCSVMJetsPF')
countCSVMJetsPF.minNumber = cms.uint32(1)

countCSVTJetsPF = countPatJets.clone()
countCSVTJetsPF.src = cms.InputTag('patCSVTJetsPF')
countCSVTJetsPF.minNumber = cms.uint32(1)

countCSVMJetsPFPt30 = countPatJets.clone()
countCSVMJetsPFPt30.src = cms.InputTag('patCSVMJetsPFPt30')
countCSVMJetsPFPt30.minNumber = cms.uint32(1)

countCSVTJetsPFPt30 = countPatJets.clone()
countCSVTJetsPFPt30.src = cms.InputTag('patCSVTJetsPFPt30')
countCSVTJetsPFPt30.minNumber = cms.uint32(1)

countCSVMJetsPFPt50Eta25 = countPatJets.clone()
countCSVMJetsPFPt50Eta25.src = cms.InputTag('patCSVMJetsPFPt50Eta25')
countCSVMJetsPFPt50Eta25.minNumber = cms.uint32(1)

countCSVTJetsPFPt50Eta25 = countPatJets.clone()
countCSVTJetsPFPt50Eta25.src = cms.InputTag('patCSVTJetsPFPt50Eta25')
countCSVTJetsPFPt50Eta25.minNumber = cms.uint32(1)


zinvBVeto = cms.Sequence(
    ~countCSVMJetsAK5PF
)
zinvBVetoPt30 = cms.Sequence(
    ~countCSVMJetsAK5PFPt30
)
zinvBVetoPt50Eta25 = cms.Sequence(
    ~countCSVMJetsAK5PFPt50Eta25
)

### create the jet collections
zinvBJetsAK5PF = cms.Sequence(
      patCSVMJetsAK5PF
    * patCSVMJetsAK5PFPt30Eta24
    * patCSVTJetsAK5PF
    * patCSVTJetsAK5PFPt30Eta24
)

zinvBJetsPF = cms.Sequence(
      patCSVMJetsPF
    * patCSVMJetsPFPt30Eta24
    * patCSVTJetsPF
    * patCSVTJetsPFPt30Eta24
)
