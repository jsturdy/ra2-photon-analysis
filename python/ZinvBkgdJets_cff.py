import FWCore.ParameterSet.Config as cms
CSVL = 0.244
CSVM = 0.679
CSVT = 0.898 
csvPoint = CSVT
# b-tagged jets
from PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi import *
from ZInvisibleBkgds.Photons.specialJetSelector_cff import selectedRA2PatJets

###patJetsPF
patCSVMJetsPF = selectedRA2PatJets.clone()
patCSVMJetsPF.src = cms.InputTag('newJetsMET')
patCSVMJetsPF.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVM))

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
patCSVTJetsPF.src = cms.InputTag('newJetsMET')
patCSVTJetsPF.cut = cms.string('bDiscriminator("combinedSecondaryVertexBJetTags") > %f'%(CSVT))

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
patCSVMVAJetsPF.src = cms.InputTag('newJetsMET')
patCSVMVAJetsPF.cut = cms.string('bDiscriminator("combinedSecondaryVertexMVABJetTags") > %f'%(csvPoint))

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
    ~countCSVMJetsPF
)
zinvBVetoPt30 = cms.Sequence(
    ~countCSVMJetsPFPt30
)
zinvBVetoPt50Eta25 = cms.Sequence(
    ~countCSVMJetsPFPt50Eta25
)

### create the jet collections
zinvBJetsPF = cms.Sequence(
      patCSVMJetsPF
    * patCSVMJetsPFPt30Eta24
    * patCSVTJetsPF
    * patCSVTJetsPFPt30Eta24
)
