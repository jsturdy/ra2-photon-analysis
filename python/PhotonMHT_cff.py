
from SandBox.Skims.mhtProducer_cfi import *

# MHT using PFJets cleaned for electrons
mhtPFchsNoPhotID = mht.clone()
mhtPFchsNoPhotID.JetCollection = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")

mhtPFchsNoPhotIDPFIso = mht.clone()
mhtPFchsNoPhotIDPFIso.JetCollection = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")

mhtPFchsNoPhotFitTemplate = mht.clone()
mhtPFchsNoPhotFitTemplate.JetCollection = cms.InputTag("patJetsPFNoPhotonFitTemplateSpecialPt30")

mhtPFchsNoPhotJetFake = mht.clone()
mhtPFchsNoPhotJetFake.JetCollection = cms.InputTag("patJetsPFNoPhotonJetFakeSpecialPt30")

from SandBox.Skims.mhtFilter_cfi import *

# filter on PFMHT
mhtPFchsNoPhotFilter = mhtFilter.clone()
mhtPFchsNoPhotFilter.MHTSource = cms.InputTag("mhtPFchsNoPhotIDPFIso")

