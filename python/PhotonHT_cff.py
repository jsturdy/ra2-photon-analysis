
from SandBox.Skims.htProducer_cfi import *

# HT using PFJets cleaned for electrons
htPFchsNoPhotID = ht.clone()
htPFchsNoPhotID.JetCollection = cms.InputTag("patJetsPFNoPhotonIDSpecialPt50Eta25")

htPFchsNoPhotIDPFIso = ht.clone()
htPFchsNoPhotIDPFIso.JetCollection = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25")

htPFchsNoPhotFitTemplate = ht.clone()
htPFchsNoPhotFitTemplate.JetCollection = cms.InputTag("patJetsPFNoPhotonFitTemplateSpecialPt50Eta25")

htPFchsNoPhotJetFake = ht.clone()
htPFchsNoPhotJetFake.JetCollection = cms.InputTag("patJetsPFNoPhotonJetFakeSpecialPt50Eta25")

from SandBox.Skims.htFilter_cfi import *

# filter on PFHT
htPFchsNoPhotFilter = htFilter.clone()
htPFchsNoPhotFilter.HTSource = cms.InputTag("htPFchsNoPhotIDPFIso")

