
from SandBox.Skims.htProducer_cfi import *

# HT using PFJets cleaned for electrons
htPFNoPhotID = ht.clone()
htPFNoPhotID.JetCollection = cms.InputTag("patJetsAK5PFNoPhotonIDPt50Eta25")

htPFchsNoPhotID = ht.clone()
htPFchsNoPhotID.JetCollection = cms.InputTag("patJetsPFNoPhotonIDSpecialPt50Eta25")

htPFNoPhotIDPFIso = ht.clone()
htPFNoPhotIDPFIso.JetCollection = cms.InputTag("patJetsAK5PFNoPhotonIDPFIsoPt50Eta25")

htPFchsNoPhotIDPFIso = ht.clone()
htPFchsNoPhotIDPFIso.JetCollection = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25")

from SandBox.Skims.htFilter_cfi import *

# filter on PFHT
htPFNoPhotFilter = htFilter.clone()
htPFNoPhotFilter.HTSource = cms.InputTag("htPFNoPhotIDPFIso")

htPFchsNoPhotFilter = htFilter.clone()
htPFchsNoPhotFilter.HTSource = cms.InputTag("htPFchsNoPhotIDPFIso")

