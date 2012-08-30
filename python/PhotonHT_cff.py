
from SandBox.Skims.htProducer_cfi import *

# HT using PFJets cleaned for electrons
htPFNoPhot = ht.clone()
htPFNoPhot.JetCollection = cms.InputTag("patJetsAK5PFPt50Eta25NoPhoton")

htPFNoPhotID = ht.clone()
htPFNoPhotID.JetCollection = cms.InputTag("patJetsAK5PFPt50Eta25NoPhotonID")

htPFNoPhotIDIso = ht.clone()
htPFNoPhotIDIso.JetCollection = cms.InputTag("patJetsAK5PFPt50Eta25NoPhotonIDIso")

htPFchsNoPhot = ht.clone()
htPFchsNoPhot.JetCollection = cms.InputTag("patJetsPFPt50Eta25NoPhoton")

htPFchsNoPhotID = ht.clone()
htPFchsNoPhotID.JetCollection = cms.InputTag("patJetsPFPt50Eta25NoPhotonID")

htPFchsNoPhotIDIso = ht.clone()
htPFchsNoPhotIDIso.JetCollection = cms.InputTag("patJetsPFPt50Eta25NoPhotonIDIso")

from SandBox.Skims.htFilter_cfi import *

# filter on PFHT
htPFNoPhotFilter = htFilter.clone()
htPFNoPhotFilter.HTSource = cms.InputTag("htPFNoPhot")

htPFNoPhotIDFilter = htFilter.clone()
htPFNoPhotIDFilter.HTSource = cms.InputTag("htPFNoPhotID")

htPFNoPhotIDIsoFilter = htFilter.clone()
htPFNoPhotIDIsoFilter.HTSource = cms.InputTag("htPFNoPhotIDIso")

htPFchsNoPhotFilter = htFilter.clone()
htPFchsNoPhotFilter.HTSource = cms.InputTag("htPFchsNoPhot")

htPFchsNoPhotIDFilter = htFilter.clone()
htPFchsNoPhotIDFilter.HTSource = cms.InputTag("htPFchsNoPhotID")

htPFchsNoPhotIDIsoFilter = htFilter.clone()
htPFchsNoPhotIDIsoFilter.HTSource = cms.InputTag("htPFchsNoPhotIDIso")

