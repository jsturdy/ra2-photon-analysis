
from SandBox.Skims.htProducer_cfi import *

# HT using PFJets cleaned for electrons
htPFNoMuon = ht.clone()
htPFNoMuon.JetCollection = cms.InputTag("patJetsAK5PFNoMuonPt50Eta25")

htPFchsNoMuon = ht.clone()
htPFchsNoMuon.JetCollection = cms.InputTag("patJetsPFNoMuonPt50Eta25")

from SandBox.Skims.htFilter_cfi import *

# filter on PFHT
htPFNoMuonFilter = htFilter.clone()
htPFNoMuonFilter.HTSource = cms.InputTag("htPFNoMuon")

htPFchsNoMuonFilter = htFilter.clone()
htPFchsNoMuonFilter.HTSource = cms.InputTag("htPFchsNoMuon")
