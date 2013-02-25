
from SandBox.Skims.htProducer_cfi import *

# HT using PFJets cleaned for electrons
htPFchsNoMuon = ht.clone()
htPFchsNoMuon.JetCollection = cms.InputTag("patJetsPFNoMuonPt50Eta25")

from SandBox.Skims.htFilter_cfi import *

# filter on PFHT
htPFchsNoMuonFilter = htFilter.clone()
htPFchsNoMuonFilter.HTSource = cms.InputTag("htPFchsNoMuon")
