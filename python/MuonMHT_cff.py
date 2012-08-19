
from SandBox.Skims.mhtProducer_cfi import *

# MHT using PFJets cleaned for electrons
mhtPFNoMuon = mht.clone()
mhtPFNoMuon.JetCollection = cms.InputTag("patJetsAK5PFNoMuonPt30")

mhtPFchsNoMuon = mht.clone()
mhtPFchsNoMuon.JetCollection = cms.InputTag("patJetsPFNoMuonPt30")

from SandBox.Skims.mhtFilter_cfi import *

# filter on PFMHT
mhtPFNoMuonFilter = mhtFilter.clone()
mhtPFNoMuonFilter.MHTSource = cms.InputTag("mhtPFNoMuon")

mhtPFchsNoMuonFilter = mhtFilter.clone()
mhtPFchsNoMuonFilter.MHTSource = cms.InputTag("mhtPFchsNoMuon")
