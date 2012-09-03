
from SandBox.Skims.htProducer_cfi import *

# HT using PFJets cleaned for electrons
htPFNoPhot = ht.clone()
htPFNoPhot.JetCollection = cms.InputTag("patJetsAK5PFNoPhotonPt50Eta25")

htPFchsNoPhot = ht.clone()
htPFchsNoPhot.JetCollection = cms.InputTag("patJetsPFNoPhotonSpecialPt50Eta25")

from SandBox.Skims.htFilter_cfi import *

# filter on PFHT
htPFNoPhotFilter = htFilter.clone()
htPFNoPhotFilter.HTSource = cms.InputTag("htPFNoPhot")

htPFchsNoPhotFilter = htFilter.clone()
htPFchsNoPhotFilter.HTSource = cms.InputTag("htPFchsNoPhot")

