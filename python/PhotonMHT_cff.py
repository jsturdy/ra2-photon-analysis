
from SandBox.Skims.mhtProducer_cfi import *

# MHT using PFJets cleaned for electrons
mhtPFNoPhot = mht.clone()
mhtPFNoPhot.JetCollection = cms.InputTag("patJetsAK5PFNoPhotonPt30")

mhtPFchsNoPhot = mht.clone()
mhtPFchsNoPhot.JetCollection = cms.InputTag("patJetsPFNoPhotonSpecialPt30")

from SandBox.Skims.mhtFilter_cfi import *

# filter on PFMHT
mhtPFNoPhotFilter = mhtFilter.clone()
mhtPFNoPhotFilter.MHTSource = cms.InputTag("mhtPFNoPhot")

mhtPFchsNoPhotFilter = mhtFilter.clone()
mhtPFchsNoPhotFilter.MHTSource = cms.InputTag("mhtPFchsNoPhot")

