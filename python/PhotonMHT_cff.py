
from SandBox.Skims.mhtProducer_cfi import *

# MHT using PFJets cleaned for electrons
mhtPFNoPhotID = mht.clone()
mhtPFNoPhotID.JetCollection = cms.InputTag("patJetsAK5PFNoPhotonIDPt30")

mhtPFchsNoPhotID = mht.clone()
mhtPFchsNoPhotID.JetCollection = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")

mhtPFNoPhotIDPFIso = mht.clone()
mhtPFNoPhotIDPFIso.JetCollection = cms.InputTag("patJetsAK5PFNoPhotonIDPFIsoPt30")

mhtPFchsNoPhotIDPFIso = mht.clone()
mhtPFchsNoPhotIDPFIso.JetCollection = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")

from SandBox.Skims.mhtFilter_cfi import *

# filter on PFMHT
mhtPFNoPhotFilter = mhtFilter.clone()
mhtPFNoPhotFilter.MHTSource = cms.InputTag("mhtPFNoPhotIDPFIso")

mhtPFchsNoPhotFilter = mhtFilter.clone()
mhtPFchsNoPhotFilter.MHTSource = cms.InputTag("mhtPFchsNoPhotIDPFIso")

