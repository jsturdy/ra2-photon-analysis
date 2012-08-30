
from SandBox.Skims.mhtProducer_cfi import *

# MHT using PFJets cleaned for electrons
mhtPFNoPhot = mht.clone()
mhtPFNoPhot.JetCollection = cms.InputTag("patJetsAK5PFPt30NoPhoton")

mhtPFNoPhotID = mht.clone()
mhtPFNoPhotID.JetCollection = cms.InputTag("patJetsAK5PFPt30NoPhotonID")

mhtPFNoPhotIDIso = mht.clone()
mhtPFNoPhotIDIso.JetCollection = cms.InputTag("patJetsAK5PFPt30NoPhotonIDIso")

mhtPFchsNoPhot = mht.clone()
mhtPFchsNoPhot.JetCollection = cms.InputTag("patJetsPFPt30NoPhoton")

mhtPFchsNoPhotID = mht.clone()
mhtPFchsNoPhotID.JetCollection = cms.InputTag("patJetsPFPt30NoPhotonID")

mhtPFchsNoPhotIDIso = mht.clone()
mhtPFchsNoPhotIDIso.JetCollection = cms.InputTag("patJetsPFPt30NoPhotonIDIso")

from SandBox.Skims.mhtFilter_cfi import *

# filter on PFMHT
mhtPFNoPhotFilter = mhtFilter.clone()
mhtPFNoPhotFilter.MHTSource = cms.InputTag("mhtPFNoPhot")

mhtPFNoPhotIDFilter = mhtFilter.clone()
mhtPFNoPhotIDFilter.MHTSource = cms.InputTag("mhtPFNoPhotID")

mhtPFNoPhotIDIsoFilter = mhtFilter.clone()
mhtPFNoPhotIDIsoFilter.MHTSource = cms.InputTag("mhtPFNoPhotIDIso")

mhtPFchsNoPhotFilter = mhtFilter.clone()
mhtPFchsNoPhotFilter.MHTSource = cms.InputTag("mhtPFchsNoPhot")

mhtPFchsNoPhotIDFilter = mhtFilter.clone()
mhtPFchsNoPhotIDFilter.MHTSource = cms.InputTag("mhtPFchsNoPhotID")

mhtPFchsNoPhotIDIsoFilter = mhtFilter.clone()
mhtPFchsNoPhotIDIsoFilter.MHTSource = cms.InputTag("mhtPFchsNoPhotIDIso")
