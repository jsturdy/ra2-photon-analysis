
from SandBox.Skims.jetMHTDPhiFilter_cfi import *

##### Photons
jetMHTPFNoPhotIDDPhiFilter = jetMHTDPhiFilter.clone()
jetMHTPFNoPhotIDDPhiFilter.JetSource = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
jetMHTPFNoPhotIDDPhiFilter.MHTSource = cms.InputTag("mhtPFchsNoPhotID")

jetMHTPFNoPhotIDPFIsoDPhiFilter = jetMHTDPhiFilter.clone()
jetMHTPFNoPhotIDPFIsoDPhiFilter.JetSource = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
jetMHTPFNoPhotIDPFIsoDPhiFilter.MHTSource = cms.InputTag("mhtPFchsNoPhotIDPFIso")
