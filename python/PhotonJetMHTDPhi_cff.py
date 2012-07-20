
from SandBox.Skims.jetMHTDPhiFilter_cfi import *

##### Photons
jetMHTPFNoPhotDPhiFilter = jetMHTDPhiFilter.clone()
jetMHTPFNoPhotDPhiFilter.JetSource = cms.InputTag("patJetsAK5PFPt30NoPhoton")
jetMHTPFNoPhotDPhiFilter.MHTSource = cms.InputTag("mhtPFNoPhot")

jetMHTPFNoPhotIDDPhiFilter = jetMHTDPhiFilter.clone()
jetMHTPFNoPhotIDDPhiFilter.JetSource = cms.InputTag("patJetsAK5PFPt30NoPhotonID")
jetMHTPFNoPhotIDDPhiFilter.MHTSource = cms.InputTag("mhtPFNoPhotID")

jetMHTPFNoPhotIDIsoDPhiFilter = jetMHTDPhiFilter.clone()
jetMHTPFNoPhotIDIsoDPhiFilter.JetSource = cms.InputTag("patJetsAK5PFPt30NoPhotonIDIso")
jetMHTPFNoPhotIDIsoDPhiFilter.MHTSource = cms.InputTag("mhtPFNoPhotIDIso")
