
from SandBox.Skims.jetMHTDPhiFilter_cfi import *

##### Muons
jetMHTPFNoMuonDPhiFilter = jetMHTDPhiFilter.clone()
jetMHTPFNoMuonDPhiFilter.JetSource = cms.InputTag("patJetsPFPt30NoMuon")
jetMHTPFNoMuonDPhiFilter.MHTSource = cms.InputTag("mhtPFchsNoMuon")

#jetMHTPFNoMuonIDDPhiFilter = jetMHTDPhiFilter.clone()
#jetMHTPFNoMuonIDDPhiFilter.JetSource = cms.InputTag("patJetsPFPt30NoMuonID")
#jetMHTPFNoMuonIDDPhiFilter.MHTSource = cms.InputTag("mhtPFchsNoMuonID")
#
#jetMHTPFNoMuonIDIsoDPhiFilter = jetMHTDPhiFilter.clone()
#jetMHTPFNoMuonIDIsoDPhiFilter.JetSource = cms.InputTag("patJetsPFPt30NoMuonIDIso")
#jetMHTPFNoMuonIDIsoDPhiFilter.MHTSource = cms.InputTag("mhtPFchsNoMuonIDIso")
