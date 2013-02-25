import FWCore.ParameterSet.Config as cms

from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
from ZInvisibleBkgds.Photons.PhotonHT_cff import *
from ZInvisibleBkgds.Photons.PhotonMHT_cff import *


photonObjectsPF = cms.Sequence(
      zinvPhotons 
    * photonCleanedPFJetsPF 
    * specialPhotonCleanedPFJetsPF 
    * htPFchsNoPhotID  
    * mhtPFchsNoPhotID 
    * htPFchsNoPhotIDPFIso  
    * mhtPFchsNoPhotIDPFIso
)
photonTemplateObjectsPF = cms.Sequence(
      photonTemplateCleanedPFJetsPF 
    * specialFitTemplatePhotonCleanedPFJetsPF 
    * htPFchsNoPhotFitTemplate  
    * mhtPFchsNoPhotFitTemplate
)
photonJetFakeObjectsPF = cms.Sequence(
      photonJetFakeCleanedPFJetsPF 
    * specialJetFakePhotonCleanedPFJetsPF 
    * htPFchsNoPhotJetFake  
    * mhtPFchsNoPhotJetFake
)

#from ZInvisibleBkgds.Photons.ZinvBkgdMuons_cff import *
from ZInvisibleBkgds.Photons.ZinvMuonJets_cff import *
from ZInvisibleBkgds.Photons.MuonHT_cff import *
from ZInvisibleBkgds.Photons.MuonMHT_cff import *

muonObjectsPF = cms.Sequence(
#    zinvMuons *
    muonCleanedPFJetsPF *
    htPFchsNoMuon  *
    mhtPFchsNoMuon
)

