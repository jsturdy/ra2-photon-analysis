import FWCore.ParameterSet.Config as cms

from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
from ZInvisibleBkgds.Photons.PhotonHT_cff import *
from ZInvisibleBkgds.Photons.PhotonMHT_cff import *

photonObjectsAK5PF = cms.Sequence(
    zinvPhotons *
    photonCleanedPFJetsAK5PF *
    htPFNoPhot  *
    mhtPFNoPhot
)

photonObjectsPF = cms.Sequence(
    zinvPhotons *
    photonCleanedPFJetsPF *
    specialPhotonCleanedPFJetsPF *
    htPFchsNoPhot  *
    mhtPFchsNoPhot
)

#from ZInvisibleBkgds.Photons.ZinvBkgdMuons_cff import *
from ZInvisibleBkgds.Photons.ZinvMuonJets_cff import *
from ZInvisibleBkgds.Photons.MuonHT_cff import *
from ZInvisibleBkgds.Photons.MuonMHT_cff import *

muonObjectsAK5PF = cms.Sequence(
#    zinvMuons *
    muonCleanedPFJetsAK5PF *
    htPFNoMuon  *
    mhtPFNoMuon
)

muonObjectsPF = cms.Sequence(
#    zinvMuons *
    muonCleanedPFJetsPF *
    htPFchsNoMuon  *
    mhtPFchsNoMuon
)

