import FWCore.ParameterSet.Config as cms

from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
from ZInvisibleBkgds.Photons.PhotonHT_cff import *
from ZInvisibleBkgds.Photons.PhotonMHT_cff import *

photonObjects = cms.Sequence(
    zinvPhotons *
    photonCleanedPFJets *
    htPFNoPhotIDIso  *
    mhtPFNoPhotIDIso
)

