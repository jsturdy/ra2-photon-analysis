import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

#===================== Message Logger =============================
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
            limit = cms.untracked.int32(10),
            reportEvery = cms.untracked.int32(100)
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 250
process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(True)
            )

#================= configure poolsource module ===================

#process.load('SusyAnalysis.PhotonAnalysis.PhotonRun2011AMay10ReReco_160404to163869_cfi');
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_98_1_YIN.root',
        '/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_99_1_Gav.root',
        '/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_96_1_KJU.root',
        '/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_94_1_9ZP.root',
        '/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_95_1_wIR.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================

from ZInvisibleBkgds.Photons.genstudytree_cfi import *
process.directPhotons = genstudytree.clone(
    debug = cms.bool(True),
    genLabel       = cms.InputTag("zinvBkgdDirectPhotons"),
    debugString = cms.string("direct photons"),
)
process.secondaryPhotons = process.directPhotons.clone(
    genLabel       = cms.InputTag("zinvBkgdSecondaryPhotons"),
    debugString = cms.string("secondary photons"),
)
process.fragmentationPhotons = process.directPhotons.clone(
    genLabel       = cms.InputTag("zinvBkgdFragmentationPhotons"),
    debugString = cms.string("fragmentation photons"),
)
process.mistagPhotons = process.directPhotons.clone(
    genLabel       = cms.InputTag("zinvBkgdMistagPhotons"),
    debugString = cms.string("mistag photons"),
)
process.st1ZBosons = genstudytree.clone(
    genLabel       = cms.InputTag("zinvBkgdst1ZBosons"),
    debugString = cms.string("status 1 z's"),
    studyAcceptance = cms.bool(False),
    studyRecoIso    = cms.bool(False),
)
process.st3ZBosons = process.st1ZBosons.clone(
    genLabel       = cms.InputTag("zinvBkgdst3ZBosons"),
    debugString = cms.string("status 3 z's"),
)
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *

process.load('ZInvisibleBkgds.Photons.adduserdata_cfi')

from ZInvisibleBkgds.Photons.photonmap_cfi import *
process.rhoToPhotonMap = photonmap.clone()
from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import *
process.patPhotonsUser1 = addphotonuserdata1.clone()
process.patPhotonsUser1.photonLabel = cms.InputTag("patPhotonsAlt")
process.patPhotonsUser1.userData.userFloats = cms.PSet(
    src = cms.VInputTag(
        cms.InputTag("rhoToPhotonMap")
    )
)
process.patPhotonsUserData = addphotonuserdata2.clone()
process.patPhotonsUserData.photonLabel = cms.InputTag("patPhotonsUser1")

process.countGenBosons  = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdDirectPhotons"))

process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.zinvBkgdGenPhotons
                                   * process.zinvBkgdGenZBosons
                                   * process.photonObjectsPF
                                   * process.zinvBJetsPFNoPhotonSpecial
                                   * process.zinvBJetsPF
                                   * process.countGenBosons
                                   * process.ra2MuonVeto
                                   * process.ra2ElectronVeto
####
                                   * process.directPhotons
                                   #* process.secondaryPhotons
                                   #* process.fragmentationPhotons
                                   #* process.mistagPhotons
                                   #* process.st1ZBosons
                                   #* process.st3ZBosons

)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonGenStudyTree.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')

process.load('SandBox.Utilities.puWeightProducer_cfi')
##process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/May10_Prompt167151_pudist.root"
#process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/Cert_160404-177515_JSON.pileup.root"

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.p1 = cms.Path(process.puWeight * process.analysisSeq )
##file = open('wtf_gentree.py','w')
##file.write(str(process.dumpPython()))
##file.close()
