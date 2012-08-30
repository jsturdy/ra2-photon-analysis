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
        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_464_1_f7L.root',
        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_463_1_KMz.root',
        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_466_1_f5b.root',
        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_465_1_ELg.root',
        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_468_1_xut.root',
        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_467_1_yXy.root',
        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_46_1_ezp.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================

from ZInvisibleBkgds.Photons.genstudy_cfi import *
process.directPhotons = genstudy.clone(
    genLabel       = cms.InputTag("zinvBkgdDirectPhotons"),
    debugString = cms.string("direct photons"),
#    genJetLabel    = cms.InputTag("ak5GenJets"),
    #genJetLabel    = cms.InputTag("ak5GenJetsNoNu"),
#    doPUReweight = cms.bool(False),
#    puWeight     = cms.InputTag("puWeight"),
#    bosonMinPt    = cms.double(50.),
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
process.st1ZBosons = genstudy.clone(
    genLabel       = cms.InputTag("zinvBkgdst1ZBosons"),
    debugString = cms.string("status 1 z's"),
)
process.st3ZBosons = process.st1ZBosons.clone(
    genLabel       = cms.InputTag("zinvBkgdst3ZBosons"),
    debugString = cms.string("status 3 z's"),
)
##different gen jet collection
process.directPhotonsNoNu = process.directPhotons.clone(
    genJetLabel    = cms.InputTag("ak5GenJetsNoNu"),
    debugString = cms.string("direct photons gen jets no nu"),
)
process.secondaryPhotonsNoNu = process.directPhotonsNoNu.clone(
    genLabel       = cms.InputTag("zinvBkgdSecondaryPhotons"),
    debugString = cms.string("secondary photons gen jets no nu"),
)
process.fragmentationPhotonsNoNu = process.directPhotonsNoNu.clone(
    genLabel       = cms.InputTag("zinvBkgdFragmentationPhotons"),
    debugString = cms.string("fragmentation photons gen jets no nu"),
)
process.mistagPhotonsNoNu = process.directPhotonsNoNu.clone(
    genLabel       = cms.InputTag("zinvBkgdMistagPhotons"),
    debugString = cms.string("mistag photons gen jets no nu"),
)
process.st1ZBosonsNoNu = genstudy.clone(
    genJetLabel    = cms.InputTag("ak5GenJetsNoNu"),
    genLabel       = cms.InputTag("zinvBkgdst1ZBosons"),
    debugString = cms.string("status 1 z's gen jets no nu"),
)
process.st3ZBosonsNoNu = process.st1ZBosonsNoNu.clone(
    genLabel       = cms.InputTag("zinvBkgdst3ZBosons"),
    debugString = cms.string("status 3 z's gen jets no nu"),
)
##Acceptance study
process.directPhotonsAcc = process.directPhotons.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("direct photons acceptance study"),
)
process.secondaryPhotonsAcc = process.secondaryPhotons.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("secondary photons acceptance study"),
)
process.fragmentationPhotonsAcc = process.fragmentationPhotons.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("fragmentation photons acceptance study"),
)
process.mistagPhotonsAcc = process.mistagPhotons.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("mistag photons acceptance study"),
)
process.st1ZBosonsAcc = process.st1ZBosons.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("status 1 z's acceptance study"),
)
process.st3ZBosonsAcc = process.st3ZBosons.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("status 3 z's acceptance study"),
)
##different gen jet collection
process.directPhotonsNoNuAcc = process.directPhotonsNoNu.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("direct photons acceptance study gen jets no nu"),
)
process.secondaryPhotonsNoNuAcc = process.secondaryPhotonsNoNu.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("secondary photons acceptance study gen jets no nu"),
)
process.fragmentationPhotonsNoNuAcc = process.fragmentationPhotonsNoNu.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("fragmentation photons acceptance study gen jets no nu"),
)
process.mistagPhotonsNoNuAcc = process.mistagPhotonsNoNu.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("mistag photons acceptance study gen jets no nu"),
)
process.st1ZBosonsNoNuAcc = process.st1ZBosonsNoNu.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("status 1 z's acceptance study gen jets no nu"),
)
process.st3ZBosonsNoNuAcc = process.st3ZBosonsNoNu.clone(
    studyAcceptance    = cms.bool(True),
    debugString = cms.string("status 3 z's acceptance study gen jets no nu"),
)

##RecoIso study
process.directPhotonsRecoIso = process.directPhotons.clone(
    studyRecoIso    = cms.bool(True),
    debugString = cms.string("direct photons reco/isolation study"),
)
process.secondaryPhotonsRecoIso = process.secondaryPhotons.clone(
    studyRecoIso    = cms.bool(True),
    debugString = cms.string("secondary photons reco/isolation study"),
)
process.fragmentationPhotonsRecoIso = process.fragmentationPhotons.clone(
    studyRecoIso    = cms.bool(True),
    debugString = cms.string("fragmentation photons reco/isolation study"),
)
process.mistagPhotonsRecoIso = process.mistagPhotons.clone(
    studyRecoIso    = cms.bool(True),
    debugString = cms.string("mistag photons reco/isolation study"),
)
##different gen jet collection
process.directPhotonsNoNuRecoIso = process.directPhotonsNoNu.clone(
    studyRecoIso    = cms.bool(True),
    debugString = cms.string("direct photons reco/isolation study gen jets no nu"),
)
process.secondaryPhotonsNoNuRecoIso = process.secondaryPhotonsNoNu.clone(
    studyRecoIso    = cms.bool(True),
    debugString = cms.string("secondary photons reco/isolation study gen jets no nu"),
)
process.fragmentationPhotonsNoNuRecoIso = process.fragmentationPhotonsNoNu.clone(
    studyRecoIso    = cms.bool(True),
    debugString = cms.string("fragmentation photons reco/isolation study gen jets no nu"),
)
process.mistagPhotonsNoNuRecoIso = process.mistagPhotonsNoNu.clone(
    studyRecoIso    = cms.bool(True),
    debugString = cms.string("mistag photons reco/isolation study gen jets no nu"),
)

#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *

process.load('ZInvisibleBkgds.Photons.adduserdata_cfi')
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

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

process.analysisSeq = cms.Sequence(  process.ra2Objects
                                   * process.kt6PFJetsForIsolation
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.zinvBkgdGenPhotons
                                   * process.zinvPhotons
                                   * process.zinvBkgdGenZBosons
                                   * process.ra2PFMuonVeto
                                   * process.ra2PFElectronVeto

                                   * process.directPhotons
                                   * process.secondaryPhotons
                                   * process.fragmentationPhotons
                                   #* process.mistagPhotons
                                   #* process.st1ZBosons
                                   #* process.st3ZBosons

                                   * process.directPhotonsNoNu
                                   * process.secondaryPhotonsNoNu
                                   * process.fragmentationPhotonsNoNu
                                   #* process.mistagPhotonsNoNu
                                   #* process.st1ZBosonsNoNu
                                   #* process.st3ZBosonsNoNu

                                   * process.directPhotonsAcc
                                   * process.secondaryPhotonsAcc
                                   * process.fragmentationPhotonsAcc
                                   #* process.mistagPhotonsAcc
                                   #* process.st1ZBosonsAcc
                                   #* process.st3ZBosonsAcc

                                   * process.directPhotonsNoNuAcc
                                   * process.secondaryPhotonsNoNuAcc
                                   * process.fragmentationPhotonsNoNuAcc
                                   #* process.mistagPhotonsNoNuAcc
                                   #* process.st1ZBosonsNoNuAcc
                                   #* process.st3ZBosonsNoNuAcc

                                   * process.directPhotonsRecoIso
                                   * process.secondaryPhotonsRecoIso
                                   * process.fragmentationPhotonsRecoIso
                                   #* process.mistagPhotonsRecoIso

                                   * process.directPhotonsNoNuRecoIso
                                   * process.secondaryPhotonsNoNuRecoIso
                                   * process.fragmentationPhotonsNoNuRecoIso
                                   #* process.mistagPhotonsNoNuRecoIso
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonGenStudyRecoIso.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')

process.load('SandBox.Utilities.puWeightProducer_cfi')
##process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/May10_Prompt167151_pudist.root"
#process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/Cert_160404-177515_JSON.pileup.root"

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.p1 = cms.Path(process.puWeight * process.analysisSeq )
