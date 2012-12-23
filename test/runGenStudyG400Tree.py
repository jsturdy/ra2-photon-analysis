import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

#===================== Message Logger =============================
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
            limit = cms.untracked.int32(10),
            reportEvery = cms.untracked.int32(100)
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(True)
            )

#================= configure poolsource module ===================

#process.load('SusyAnalysis.PhotonAnalysis.PhotonRun2011AMay10ReReco_160404to163869_cfi');
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/lpcsusyhad/53X_ntuples/GJets_HT-400ToInf_8TeV_madgraph_v2_Summer12/dhare/GJets_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_1849_1_cum.root',
        '/store/user/lpcsusyhad/53X_ntuples/GJets_HT-400ToInf_8TeV_madgraph_v2_Summer12/dhare/GJets_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_184_1_VO8.root',
        '/store/user/lpcsusyhad/53X_ntuples/GJets_HT-400ToInf_8TeV_madgraph_v2_Summer12/dhare/GJets_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_1850_1_IK2.root',
        '/store/user/lpcsusyhad/53X_ntuples/GJets_HT-400ToInf_8TeV_madgraph_v2_Summer12/dhare/GJets_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_1851_1_r0Z.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source.skipEvents = cms.untracked.uint32(0)
#========================= analysis module =====================================

scaleF = 107.5*10*1000/9534744.
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.genstudytree_cfi import *
process.directPhotonsID = genstudytree.clone(
    debug = cms.bool(False),
    genSrc         = cms.InputTag("zinvBkgdDirectPhotons"),
    debugString = cms.string("direct photons"),
    ScaleFactor     = cms.double(scaleF),
    recoPhotonSrc = cms.InputTag("patPhotonsID"),
    recoJetSrc    = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30"),
    htJetSrc      = cms.InputTag("patJetsPFNoPhotonIDSpecialPt50Eta25"),
    bJetSrc       = cms.InputTag("patCSVTJetsPFNoPhotonIDSpecialPt30Eta24"),
    htNoBosonSource  = cms.InputTag("htPFchsNoPhotID"),
    mhtNoBosonSource = cms.InputTag("mhtPFchsNoPhotID"),
    storeExtraVetos    = cms.bool(True),
    muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoZInv"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVeto"),
)
process.directPhotonsIDPFIso = process.directPhotonsID.clone(
    recoPhotonSrc = cms.InputTag("patPhotonsIDPFIso"),
    recoJetSrc    = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30"),
    htJetSrc      = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25"),
    bJetSrc       = cms.InputTag("patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    htNoBosonSource  = cms.InputTag("htPFchsNoPhotIDPFIso"),
    mhtNoBosonSource = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
)
process.secondaryPhotonsID = process.directPhotonsID.clone(
    genSrc           = cms.InputTag("zinvBkgdSecondaryPhotons"),
    debugString      = cms.string("secondary photons"))
process.secondaryPhotonsIDPFIso = process.directPhotonsIDPFIso.clone(
    genSrc           = cms.InputTag("zinvBkgdSecondaryPhotons"),
    debugString      = cms.string("secondary photons"))
process.fragmentationPhotonsID = process.directPhotonsID.clone(
    genSrc           = cms.InputTag("zinvBkgdFragmentationPhotons"),
    debugString      = cms.string("fragmentation photons"))
process.fragmentationPhotonsIDPFIso = process.directPhotonsIDPFIso.clone(
    genSrc           = cms.InputTag("zinvBkgdFragmentationPhotons"),
    debugString      = cms.string("fragmentation photons"))
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
#from ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff import *
#from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
#from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
#from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *

process.load('ZInvisibleBkgds.Photons.adduserdata_cfi')

from ZInvisibleBkgds.Photons.photonmap_cfi import *
process.rhoToPhotonMap = photonmap.clone(photonSrc = cms.InputTag("patPhotonsRA2"))
from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import *
process.patPhotonsUser1 = addphotonuserdata1.clone()
process.patPhotonsUser1.photonLabel = cms.InputTag("patPhotonsRA2")
process.patPhotonsUser1.userData.userFloats = cms.PSet(
    src = cms.VInputTag(
        cms.InputTag("rhoToPhotonMap")
    )
)
process.patPhotonsUserData = addphotonuserdata2.clone()
process.patPhotonsUserData.photonLabel = cms.InputTag("patPhotonsUser1")

process.countDirectPhotons        = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdDirectPhotons"))
process.countSecondaryPhotons     = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdSecondaryPhotons"))
process.countFragmentationPhotons = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdFragmentationPhotons"))

process.load('ZInvisibleBkgds.Photons.ZinvVetos_cff')

process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.zinvBkgdGenPhotons
                                   * process.photonObjectsPF
                                   * process.photonVetos
                                   * process.zinvBJetsPFNoPhotonIDSpecial
                                   * process.zinvBJetsPFNoPhotonIDPFIsoSpecial
                                   * process.zinvBJetsPF
)

process.directAnalysisSeq = cms.Sequence(process.countDirectPhotons
                                       * process.directPhotonsID
                                       * process.directPhotonsIDPFIso
)

process.secondaryAnalysisSeq = cms.Sequence(process.countSecondaryPhotons
                                          * process.secondaryPhotonsID
                                          * process.secondaryPhotonsIDPFIso
)

process.fragmentationAnalysisSeq = cms.Sequence(process.countFragmentationPhotons
                                              * process.fragmentationPhotonsID
                                              * process.fragmentationPhotonsIDPFIso
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('gjetsht400toinf_gen_tree.root')
)

#============================== configure paths ===============================
process.p1 = cms.Path(process.puWeight
                    * process.eventWeight
                    * process.analysisSeq )
process.pdirect        = cms.Path(process.directAnalysisSeq)
process.psecondary     = cms.Path(process.secondaryAnalysisSeq)
process.pfragmentation = cms.Path(process.fragmentationAnalysisSeq)
