import FWCore.ParameterSet.Config as cms

process = cms.Process("TemplateMaker")

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

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
##process.GlobalTag.globaltag = "START52_V5::All"
##if runningOnMC == False:
#process.GlobalTag.globaltag = "GR_R_52_V9D::All"

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
scaleF = 107.5*10*1000/1611963.
process.analysis = cms.EDAnalyzer('RA2ZInvPhotonTemplateMaker',
                                  Debug           = cms.bool(False),
                                  Data            = cms.bool(False),
                                  ScaleFactor     = cms.double(scaleF),
                                  PhotonSrc       = cms.InputTag("patPhotonsUserData"),
                                  JetSrc          = cms.InputTag("patJetsAK5PFPt30"),
                                  bJetSrc         = cms.InputTag("patCSVJetsAK5PFPt30Eta24"),
                                  JetHTSource     = cms.InputTag("patJetsAK5PFPt50Eta25"),
#                                  RA2NJets        = cms.uint32(3),
#                                  RA2HT           = cms.double(350.0),
#                                  RA2MHT          = cms.double(200.0),
#                                  RA2ApplyDphiCuts= cms.bool(True),
                                  DoPUReweight    = cms.bool(True),
                                  PUWeightSource  = cms.InputTag("puWeight"),
)

#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.addphotonuserdata_cfi')
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

process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                   process.ra2Objects
                                   * process.kt6PFJetsForIsolation
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjects
                                   * process.zinvBJets
                                   * process.ra2MuonVeto
                                   * process.ra2ElectronVeto
                                   * process.analysis
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonMCTemplate.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')
process.load('SandBox.Utilities.puWeightProducer_cfi')
#============================== configure paths ===============================
process.p1 = cms.Path(process.puWeight * process.analysisSeq )
##process.p1 = cms.Path( process.runRangeFilter1 * process.analysisSeq )  #160431 - 161176
