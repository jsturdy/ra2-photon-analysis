import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMaker")

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

#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(
#        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_464_1_f7L.root',
#        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_463_1_KMz.root',
#        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_466_1_f5b.root',
#        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_465_1_ELg.root',
#        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_468_1_xut.root',
#        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_467_1_yXy.root',
#        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_46_1_ezp.root'
#    )
#)
#
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#
#process.source.skipEvents = cms.untracked.uint32(0)
#
#FILELIST = ['file:/eos/uscms/store/user/seema/SusyRA2Analysis2012/06Aug2012_ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12V3/susypat_34_1_4f8.root']
#MAXEVENTS = 1000
#SKIPEVENTS = 0
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )

process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)  

#========================= analysis module =====================================

scaleF = 5.274*10*1000/1006928.
process.analysis = cms.EDAnalyzer('RA2ZInvTreeMaker',
                                  Debug           = cms.bool(False),
                                  Data            = cms.bool(False),
                                  ScaleFactor     = cms.double(scaleF),
                                  genLabel        = cms.InputTag("zinvBkgdst3ZBosons"),
                                  PhotonSrc       = cms.InputTag("patPhotonsUserData"),
                                  VertexSrc       = cms.InputTag("goodVertices"),
                                  JetSrc          = cms.InputTag("patJetsPFPt30"),
                                  bJetSrc         = cms.InputTag("patCSVJetsPFPt30Eta24"),
                                  JetHTSource     = cms.InputTag("patJetsPFPt50Eta25"),
                                  DoPUReweight    = cms.bool(True),
                                  PUWeightSource  = cms.InputTag("puWeight")
)

#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *

process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                   process.ra2ObjectsPF
                                   * process.kt6PFJetsForIsolation
                                   * process.zinvBkgdGenZBosons
                                   * process.zinvBJetsPF
                                   * process.ra2MuonVeto
                                   * process.ra2ElectronVeto
                                   * process.analysis
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zinvisible_HT_400toInf_JOBID.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')

process.load('SandBox.Utilities.puWeightProducer_cfi')
##process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/May10_Prompt167151_pudist.root"
#process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/Cert_160404-177515_JSON.pileup.root"

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.p1 = cms.Path(process.puWeight * process.analysisSeq )
