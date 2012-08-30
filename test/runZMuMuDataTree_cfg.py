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

##process.source = cms.Source("PoolSource",
##    fileNames = cms.untracked.vstring(
##        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_464_1_f7L.root',
##        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_463_1_KMz.root',
##        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_466_1_f5b.root',
##        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_465_1_ELg.root',
##        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_468_1_xut.root',
##        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_467_1_yXy.root',
##        '/store/user/sturdy07/RA2_525_Skims/GJets_HT400_cmslpc/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims_GJets_HT400_cmslpc/7d4ef27531e2177d5832a38a4c4fa602/susypat_mc_46_1_ezp.root'
##    )
##)
##
##process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
##
##process.source.skipEvents = cms.untracked.uint32(0)
##
##FILELIST = ['/store/user/lpcsusyhad/kasmi/2012AUG16/kasmi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1584_1_Mo2.root']
##MAXEVENTS = -1
##SKIPEVENTS = 0
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )

process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)  

#========================= analysis module =====================================
process.zToMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string('patMuonsIDIso@+ patMuonsIDIso@-'),
                                 cut = cms.string('50.0 < mass < 120.0'),
                                 name = cms.string('zToMuMu'),
                                 roles = cms.vstring('muon1', 'muon2')
                                 )

scaleF = 2950*10*1000/30460994.
process.analysis = cms.EDAnalyzer('RA2ZInvDiMuonTreeMaker',
                                  Debug           = cms.bool(False),
                                  Data            = cms.bool(False),
                                  ScaleFactor     = cms.double(1.0),
                                  MuonSrc         = cms.InputTag("specialMuonCollection"),
                                  VertexSrc       = cms.InputTag("goodVertices"),
                                  JetSrc          = cms.InputTag("patJetsPFNoMuonPt30"),
                                  htJetSrc        = cms.InputTag("patJetsPFNoMuonPt50Eta25"),
                                  bJetSrc         = cms.InputTag("patCSVJetsPFNoMuonPt30Eta24"),
                                  htSource        = cms.InputTag("htPFchsNoMuon"),
                                  mhtSource       = cms.InputTag("mhtPFchsNoMuon"),
                                  DoPUReweight    = cms.bool(False),
                                  PUWeightSource  = cms.InputTag("puWeight")
)

#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvMuonJets_cff')
process.load('ZInvisibleBkgds.Photons.zCandFilter_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *
process.load('ZInvisibleBkgds.Photons.specialMuonCollection_cff')
from ZInvisibleBkgds.Photons.specialMuonCollection_cff import specialMuonCollection
process.specialMuonCollection.candidateLabel = cms.InputTag("zToMuMu")
process.load('ZInvisibleBkgds.Photons.MuonHT_cff')
process.load('ZInvisibleBkgds.Photons.MuonMHT_cff')

process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                   process.ra2ObjectsPF
                                   * process.zToMuMu
                                   * process.zCandFilter
                                   * process.specialMuonCollection
                                   * process.muonCleanedPFJetsPF
                                   * process.zinvBJetsPFNoMuon
                                   * process.htPFchsNoMuon
                                   * process.mhtPFchsNoMuon
                                   * process.ra2ElectronVeto
                                   * process.analysis
)

#======================= output module configuration ===========================
##process.out = cms.OutputModule("PoolOutputModule",
##    fileName = cms.untracked.string('outputcolls.root'),
##    SelectEvents = cms.untracked.PSet(
##        SelectEvents = cms.vstring('p1')
##    ),
##    outputCommands = cms.untracked.vstring('keep *'),
##    dropMetaData = cms.untracked.string('DROPPED')
##)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('dimuonDataTree_JOBID.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')

process.load('SandBox.Utilities.puWeightProducer_cfi')
##process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/May10_Prompt167151_pudist.root"
#process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/Cert_160404-177515_JSON.pileup.root"

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.p1 = cms.Path(process.puWeight * process.analysisSeq )
#process.outpath = cms.EndPath(process.out)
