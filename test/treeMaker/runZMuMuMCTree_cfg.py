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

FILELIST = ['/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_473_1_xt4.root',
            '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_100_1_pAN.root',
            '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_101_1_fEs.root',
            '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_102_1_kr1.root',
            '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_103_1_3Ho.root',
            '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_104_1_LRu.root',
            '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_105_1_ceN.root',
            '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_107_1_qGd.root',
            '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_108_1_h7f.root',
            '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_110_1_AjL.root',
            ]
MAXEVENTS = -1
SKIPEVENTS = 0
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )

process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)  

#========================= analysis module =====================================
process.zToMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string('patMuonsPFIDIso@+ patMuonsPFIDIso@-'),
                                 cut = cms.string('50.0 < mass < 120.0'),
                                 name = cms.string('zToMuMu'),
                                 roles = cms.vstring('muon1', 'muon2')
                                 )

scaleF = 2950*10*1000/30460994.
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.treemaker_cfi import photonTree
from ZInvisibleBkgds.Photons.treemaker_cfi import dimuonTree
process.analysisNoMuon = dimuonTree.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),
)
process.analysisWithMuon = process.analysisNoMuon.clone()
process.analysisWithMuon.JetSrc    = cms.InputTag("patJetsPFchsPt30")
process.analysisWithMuon.htJetSrc  = cms.InputTag("patJetsPFchsPt50Eta25")
process.analysisWithMuon.bJetSrc   = cms.InputTag("patCSVJetsPFPt30Eta24")
process.analysisWithMuon.htSource  = cms.InputTag("htPFchs")
process.analysisWithMuon.mhtSource = cms.InputTag("mhtPFchs")

#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvMuonJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.zCandFilter_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *
process.load('ZInvisibleBkgds.Photons.specialMuonCollection_cff')
from ZInvisibleBkgds.Photons.specialMuonCollection_cff import specialMuonCollection
process.specialMuonCollection.candidateLabel = cms.InputTag("zToMuMu")
process.load('ZInvisibleBkgds.Photons.MuonHT_cff')
process.load('ZInvisibleBkgds.Photons.MuonMHT_cff')

from SandBox.Skims.RA2Objects_cff import countPFMuonsIDIso
process.countPFMuonsIDIsoForZ  = countPFMuonsIDIso.clone(minNumber = cms.uint32(2))
process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                     process.ra2PFchsJets
                                   * process.countPFMuonsIDIsoForZ
                                   * process.zToMuMu
                                   * process.zCandFilter
                                   * process.specialMuonCollection
                                   * process.muonCleanedPFJetsPF
                                   * process.zinvBJetsPFNoMuon
                                   * process.zinvBJetsPF
                                   * process.htPFchsNoMuon
                                   * process.mhtPFchsNoMuon
                                   * process.ra2ElectronVeto
                                   * process.analysisNoMuon
                                   * process.analysisWithMuon
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
    fileName = cms.string('dimuonMCTree_JOBID.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')

process.load('SandBox.Utilities.puWeightProducer_cfi')
##process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/May10_Prompt167151_pudist.root"
#process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/Cert_160404-177515_JSON.pileup.root"

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.p1 = cms.Path(process.eventWeight
                    * process.puWeight
                    * process.analysisSeq )
#process.outpath = cms.EndPath(process.out)
