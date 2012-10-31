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
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/lpcsusyhad/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12/samantha/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_100_1_yAC.root',
        '/store/user/lpcsusyhad/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12/samantha/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_101_1_iUW.root',
        '/store/user/lpcsusyhad/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12/samantha/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_102_1_Q7x.root',
        '/store/user/lpcsusyhad/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12/samantha/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_107_1_flB.root',
        '/store/user/lpcsusyhad/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12/samantha/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_103_1_ErO.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================

scaleF = 5.274*10*1000/1006928.
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.treemaker_cfi import photonTree
from ZInvisibleBkgds.Photons.treemaker_cfi import zvvTree
process.analysis = zvvTree.clone(
    Debug           = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    topTaggerSource = cms.string("myTopTagger"),
)
process.analysis4M = process.analysis.clone(
    topTaggerSource = cms.string("myTopTagger4M"),
)
process.analysis5M = process.analysis4M.clone(
    topTaggerSource = cms.string("myTopTagger5M"),
)
process.analysis6M = process.analysis4M.clone(
    topTaggerSource = cms.string("myTopTagger6M"),
)
process.analysis4T = process.analysis.clone(
    bJetSrc         = cms.InputTag("patCSVTJetsPFPt30Eta24"),
    topTaggerSource = cms.string("myTopTagger4T"),
)
process.analysis5T = process.analysis4T.clone(
    topTaggerSource = cms.string("myTopTagger5T"),
)
process.analysis6T = process.analysis4T.clone(
    topTaggerSource = cms.string("myTopTagger6T"),
)
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *

##top tagger
from UserCode.TopTagger.topTagger_cfi import *
process.load("UserCode.TopTagger.topTagger_cfi")
process.myTopTagger = topTagger.clone()
process.myTopTagger.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTagger4M = topTagger4M.clone()
process.myTopTagger4M.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTagger5M = topTagger5M.clone()
process.myTopTagger5M.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTagger6M = topTagger6M.clone()
process.myTopTagger6M.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTagger4T = topTagger4T.clone()
process.myTopTagger4T.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTagger5T = topTagger5T.clone()
process.myTopTagger5T.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTagger6T = topTagger6T.clone()
process.myTopTagger6T.jetSrc = cms.InputTag("patJetsPFchsPt30")
      
process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                     process.ra2MuonVeto
                                   * process.ra2ElectronVeto
                                   * process.ra2PFchsJets
                                   * process.zinvBkgdGenZBosons
                                   * process.zinvBJetsPF
                                   * process.myTopTagger4M
                                   * process.myTopTagger5M
                                   * process.myTopTagger6M
                                   * process.myTopTagger4T
                                   * process.myTopTagger5T
                                   * process.myTopTagger6T
                                   * process.analysis4M
                                   * process.analysis5M
                                   * process.analysis6M
                                   * process.analysis4T
                                   * process.analysis5T
                                   * process.analysis6T
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zinvisible_HT_400toInf.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')

#process.load('SandBox.Utilities.puWeightProducer_cfi')
##process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/May10_Prompt167151_pudist.root"
#process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/Cert_160404-177515_JSON.pileup.root"

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.p1 = cms.Path(process.eventWeight
                    * process.puWeight
                    * process.analysisSeq )
#file = open('zinv400tree.py','w')
#file.write(str(process.dumpPython()))
#file.close()
