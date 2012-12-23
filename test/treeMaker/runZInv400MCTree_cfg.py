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
FILELIST = ['file:/eos/uscms/store/user/seema/SusyRA2Analysis2012/06Aug2012_ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12V3/susypat_34_1_4f8.root']
MAXEVENTS = 1000
SKIPEVENTS = 0
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )

process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)  

#========================= analysis module =====================================

scaleF = 5.274*10*1000/1006928.
from ZInvisibleBkgds.Photons.treemaker_cfi import zvvTree
process.analysis = zvvTree.clone(
    Debug           = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
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

process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                     process.ra2PFchsJets
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
