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
        "file:/eos/uscms/store/user/seema/SusyRA2Analysis2012/06Aug2012_ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12V3/susypat_1_1_sVd.root"
       ,"file:/eos/uscms/store/user/seema/SusyRA2Analysis2012/06Aug2012_ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12V3/susypat_2_1_I3m.root"
       ,"file:/eos/uscms/store/user/seema/SusyRA2Analysis2012/06Aug2012_ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12V3/susypat_3_1_vSk.root"
       ,"file:/eos/uscms/store/user/seema/SusyRA2Analysis2012/06Aug2012_ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12V3/susypat_4_1_FFn.root"
       ,"file:/eos/uscms/store/user/seema/SusyRA2Analysis2012/06Aug2012_ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12V3/susypat_5_1_m1W.root"
       ,"file:/eos/uscms/store/user/seema/SusyRA2Analysis2012/06Aug2012_ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12V3/susypat_6_1_2on.root"
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source.skipEvents = cms.untracked.uint32(0)

#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))
#
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )
#
#process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)  

#========================= analysis module =====================================

from ZInvisibleBkgds.Photons.genstudytree_cfi import *
scaleF = 41.49*10*1000/5066608.
process.st1ZBosons = genstudytree.clone(
    genSrc          = cms.InputTag("zinvBkgdst1ZBosons"),
    debugString     = cms.string("status 1 z's"),
    ScaleFactor     = cms.double(scaleF),
    studyAcceptance = cms.bool(False),
    studyRecoIso    = cms.bool(False),
    recoJetSrc    = cms.InputTag("patJetsPFchsPt30"),
    htJetSrc      = cms.InputTag("patJetsPFchsPt50Eta25"),
    bJetSrc       = cms.InputTag("patCSVJetsPFPt30Eta24"),
    htSource      = cms.InputTag("htPFchs"),
    mhtSource     = cms.InputTag("mhtPFchs"),
    htNoBosonSource  = cms.InputTag("htPFchs"),
    mhtNoBosonSource = cms.InputTag("mhtPFchs"),
)
process.st3ZBosons = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZBosons"),
    debugString = cms.string("status 3 z's"),
)
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import countPhotonsIDPFIso

process.countGenBosons  = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdst3ZBosons"))

process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
                                   * process.zinvBkgdGenZBosons
                                   * process.zinvBJetsPF
                                   * process.countGenBosons
                                   * process.ra2MuonVeto
                                   * process.ra2ElectronVeto
####
                                   #* process.directPhotons
                                   #* process.secondaryPhotons
                                   #* process.fragmentationPhotons
                                   #* process.mistagPhotons
                                   #* process.st1ZBosons
                                   * process.st3ZBosons

)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zinvht200to400_gen_tree_JOBID.root')
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
