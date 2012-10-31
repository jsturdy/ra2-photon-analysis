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
        '/store/user/lpcsusyhad/kasmi/kasmi/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1000_1_SQV.root',
        '/store/user/lpcsusyhad/kasmi/kasmi/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1001_1_InC.root',
        '/store/user/lpcsusyhad/kasmi/kasmi/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1002_1_RRw.root',
        '/store/user/lpcsusyhad/kasmi/kasmi/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1003_1_rat.root',
        '/store/user/lpcsusyhad/kasmi/kasmi/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1004_1_wfF.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source.skipEvents = cms.untracked.uint32(0)

##process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))
##
##process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )
##
##process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)  

#========================= analysis module =====================================

scaleF = 41.49*10*1000/5066608.
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.genstudytree_cfi import *
process.st1ZBosons = genstudytree.clone(
#    debug           = cms.bool(True),
    genSrc          = cms.InputTag("zinvBkgdst1ZBosons"),
    debugString     = cms.string("status 1 z's"),
    ScaleFactor     = cms.double(scaleF),
    studyAcceptance = cms.bool(False),
    studyRecoIso    = cms.bool(False),
    recoJetSrc      = cms.InputTag("patJetsPFchsPt30"),
    htJetSrc        = cms.InputTag("patJetsPFchsPt50Eta25"),
    bJetSrc         = cms.InputTag("patCSVTJetsPFPt30Eta24"),
    htNoBosonSource  = cms.InputTag("htPFchs"),
    mhtNoBosonSource = cms.InputTag("mhtPFchs"),
)
process.st3ZBosons = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZBosons"),
    debugString = cms.string("status 3 z's"),
)
process.st3ZMuBosonsNoVeto = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZMuMuBosons"),
    debugString = cms.string("status 3 z's from di-muons, no reco vetos"),
)
process.st3ZMuBosons = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZMuMuBosons"),
    debugString = cms.string("status 3 z's from di-muons"),
)
process.st3ZElBosonsNoVeto = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZElElBosons"),
    debugString = cms.string("status 3 z's from di-electrons, no reco vetos"),
)
process.st3ZElBosons = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZElElBosons"),
    debugString = cms.string("status 3 z's from di-electrons"),
)
process.st3ZTauBosonsNoVeto = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZTauTauBosons"),
    debugString = cms.string("status 3 z's from di-taus, no reco vetos"),
)
process.st3ZTauBosons = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZTauTauBosons"),
    debugString = cms.string("status 3 z's from di-taus"),
)
process.st3ZNuBosonsNoVeto = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZNuNuBosons"),
    debugString = cms.string("status 3 z's from di-nus, no reco vetos"),
)
process.st3ZNuBosons = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZNuNuBosons"),
    debugString = cms.string("status 3 z's from di-nus"),
)

process.allGens = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("genParticles"),
    debugString = cms.string("all gen's"),
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

process.countGenBosons = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZBosons"))
process.countGenNuNu   = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZNuNuBosons"),
    minNumber = cms.uint32(2))
process.countGenMuMu   = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZMuMuBosons"),
    minNumber = cms.uint32(2))
process.countGenElEl   = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZElElBosons"),
    minNumber = cms.uint32(2))
process.countGenTauTau = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZTauTauBosons"),
    minNumber = cms.uint32(2))

process.preAnalysisSeq = cms.Sequence(  process.ra2PFchsJets
                                      * process.zinvBkgdGenZBosons
                                      * process.zinvBJetsPF
                                      #* process.countGenBosons
                                      #* process.ra2MuonVeto
                                      #* process.ra2ElectronVeto
                                      #* process.st3ZBosons
)

process.muAnalysisSeq = cms.Sequence(process.countGenMuMu
                                   * process.st3ZMuBosonsNoVeto
                                   * process.ra2ElectronVeto
                                   * process.st3ZMuBosons
)
process.elAnalysisSeq = cms.Sequence(process.countGenElEl
                                   * process.st3ZElBosonsNoVeto
                                   * process.ra2MuonVeto
                                   * process.st3ZElBosons
)
process.tauAnalysisSeq = cms.Sequence(process.countGenTauTau
                                    * process.st3ZTauBosonsNoVeto
                                    * process.ra2MuonVeto
                                    * process.ra2ElectronVeto
                                    * process.st3ZTauBosons
)
process.nuAnalysisSeq = cms.Sequence(process.countGenNuNu
                                   * process.st3ZNuBosonsNoVeto
                                   * process.ra2MuonVeto
                                   * process.ra2ElectronVeto
                                   * process.st3ZNuBosons
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zinvht200to400_gen_tree.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')

#process.load('SandBox.Utilities.puWeightProducer_cfi')
##process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/May10_Prompt167151_pudist.root"
#process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/Cert_160404-177515_JSON.pileup.root"

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.p1 = cms.Path(process.puWeight
                    * process.eventWeight
                    * process.preAnalysisSeq
)
#process.pmu = cms.Path(process.muAnalysisSeq)
#process.pel = cms.Path(process.elAnalysisSeq)
#process.ptau = cms.Path(process.tauAnalysisSeq)
process.pnu = cms.Path(process.nuAnalysisSeq)

##file = open('wtf_gentree.py','w')
##file.write(str(process.dumpPython()))
##file.close()
