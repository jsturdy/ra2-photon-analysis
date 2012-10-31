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
        '/store/user/lpcsusyhad/kasmi/2012AUG16/kasmi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1001_1_bVY.root',
        '/store/user/lpcsusyhad/kasmi/2012AUG16/kasmi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1000_1_vgx.root',
        '/store/user/lpcsusyhad/kasmi/2012AUG16/kasmi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1004_1_k6m.root',
        '/store/user/lpcsusyhad/kasmi/2012AUG16/kasmi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1003_1_pga.root',
        '/store/user/lpcsusyhad/kasmi/2012AUG16/kasmi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1002_1_mnK.root',
        '/store/user/lpcsusyhad/kasmi/2012AUG16/kasmi/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/Summer12-PU_S7_START52_V9-v1_NOCUTS_SkimsCode_09Aug2012V1/30d962f2384a73745773eb8ebda4b94d/SUSYPAT_1007_1_OcL.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================
process.zToMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string('patMuonsPFIDIso@+ patMuonsPFIDIso@-'),
                                 cut = cms.string('50.0 < mass < 120.0'),
                                 name = cms.string('zToMuMu'),
                                 roles = cms.vstring('muon1', 'muon2')
                                 )

from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(1.0),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.treemaker_cfi import photonTree
from ZInvisibleBkgds.Photons.treemaker_cfi import dimuonTree
process.analysisNoMuon = dimuonTree.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(True),
    ScaleFactor     = cms.double(1.0),
    DoPUReweight    = cms.bool(False),
    topTaggerSource = cms.string("myTopTaggerNoMuon"),
)
process.analysisNoMuon4M = process.analysisNoMuon.clone(
    topTaggerSource = cms.string("myTopTaggerNoMuon4M"),
)
process.analysisNoMuon5M = process.analysisNoMuon4M.clone(
    topTaggerSource = cms.string("myTopTaggerNoMuon5M"),
)
process.analysisNoMuon6M = process.analysisNoMuon4M.clone(
    topTaggerSource = cms.string("myTopTaggerNoMuon6M"),
)
process.analysisNoMuon4T = process.analysisNoMuon.clone(
    bJetSrc         = cms.InputTag("patCSVTJetsPFNoMuonPt30Eta24"),
    topTaggerSource = cms.string("myTopTaggerNoMuon4T"),
)
process.analysisNoMuon5T = process.analysisNoMuon4T.clone(
    topTaggerSource = cms.string("myTopTaggerNoMuon5T"),
)
process.analysisNoMuon6T = process.analysisNoMuon4T.clone(
    topTaggerSource = cms.string("myTopTaggerNoMuon6T"),
)

process.analysisWithMuon = process.analysisNoMuon.clone()
process.analysisWithMuon.JetSrc    = cms.InputTag("patJetsPFchsPt30")
process.analysisWithMuon.htJetSrc  = cms.InputTag("patJetsPFchsPt50Eta25")
process.analysisWithMuon.bJetSrc   = cms.InputTag("patCSVMJetsPFPt30Eta24")
process.analysisWithMuon.htSource  = cms.InputTag("htPFchs")
process.analysisWithMuon.mhtSource = cms.InputTag("mhtPFchs")
process.analysisWithMuon.topTaggerSource = cms.string("myTopTaggerWithMuon")
process.analysisWithMuon4M = process.analysisWithMuon.clone(
    topTaggerSource = cms.string("myTopTaggerWithMuon4M"),
)
process.analysisWithMuon5M = process.analysisWithMuon4M.clone(
    topTaggerSource = cms.string("myTopTaggerWithMuon5M"),
)
process.analysisWithMuon6M = process.analysisWithMuon4M.clone(
    topTaggerSource = cms.string("myTopTaggerWithMuon6M"),
)
process.analysisWithMuon4T = process.analysisWithMuon.clone(
    bJetSrc         = cms.InputTag("patCSVTJetsPFPt30Eta24"),
    topTaggerSource = cms.string("myTopTaggerWithMuon4T"),
)
process.analysisWithMuon5T = process.analysisWithMuon4T.clone(
    topTaggerSource = cms.string("myTopTaggerWithMuon5T"),
)
process.analysisWithMuon6T = process.analysisWithMuon4T.clone(
    topTaggerSource = cms.string("myTopTaggerWithMuon6T"),
)
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

from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import addphotonuserbasic
process.patPhotonsAlt = addphotonuserbasic.clone()

from SandBox.Skims.RA2Objects_cff import countPFMuonsIDIso
process.countPFMuonsIDIsoForZ  = countPFMuonsIDIso.clone(minNumber = cms.uint32(2))

##top tagger
from UserCode.TopTagger.topTagger_cfi import *
process.load("UserCode.TopTagger.topTagger_cfi")
process.myTopTaggerNoMuon = topTagger.clone()
process.myTopTaggerNoMuon.jetSrc = cms.InputTag("patJetsPFNoMuonPt30")
process.myTopTaggerNoMuon4M = topTagger4M.clone()
process.myTopTaggerNoMuon4M.jetSrc = cms.InputTag("patJetsPFNoMuonPt30")
process.myTopTaggerNoMuon5M = topTagger5M.clone()
process.myTopTaggerNoMuon5M.jetSrc = cms.InputTag("patJetsPFNoMuonPt30")
process.myTopTaggerNoMuon6M = topTagger6M.clone()
process.myTopTaggerNoMuon6M.jetSrc = cms.InputTag("patJetsPFNoMuonPt30")
process.myTopTaggerNoMuon4T = topTagger4T.clone()
process.myTopTaggerNoMuon4T.jetSrc = cms.InputTag("patJetsPFNoMuonPt30")
process.myTopTaggerNoMuon5T = topTagger5T.clone()
process.myTopTaggerNoMuon5T.jetSrc = cms.InputTag("patJetsPFNoMuonPt30")
process.myTopTaggerNoMuon6T = topTagger6T.clone()
process.myTopTaggerNoMuon6T.jetSrc = cms.InputTag("patJetsPFNoMuonPt30")
      
process.myTopTaggerWithMuon = topTagger.clone()
process.myTopTaggerWithMuon.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTaggerWithMuon4M = topTagger4M.clone()
process.myTopTaggerWithMuon4M.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTaggerWithMuon5M = topTagger5M.clone()
process.myTopTaggerWithMuon5M.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTaggerWithMuon6M = topTagger6M.clone()
process.myTopTaggerWithMuon6M.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTaggerWithMuon4T = topTagger4T.clone()
process.myTopTaggerWithMuon4T.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTaggerWithMuon5T = topTagger5T.clone()
process.myTopTaggerWithMuon5T.jetSrc = cms.InputTag("patJetsPFchsPt30")
process.myTopTaggerWithMuon6T = topTagger6T.clone()
process.myTopTaggerWithMuon6T.jetSrc = cms.InputTag("patJetsPFchsPt30")
      
process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                     process.ra2ElectronVeto
                                   * process.ra2PFchsJets
                                   * process.patPhotonsAlt
                                   * process.countPFMuonsIDIsoForZ
                                   * process.zToMuMu
                                   * process.zCandFilter
                                   * process.specialMuonCollection
                                   * process.muonCleanedPFJetsPF
                                   * process.zinvBJetsPFNoMuon
                                   * process.zinvBJetsPF
                                   * process.htPFchsNoMuon
                                   * process.mhtPFchsNoMuon
                                   * process.myTopTaggerNoMuon4M
                                   * process.analysisNoMuon4M
                                   * process.analysisWithMuon4M
                                   * process.myTopTaggerNoMuon5M
                                   * process.analysisNoMuon5M
                                   * process.analysisWithMuon5M
                                   * process.myTopTaggerNoMuon6M
                                   * process.analysisNoMuon6M
                                   * process.analysisWithMuon6M
                                   * process.myTopTaggerNoMuon4T
                                   * process.analysisNoMuon4T
                                   * process.analysisWithMuon4T
                                   * process.myTopTaggerNoMuon5T
                                   * process.analysisNoMuon5T
                                   * process.analysisWithMuon5T
                                   * process.myTopTaggerNoMuon6T
                                   * process.analysisNoMuon6T
                                   * process.analysisWithMuon6T
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('dimuonMCTree.root')
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
