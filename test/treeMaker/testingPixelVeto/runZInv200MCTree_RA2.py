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

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#================= configure poolsource module ===================
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:///eos/uscms/store/user/bellan/PAT/Summer12_DR53X-PU_S10_START53_V7A-v1/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/v1/patuple_100_1_TSk.root',
        #'file:/tmp/sturdy/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph.root'
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_400_HT_inf_ext_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_9_1_gtL.root',
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_400_HT_inf_ext_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_99_1_MVz.root',
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_400_HT_inf_ext_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_999_1_zY9.root',
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_400_HT_inf_ext_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_985_1_8cs.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.source.skipEvents = cms.untracked.uint32(0)
process.GlobalTag.globaltag = "START53_V7F::All"
#========================= analysis module =====================================

scaleF = 49.28*10*1000/5055885.
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.treemaker_cfi import zvvTree
process.analysis = zvvTree.clone(
#    Debug           = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),

    metSource       = cms.InputTag("patMETsPF"),

    runTopTagger           = cms.bool(True),
    looseTopTaggerSource   = cms.string("zInvTopTagger5Loose"),
    nominalTopTaggerSource = cms.string("zInvTopTagger5M"),

    storeExtraVetos    = cms.bool(True),
    muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoZInv"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVeto"),
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

process.load('ZInvisibleBkgds.Photons.ZinvMETProducers_cff')
process.load('ZInvisibleBkgds.Photons.ZinvVetos_cff')
process.load('ZInvisibleBkgds.Photons.ZinvTopTaggers_cff')

process.load('SandBox.Skims.RA2CleaningFilterResults_cfg')
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
from SandBox.Skims.htFilter_cfi  import *
process.zinvHTFilter      = htFilter.clone(HTSource = cms.InputTag("htPFchs"),MinHT = cms.double(250))
from SandBox.Skims.mhtFilter_cfi import *
process.zinvMHTFilter      = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchs"),MinMHT = cms.double(100))

####
process.analysisSeq = cms.Sequence(process.ra2PFchsJets
                                 * process.htPFchs
                                 * process.mhtPFchs
                                 * process.ecalLaserCorrFilter
                                 * process.cleaningOnFilterResults
                                 * process.zinvHTFilter
                                 #* process.zinvMHTFilter
                                 * process.zinvVetos
                                 * process.zinvTopTaggers
                                 * process.zinvBkgdGenZBosons
                                 * process.zinvBJetsPF
                                 * process.analysis
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zinvisibleMC200Tree.root')
)

#============================== configure paths ===============================
process.p1 = cms.Path(process.eventWeight
                    * process.puWeight
                    * process.analysisSeq )
#file = open('zinv400tree.py','w')
#file.write(str(process.dumpPython()))
#file.close()
