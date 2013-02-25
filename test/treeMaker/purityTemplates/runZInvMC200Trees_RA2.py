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
        #'file:/uscms_data/d2/sturdy07/SUSY/RA2/newRelease/CMSSW_5_3_5/src/SandBox/Skims/test/susypat_mc_zinv.root'
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_1000_1_LpL.root',
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_1001_1_aWU.root',
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_100_1_lty.root',
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_104_1_yxr.root',
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_105_1_v3u.root',
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_101_1_itH.root',
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_108_1_XlH.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(25000) )
process.source.skipEvents = cms.untracked.uint32(0)
process.GlobalTag.globaltag = "START53_V7G::All"
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
    #Debug           = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),

    metSource       = cms.InputTag("newMETwPhiCorr"),

    runTopTagger           = cms.bool(True),
    looseTopTaggerSource   = cms.string("zInvTopTagger5Loose"),
    nominalTopTaggerSource = cms.string("zInvTopTagger5M"),

    storeExtraVetos    = cms.bool(True),
    muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoZInv"),
    isoTrkVetoSource   = cms.InputTag("sTopTrkIsolationMaker","trkIsoVeto"),
)
process.analysisGEN = process.analysis.clone(
    runTopTagger           = cms.bool(False),
)
#================ configure filters and analysis sequence=======================

process.load("SandBox.Skims.RA2Leptons_cff")
process.load("SandBox.Skims.jesChange_cfi")
process.newJetsMET.JECLevel = cms.string('ak5PFchsL1FastL2L3')
process.patMETPhiCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc
                        
process.load('SandBox.Skims.RA2Objects_cff')
process.patJetsPFchsPt30.src      = cms.InputTag('newJetsMET')
process.load('SandBox.Skims.RA2Selection_cff')
from SandBox.Skims.RA2Objects_cff import countJetsPFchsPt50Eta25
process.countJetsPFchsPt50Eta25.src = cms.InputTag('patJetsPFchsPt50Eta25')
process.countJetsPFchsPt50Eta25.minNumber = cms.uint32(2)

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')

process.load('ZInvisibleBkgds.Photons.ZinvMETProducers_cff')
process.load('ZInvisibleBkgds.Photons.ZinvVetos_cff')
process.load('ZInvisibleBkgds.Photons.ZinvTopTaggers_cff')

process.load('SandBox.Skims.RA2CleaningFilterResults_cfg')
process.load('SandBox.Skims.RA2CaloVsPFMHTFilterSequence_cff')
process.RA2CaloVsPFMHTFilter.TaggingMode = cms.bool(False)
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')

process.countZsGEN = process.countPhotonsID.clone()
process.countZsGEN.src = cms.InputTag("zinvBkgdst3ZBosons")
process.countZsGEN.minNumber = cms.uint32(1)

process.countMaxZsGEN = process.countZsGEN.clone()
process.countMaxZsGEN.maxNumber = cms.uint32(1)

from SandBox.Skims.jetMHTDPhiFilter_cfi  import *
process.zinvDPhiFilter   = jetMHTDPhiFilter.clone(MHTSource = cms.InputTag("mhtPFchs"),
                                                  JetSource = cms.InputTag("patJetsPFchsPt30"))
from SandBox.Skims.htFilter_cfi  import *
process.zinvHTPreFilter   = htFilter.clone(HTSource = cms.InputTag("htPFchs"),MinHT = cms.double(300.0))
process.zinvHTFilter      = htFilter.clone(HTSource = cms.InputTag("htPFchs"),MinHT = cms.double(500.0))
from SandBox.Skims.mhtFilter_cfi import *
process.zinvMHTPreFilter   = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchs"),MinMHT = cms.double(100))
process.zinvMHTFilter      = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchs"),MinMHT = cms.double(200))

####
process.setupSeq = cms.Sequence(
    process.newra2PFchsJets
    * process.ra2Electrons
    * process.ra2PFchsJets
    * process.htPFchs
    * process.mhtPFchs
    * process.zinvVetos
    * process.zinvTopTaggers
    * process.zinvBkgdGenZBosons
    * process.zinvBJetsPF
)
process.cleaningSeq = cms.Sequence(
    process.ecalLaserCorrFilter
    * process.cleaningOnFilterResults
    #* process.RA2CaloVsPFMHTFilterSequence
    #                                   * process.countRA2ElectronsIDIso
    #                                   * process.countRA2MuonsPFIDIso
)
process.analysisSeq = cms.Sequence(
    process.countJetsPFchsPt50Eta25
    * process.zinvHTPreFilter
    * process.zinvMHTPreFilter
    * process.analysis
    #* process.zinvDPhiFilter
    #* process.zinvHTFilter
    #* process.zinvMHTFilter
)
process.genAnalysisSeq = cms.Sequence(
    process.countZsGEN
    * process.analysisGEN
    #* process.countJetsPFchsPt50Eta25
    #* process.zinvHTPreFilter
    #* process.zinvMHTPreFilter
    #* process.zinvDPhiFilter
    #* process.zinvHTFilter
    #* process.zinvMHTFilter
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zinvisibleMC200Trees.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')

#process.load('SandBox.Utilities.puWeightProducer_cfi')
##process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/May10_Prompt167151_pudist.root"
#process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/Cert_160404-177515_JSON.pileup.root"

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.recoAnalysis = cms.Path(process.eventWeight
                             * process.puWeight
                             * process.cleaningSeq
                             * process.setupSeq
                             * process.analysisSeq
)
process.genAnalysis = cms.Path(process.eventWeight
                             * process.puWeight
                             * process.cleaningSeq
                             * process.setupSeq
                             * process.genAnalysisSeq
                               )
#process.withCaloPFMHT = cms.Path(process.eventWeight
#                               * process.puWeight
#                               * process.cleaningSeq
#                               * process.RA2CaloVsPFMHTFilterSequence
#                               * process.setupSeq
#                               * process.analysisSeq
#                               )
#file = open('zinv400tree.py','w')
#file.write(str(process.dumpPython()))
#file.close()
