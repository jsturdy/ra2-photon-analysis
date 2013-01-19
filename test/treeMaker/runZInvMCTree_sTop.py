import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMaker")

#===================== Message Logger =============================
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
            limit = cms.untracked.int32(10),
            reportEvery = cms.untracked.int32(25)
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 25
process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(True)
            )

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#================= configure poolsource module ===================
USECONDOR = False
condorMode = USECONDOR
if condorMode:
    process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )
    process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)
    process.GlobalTag.globaltag = "GLOBALTAG"
    scaleF = CROSSSEC*TOTALLUMI*1000/PROCESSEDEVENTS

else:
    process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(
            'file:/tmp/sturdy/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph.root'
        )
    )
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(250) )
    process.source.skipEvents = cms.untracked.uint32(0)
    process.GlobalTag.globaltag = "START53_V7F::All"
    scaleF = 5.274*10*1000/5095710.

#========================= analysis module =====================================

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
    metSource       = cms.InputTag("pfMetType1"),
    topTaggerSource = cms.string("myTopTagger"),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVeto"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVeto"),
)
process.analysis4Loose = process.analysis.clone(
    topTaggerSource = cms.string("myTopTagger4Loose"),
)
process.analysis5Loose = process.analysis4Loose.clone(
    topTaggerSource = cms.string("myTopTagger5Loose"),
)
process.analysis6Loose = process.analysis4Loose.clone(
    topTaggerSource = cms.string("myTopTagger6Loose"),
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

#================ analysis sequence=======================#

##top tagger
from UserCode.TopTagger.topTagger_cfi import *
process.load("UserCode.TopTagger.topTagger_cfi")
process.myTopTagger = topTagger.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsAK5PFchs"),
)
process.myTopTagger4Loose = topTagger4Loose.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsAK5PFchs"),
)
process.myTopTagger5Loose = topTagger5Loose.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsAK5PFchs"),
)
process.myTopTagger6Loose = topTagger6Loose.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsAK5PFchs"),
)
process.myTopTagger4M = topTagger4M.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsAK5PFchs"),
)
process.myTopTagger5M = topTagger5M.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsAK5PFchs"),
)
process.myTopTagger6M = topTagger6M.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsAK5PFchs"),
)

##
process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.patJetsPFchsPt30.src = cms.InputTag("patJetsAK5PFchs")
process.patCSVMJetsPF.src = cms.InputTag('patJetsAK5PFchs')
process.patCSVTJetsPF.src = cms.InputTag('patJetsAK5PFchs')

####
from SandBox.Skims.basicJetSelector_cfi import selectedPatJets
process.patJetsForIndirectTauVeto = selectedPatJets.clone()
process.patJetsForIndirectTauVeto.src = cms.InputTag("patJetsAK5PFchs")
process.patJetsForIndirectTauVeto.cut = cms.string('pt > 15 && abs(eta) < 2.4 && bDiscriminator("combinedSecondaryVertexBJetTags") <= 0.898')

#================ configure filters =======================#

from SandBox.HadSTopSkims.muonVetoFilter_cfi import muonVetoSTop
process.sTopPFMuonVeto = muonVetoSTop.clone()
from SandBox.HadSTopSkims.electronVetoFilter_cfi import electronVetoSTop
process.sTopPFElectronVeto = electronVetoSTop.clone(
    ElectronSource = cms.InputTag('patElectronsPFchs')
    #ElectronSource = cms.InputTag('gsfElectrons')
)
from SandBox.HadSTopSkims.trackIsolationMaker_cfi import trackIsolationMaker
process.sTopTrkIsolationMaker = trackIsolationMaker.clone(
    pfCandidatesTag     = cms.InputTag("pfNoPileUpPFchs"),
    vertexInputTag      = cms.InputTag("offlinePrimaryVertices"),
    dR_ConeSize         = cms.double(0.3),
    dz_CutValue         = cms.double(0.05),
    minPt_PFCandidate   = cms.double(10.0)
)
from SandBox.HadSTopSkims.isolatedTrackVeto_cfi import isolatedTrackVeto
process.sTopIsoTrkVeto = isolatedTrackVeto.clone(
    CandidatesCharge     = cms.InputTag("sTopTrkIsolationMaker","pfcandschg"),
    CandidatesTrkIso     = cms.InputTag("sTopTrkIsolationMaker","pfcandstrkiso"),
    CandidatesPt         = cms.InputTag("sTopTrkIsolationMaker","pfcandspt"),
    CandidatesDZ         = cms.InputTag("sTopTrkIsolationMaker","pfcandsdzpv"),
    MaxChargedCandidates = cms.int32(0),
)
from SandBox.HadSTopSkims.indirectTauVeto_cfi import indirectTauVeto
process.sTopTauVeto = indirectTauVeto.clone()

process.load('SandBox.Skims.RA2Cleaning_cff')
process.load('SandBox.Skims.RA2CleaningFilterResults_cfg')
process.prefilterCounter        = cms.EDProducer("EventCountProducer")
process.postStdCleaningCounter  = cms.EDProducer("EventCountProducer")

process.RA2_NoScraping = process.booleanFilter.clone(ResultSource = cms.InputTag("noscraping"))
#process.RA2_beamHaloFilter.ResultSource        = cms.InputTag("passCSCTightHaloFilter")
## Don't use this filter for now, needs to understand performance on physics events
#process.RA2_eeNoiseFilter.ResultSource         = cms.InputTag("passEENoiseFilter")
#process.RA2_trackingFailureFilter.ResultSource = cms.InputTag("passTrackingFailureFilter")
#process.RA2_inconsistentMuons.ResultSource     = cms.InputTag("passInconsistentMuonFilter")
#process.RA2_greedyMuons.ResultSource           = cms.InputTag("passGreedyMuonFilter")
#process.HcalLaserEventFilter.ResultSource      = cms.InputTag("passHCALLaserEventFilter")
#process.EEBadScFilter.ResultSource             = cms.InputTag("passEESuperCrystalFilter")
#process.RA2_EcalTPFilter.ResultSource          = cms.InputTag("passECALDeadCellTPFilter")
#process.RA2_HBHENoiseFilter = process.booleanFilter.clone(ResultSource = cms.InputTag("passHBHENoiseFilter"))
#process.RA2_BadPFMuonFilter = process.booleanFilter.clone(ResultSource = cms.InputTag("passBadPFMuonFilter"))
#process.EcalLaserFilter.ResultSource  = cms.InputTag("ecalLaserCorrFilter")
process.RA2_HBHENoiseFilterRA2.ResultSource = cms.InputTag("HBHENoiseFilterRA2","HBHENoiseFilterResult","TreeMaker")
#
process.muonPFCandidateProducer.PATMuonCollection = cms.InputTag("patMuonsPFchs")
process.eeNoiseFilter.taggingMode         = True
process.trackingFailureFilter.taggingMode = True
process.beamHaloFilter.taggingMode        = True
process.ra2NoiseCleaning.remove(process.HBHENoiseFilter)
#process.HBHENoiseFilterRA2.taggingMode    = True
process.inconsistentMuons.taggingMode     = True
process.greedyMuons.taggingMode           = True
process.ra2EcalTPFilter.taggingMode       = True
process.ra2EcalBEFilter.taggingMode       = True
process.hcalLaserEventFilter.taggingMode  = True
process.eeBadScFilter.taggingMode         = True
process.ecalLaserCorrFilter.taggingMode   = True

#process.ra2PostCleaning.remove(process.greedyMuons)
#process.ra2PostCleaning.remove(process.inconsistentMuons)
#process.ra2PostCleaning.remove(process.beamHaloFilter)
#process.ra2PostCleaning.remove(process.eeNoiseFilter)
process.ra2PostCleaning.remove(process.trackingFailureFilter)
#process.ra2PostCleaning.remove(process.hcalLaserEventFilter)
#process.ra2PostCleaning.remove(process.ra2EcalTPFilter)
#process.ra2PostCleaning.remove(process.eeBadScFilter)
#process.ra2PostCleaning.remove(process.HBHENoiseFilter)
#process.ra2PostCleaning.remove(process.HBHENoiseFilterRA2)
process.ra2PostCleaning.remove(process.ra2EcalBEFilter)
#process.ra2PostCleaning.remove(process.ra2NoiseCleaning)

#process.RA2_EcalBEFilter.ResultSource = cms.InputTag("ra2EcalBEFilter")

#process.ra2StdCleaning.remove()
#process.ra2NoiseCleaning.remove(process.ecalLaserCorrFilter)
process.ra2PostCleaning.remove(process.ecalLaserCorrFilter)

#process.RA2_EcalBEFilter.ResultSource = cms.InputTag("ra2EcalBEFilter")
process.cleaningSeq = cms.Sequence(process.oneGoodVertex
                                 #* process.RA2_NoScraping
                                 * process.RA2_beamHaloFilter
                                 # Don't use this filter for now, needs to understand performance on physics events
                                 #* process.RA2_eeNoiseFilter
                                 #* process.RA2_trackingFailureFilter
                                 * process.RA2_inconsistentMuons
                                 * process.RA2_greedyMuons
                                 * process.HcalLaserEventFilter
                                 * process.EEBadScFilter
                                 * process.RA2_EcalTPFilter
                                 #* process.RA2_EcalBEFilter
                                 #* process.RA2_HBHENoiseFilter
                                 * process.RA2_HBHENoiseFilterRA2
                                 #* process.RA2_BadPFMuonFilter
#                                 * process.EcalLaserFilter
)

######
process.analysisSeq = cms.Sequence(process.ra2PFchsJets
                                 * process.prefilterCounter
                                 * process.ra2StdCleaning
                                 * process.ra2PostCleaning
                                 * process.cleaningSeq
                                 * process.sTopTrkIsolationMaker
                                 * process.sTopPFMuonVeto
                                 * process.sTopPFElectronVeto
                                 * process.postStdCleaningCounter
                                 * process.htPFchs
                                 * process.mhtPFchs
                                 * process.patJetsForIndirectTauVeto
                                 * process.zinvBkgdGenZBosons
                                 * process.zinvBJetsPF
                                 * process.sTopIsoTrkVeto
                                 * process.sTopTauVeto
                                 ###compute top-tagger variables
                                 * process.myTopTagger4Loose
                                 * process.myTopTagger5Loose
                                 * process.myTopTagger6Loose
                                 * process.myTopTagger4M
                                 * process.myTopTagger5M
                                 * process.myTopTagger6M
                                 * process.analysis4Loose
                                 * process.analysis5Loose
                                 * process.analysis6Loose
                                 * process.analysis4M
                                 * process.analysis5M
                                 * process.analysis6M
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zinvisible_HT_SAMPLENAME_sTop_JOBID.root')
)

#============================== configure paths ===============================
process.p1 = cms.Path(process.eventWeight
                    * process.puWeight
                    * process.analysisSeq )
#file = open('zinv400tree.py','w')
#file.write(str(process.dumpPython()))
#file.close()
