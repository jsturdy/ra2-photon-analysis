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
            'file:/tmp/sturdy/GJets_HT-400ToInf_8TeV-madgraph.root'
        )
    )
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
    process.source.skipEvents = cms.untracked.uint32(0)
    process.GlobalTag.globaltag = "START53_V7F::All"
    scaleF = 107.5*10*1000/9539562.

###========================= analysis module =====================================

from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.treemaker_cfi import photonTree
process.analysisID = photonTree.clone(
    DebugString     = cms.string("photonsID"),
    PhotonSrc       = cms.InputTag("patPhotonsID"),
    Debug           = cms.bool(False),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),
    metSource       = cms.InputTag("pfType1MetNoPhotonID","selected"),
    TriggerResults  = cms.InputTag("TriggerResults"),
    topTaggerSource = cms.string("myTopTaggerID"),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoID"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVeto"),
)
process.analysisID4Loose = process.analysisID.clone(
    DebugString     = cms.string("photonsID4Loose"),
    PhotonSrc       = cms.InputTag("patPhotonsID"),
    topTaggerSource = cms.string("myTopTaggerID4Loose"),
)
process.analysisID5Loose = process.analysisID4Loose.clone(
    DebugString     = cms.string("photonsID5Loose"),
    topTaggerSource = cms.string("myTopTaggerID5Loose"),
)
process.analysisID6Loose = process.analysisID4Loose.clone(
    DebugString     = cms.string("photonsID6Loose"),
    topTaggerSource = cms.string("myTopTaggerID6Loose"),
)

process.analysisIDPFIso = process.analysisID.clone(
    DebugString     = cms.string("photonsIDPFIso"),
    PhotonSrc       = cms.InputTag("patPhotonsIDPFIso"),
    JetSrc          = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30"),
    htJetSrc        = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25"),
    bJetSrc         = cms.InputTag("patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    htSource        = cms.InputTag("htPFchsNoPhotIDPFIso"),
    mhtSource       = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
    metSource       = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected"),
    topTaggerSource = cms.string("myTopTaggerIDPFIso5M"),
    tauVetoSource   = cms.InputTag("sTopTauVetoIDPFIso"),
)
process.analysisIDPFIso4Loose = process.analysisIDPFIso.clone(
    DebugString     = cms.string("photonsIDPFIso4Loose"),
    topTaggerSource = cms.string("myTopTaggerIDPFIso4Loose"),
)
process.analysisIDPFIso5Loose = process.analysisIDPFIso4Loose.clone(
    DebugString     = cms.string("photonsIDPFIso5Loose"),
    topTaggerSource = cms.string("myTopTaggerIDPFIso5Loose"),
)
process.analysisIDPFIso6Loose = process.analysisIDPFIso4Loose.clone(
    DebugString     = cms.string("photonsIDPFIso6Loose"),
    topTaggerSource = cms.string("myTopTaggerIDPFIso6Loose"),
)

process.analysisID4M = process.analysisID4Loose.clone(
    DebugString     = cms.string("photonsID4M"),
    topTaggerSource = cms.string("myTopTaggerID4M"),
)
process.analysisID5M = process.analysisID5Loose.clone(
    DebugString     = cms.string("photonsID6M"),
    topTaggerSource = cms.string("myTopTaggerID5M"),
)
process.analysisID6M = process.analysisID6Loose.clone(
    DebugString     = cms.string("photonsID6M"),
    topTaggerSource = cms.string("myTopTaggerID6M"),
)

process.analysisIDPFIso4M = process.analysisIDPFIso4Loose.clone(
    DebugString     = cms.string("photonsIDPFIso4M"),
    topTaggerSource = cms.string("myTopTaggerIDPFIso4M"),
)
process.analysisIDPFIso5M = process.analysisIDPFIso5Loose.clone(
    DebugString     = cms.string("photonsIDPFIso5M"),
    topTaggerSource = cms.string("myTopTaggerIDPFIso5M"),
)
process.analysisIDPFIso6M = process.analysisIDPFIso6Loose.clone(
    DebugString     = cms.string("photonsIDPFIso6M"),
    topTaggerSource = cms.string("myTopTaggerIDPFIso6M"),
)

#================ analysis sequence =======================#

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvPhotonJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')

process.patJetsPFchsPt30.src = cms.InputTag("patJetsAK5PFchs")
process.patCSVMJetsPF.src = cms.InputTag('patJetsAK5PFchs')
process.patCSVTJetsPF.src = cms.InputTag('patJetsAK5PFchs')

process.patJetsPFNoPhotonID.src = cms.InputTag("patJetsAK5PFchs")
process.patJetsPFNoPhotonID.checkOverlaps.taus.src                = cms.InputTag('patTausPFchs')
process.patJetsPFNoPhotonID.checkOverlaps.electrons.src           = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoPhotonID.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoPhotonID.checkOverlaps.muons.src               = cms.InputTag('patMuonsPFchs')

process.patJetsPFNoPhotonIDPFIso.src = cms.InputTag("patJetsAK5PFchs")
process.patJetsPFNoPhotonIDPFIso.checkOverlaps.taus.src                = cms.InputTag('patTausPFchs')
process.patJetsPFNoPhotonIDPFIso.checkOverlaps.electrons.src           = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoPhotonIDPFIso.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoPhotonIDPFIso.checkOverlaps.muons.src               = cms.InputTag('patMuonsPFchs')

process.patJetsPFNoPhotonIDSpecial.jetLabel      = cms.InputTag("patJetsAK5PFchs")
process.patJetsPFNoPhotonIDPFIsoSpecial.jetLabel = cms.InputTag("patJetsAK5PFchs")

from ZInvisibleBkgds.Photons.photonmap_cfi import *
process.rhoToPhotonMap = photonmap.clone(photonSrc = cms.InputTag("patPhotons"))
from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import *
process.patPhotonsUser1 = addphotonuserdata1.clone(debug = cms.bool(False))
process.patPhotonsUser1.photonLabel = cms.InputTag("patPhotons")
process.patPhotonsUser1.candidateLabel   = cms.InputTag("pfNoPileUpPFchs")
process.patPhotonsUser1.userData.userFloats = cms.PSet(
    src = cms.VInputTag(
        cms.InputTag("rhoToPhotonMap")
    )
)
process.patPhotonsUserData = addphotonuserdata2.clone(debug = cms.bool(False))
process.patPhotonsUserData.photonLabel = cms.InputTag("patPhotonsUser1")

process.countPhotonsID.src = cms.InputTag("patPhotonsID")
process.countPhotonsIDPFIso.src = cms.InputTag("patPhotonsIDPFIso")

from ZInvisibleBkgds.Photons.specialMETCleaner_cff import specialPhotonCleanedMET
process.pfType1MetNoPhotonIDPFIso = specialPhotonCleanedMET.clone()
process.pfType1MetNoPhotonIDPFIso.inputObjects = cms.InputTag("patPhotonsIDPFIso")
process.pfType1MetNoPhotonID = specialPhotonCleanedMET.clone()
process.pfType1MetNoPhotonID.inputObjects = cms.InputTag("patPhotonsID")

##top tagger
from UserCode.TopTagger.topTagger_cfi import *
process.load("UserCode.TopTagger.topTagger_cfi")

process.myTopTaggerID4Loose = topTagger4Loose.clone()
process.myTopTaggerID4Loose.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID4Loose.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
process.myTopTaggerID5Loose = topTagger5Loose.clone()
process.myTopTaggerID5Loose.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID5Loose.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
process.myTopTaggerID6Loose = topTagger6Loose.clone()
process.myTopTaggerID6Loose.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID6Loose.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
      
process.myTopTaggerIDPFIso4Loose = topTagger4Loose.clone()
process.myTopTaggerIDPFIso4Loose.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso4Loose.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")
process.myTopTaggerIDPFIso5Loose = topTagger5Loose.clone()
process.myTopTaggerIDPFIso5Loose.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso5Loose.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")
process.myTopTaggerIDPFIso6Loose = topTagger6Loose.clone()
process.myTopTaggerIDPFIso6Loose.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso6Loose.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")

process.myTopTaggerID4M = topTagger4M.clone()
process.myTopTaggerID4M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID4M.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
process.myTopTaggerID5M = topTagger5M.clone()
process.myTopTaggerID5M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID5M.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
process.myTopTaggerID6M = topTagger6M.clone()
process.myTopTaggerID6M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID6M.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
      
process.myTopTaggerIDPFIso4M = topTagger4M.clone()
process.myTopTaggerIDPFIso4M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso4M.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")
process.myTopTaggerIDPFIso5M = topTagger5M.clone()
process.myTopTaggerIDPFIso5M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso5M.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")
process.myTopTaggerIDPFIso6M = topTagger6M.clone()
process.myTopTaggerIDPFIso6M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso6M.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")

####
from SandBox.Skims.basicJetSelector_cfi import selectedPatJets
process.patJetsForIndirectTauVetoID = selectedPatJets.clone()
process.patJetsForIndirectTauVetoID.src = cms.InputTag("patJetsPFNoPhotonIDSpecial")
process.patJetsForIndirectTauVetoID.cut = cms.string('pt > 15 && abs(eta) < 2.4 && bDiscriminator("combinedSecondaryVertexBJetTags") <= 0.898')
process.patJetsForIndirectTauVetoIDPFIso = process.patJetsForIndirectTauVetoID.clone()
process.patJetsForIndirectTauVetoIDPFIso.src = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecial")

#================ configure filters =======================#

from SandBox.HadSTopSkims.muonVetoFilter_cfi import muonVetoSTop
process.sTopPFMuonVeto = muonVetoSTop.clone()
process.sTopPFMuonVeto = muonVetoSTop.clone(DoMuonVeto = cms.bool(True))
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
process.sTopTauVetoID = indirectTauVeto.clone(
    JetSource     = cms.InputTag("patJetsForIndirectTauVetoID"),
    METSource     = cms.InputTag("pfType1MetNoPhotonID","selected"),
)
process.sTopTauVetoIDPFIso = indirectTauVeto.clone(
    JetSource     = cms.InputTag("patJetsForIndirectTauVetoIDPFIso"),
    METSource     = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected"),
)

####
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
process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
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
#                                   * process.patJetsForIndirectTauVeto
                                   * process.zinvBJetsPF
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjectsPF
)

process.phoIDSeq = cms.Sequence(  process.countPhotonsID
                                * process.pfType1MetNoPhotonID
                                * process.patJetsForIndirectTauVetoID
                                * process.zinvBJetsPFNoPhotonIDSpecial
                                * process.sTopIsoTrkVeto
                                * process.sTopTauVetoID
                                * process.myTopTaggerID4Loose
                                * process.myTopTaggerID5Loose
                                * process.myTopTaggerID6Loose
                                * process.myTopTaggerID4M
                                * process.myTopTaggerID5M
                                * process.myTopTaggerID6M
                                * process.analysisID4Loose
                                * process.analysisID5Loose
                                * process.analysisID6Loose
                                * process.analysisID4M
                                * process.analysisID5M
                                * process.analysisID6M
)
process.phoIDPFIsoSeq = cms.Sequence( process.countPhotonsIDPFIso
                                    * process.pfType1MetNoPhotonIDPFIso
                                    * process.patJetsForIndirectTauVetoIDPFIso
                                    * process.zinvBJetsPFNoPhotonIDPFIsoSpecial
                                    * process.sTopIsoTrkVeto
                                    * process.sTopTauVetoIDPFIso
                                    * process.myTopTaggerIDPFIso5M
                                    * process.myTopTaggerIDPFIso4Loose
                                    * process.myTopTaggerIDPFIso5Loose
                                    * process.myTopTaggerIDPFIso6Loose
                                    * process.myTopTaggerIDPFIso4M
                                    * process.myTopTaggerIDPFIso5M
                                    * process.myTopTaggerIDPFIso6M
                                    * process.analysisIDPFIso4Loose
                                    * process.analysisIDPFIso5Loose
                                    * process.analysisIDPFIso6Loose
                                    * process.analysisIDPFIso4M
                                    * process.analysisIDPFIso5M
                                    * process.analysisIDPFIso6M
)

#======================= output module configuration ===========================
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('photon.mc.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('idiso')
    ),
    outputCommands = cms.untracked.vstring('drop *',
                                           'keep *_*pfMetType*_*_*',
                                           'keep *_*pfType1Met*_*_*',
                                           ),
    dropMetaData = cms.untracked.string('DROPPED')
)
                                                                      
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('gjets_HT_SAMPLENAME_sTop_JOBID.root')
)

#============================== configure paths ===============================
#process.p1 = cms.Path(process.eventWeight
#                    * process.puWeight
#                    * process.analysisSeq )
process.id    = cms.Path(process.eventWeight
                       * process.puWeight
                       * process.analysisSeq 
                       * process.phoIDSeq)
process.idiso = cms.Path(process.eventWeight
                       * process.puWeight
                       * process.analysisSeq 
                       * process.phoIDPFIsoSeq)
process.outpath = cms.EndPath(process.out)
