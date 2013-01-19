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
            'file:/tmp/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph.root'
        )
    )
    process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
    process.source.skipEvents = cms.untracked.uint32(0)
    process.GlobalTag.globaltag = "START53_V7F::All"
    scaleF =  2.826*10*1000/1703863.

#========================= analysis module =====================================
process.zToMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string('patMuonsPFIDIso@+ patMuonsPFIDIso@-'),
                                 cut = cms.string('50.0 < mass < 120.0'),
                                 name = cms.string('zToMuMu'),
                                 roles = cms.vstring('muon1', 'muon2')
                                 )

from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.treemaker_cfi import dimuonTree
process.analysis = dimuonTree.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),
    metSource       = cms.InputTag("pfType1MetNoMuon","selected"),
    topTaggerSource = cms.string("myTopTagger"),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
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


from ZInvisibleBkgds.Photons.specialMETCleaner_cff import specialMuonCleanedMET
process.pfType1MetNoMuon = specialMuonCleanedMET.clone()
process.pfType1MetNoMuon.inputObjects = cms.InputTag("specialMuonCollection")

#process.pfMetType1NoMuon = cms.EDProducer('SpecialPATMuonCleanedMETProducer',
#    inputMET     = cms.InputTag('pfMetType1'),
#    inputObjects = cms.InputTag('specialMuonCollection'),
#                                          )

#================ analysis sequence =======================#

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.load('ZInvisibleBkgds.Photons.ZinvMuonJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')

process.patJetsPFchsPt30.src = cms.InputTag("patJetsAK5PFchs")
process.patJetsPFchsPt50Eta25.src = cms.InputTag("patJetsAK5PFchs")
process.patCSVMJetsPF.src = cms.InputTag('patJetsAK5PFchs')
process.patCSVTJetsPF.src = cms.InputTag('patJetsAK5PFchs')

process.patJetsPFNoMuon.src = cms.InputTag("patJetsAK5PFchs")
process.patJetsPFNoMuon.checkOverlaps.taus.src                = cms.InputTag('patTausPFchs')
process.patJetsPFNoMuon.checkOverlaps.electrons.src           = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoMuon.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoMuon.checkOverlaps.photons.src             = cms.InputTag('patPhotons')

process.patJetsPFNoMuonPt30.src = cms.InputTag("patJetsAK5PFchs")
process.patJetsPFNoMuonPt30.checkOverlaps.taus.src                = cms.InputTag('patTausPFchs')
process.patJetsPFNoMuonPt30.checkOverlaps.electrons.src           = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoMuonPt30.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoMuonPt30.checkOverlaps.photons.src             = cms.InputTag('patPhotons')

process.patJetsPFNoMuonPt50Eta25.src = cms.InputTag("patJetsAK5PFchs")
process.patJetsPFNoMuonPt50Eta25.checkOverlaps.taus.src                = cms.InputTag('patTausPFchs')
process.patJetsPFNoMuonPt50Eta25.checkOverlaps.electrons.src           = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoMuonPt50Eta25.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoMuonPt50Eta25.checkOverlaps.photons.src             = cms.InputTag('patPhotons')

process.load('ZInvisibleBkgds.Photons.zCandFilter_cff')
process.load('ZInvisibleBkgds.Photons.specialMuonCollection_cff')
process.specialMuonCollection.candidateLabel = cms.InputTag("zToMuMu")
process.patMuonsPFID.MuonSource    = cms.InputTag("patMuonsPFchs")
process.patMuonsPFIDIso.MuonSource = cms.InputTag("patMuonsPFchs")

from SandBox.Skims.RA2Objects_cff import countPFMuonsIDIso
process.countPFMuonsIDIsoForZ  = countPFMuonsIDIso.clone(minNumber = cms.uint32(2))

##top tagger
from UserCode.TopTagger.topTagger_cfi import *
process.load("UserCode.TopTagger.topTagger_cfi")
process.myTopTagger = topTagger.clone(
    metSrc = cms.InputTag("pfType1MetNoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)

process.myTopTagger4Loose = topTagger4Loose.clone(
    metSrc = cms.InputTag("pfType1MetNoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTagger5Loose = topTagger5Loose.clone(
    metSrc = cms.InputTag("pfType1MetNoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTagger6Loose = topTagger6Loose.clone(
    metSrc = cms.InputTag("pfType1MetNoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTagger4M = topTagger4M.clone(
    metSrc = cms.InputTag("pfType1MetNoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTagger5M = topTagger5M.clone(
    metSrc = cms.InputTag("pfType1MetNoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTagger6M = topTagger6M.clone(
    metSrc = cms.InputTag("pfType1MetNoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)

####
from SandBox.Skims.basicJetSelector_cfi import selectedPatJets
process.patJetsForIndirectTauVeto = selectedPatJets.clone()
process.patJetsForIndirectTauVeto.src = cms.InputTag("patJetsPFNoMuon")
process.patJetsForIndirectTauVeto.cut = cms.string('pt > 15 && abs(eta) < 2.4 && bDiscriminator("combinedSecondaryVertexBJetTags") <= 0.898')

#================ configure filters =======================#

from SandBox.HadSTopSkims.electronVetoFilter_cfi import electronVetoSTop
process.sTopPFElectronVeto = electronVetoSTop.clone(
    ElectronSource = cms.InputTag('patElectronsPFchs'),
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
    MaxChargedCandidates = cms.int32(2),
)
from SandBox.HadSTopSkims.indirectTauVeto_cfi import indirectTauVeto
process.sTopTauVeto = indirectTauVeto.clone(
    JetSource     = cms.InputTag("patJetsForIndirectTauVeto"),
    METSource     = cms.InputTag("pfType1MetNoMuon","selected"),
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
                                   * process.postStdCleaningCounter
                                   * process.htPFchs
                                   * process.mhtPFchs
                                   * process.zinvBJetsPF
                                   * process.patMuonsPFIDIso
                                   * process.countPFMuonsIDIsoForZ
                                   * process.zToMuMu
                                   * process.zCandFilter
                                   * process.specialMuonCollection
                                   * process.muonObjectsPF
                                   * process.patJetsForIndirectTauVeto
                                   * process.pfType1MetNoMuon
                                   * process.zinvBJetsPFNoMuon
                                   * process.sTopTrkIsolationMaker
                                   * process.sTopPFElectronVeto
                                   * process.sTopIsoTrkVeto
                                   * process.sTopTauVeto
                                   ####Compute the top-tagger variables
                                   * process.myTopTagger4Loose
                                   * process.myTopTagger5Loose
                                   * process.myTopTagger6Loose
                                   * process.myTopTagger4M
                                   * process.myTopTagger5M
                                   * process.myTopTagger6M

                                   ####Run the analysis with muon-cleaned collections and no lepton vetos
                                   * process.analysis4Loose
                                   * process.analysis5Loose
                                   * process.analysis6Loose
                                   * process.analysis4M
                                   * process.analysis5M
                                   * process.analysis6M
)


#======================= output module configuration ===========================
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('dimuon.mc.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1')
    ),
    outputCommands = cms.untracked.vstring('drop *',
                                           'keep *_*pfMetType*_*_*',
                                           'keep *_*pfType1Met*_*_*',
                                           ),
    dropMetaData = cms.untracked.string('DROPPED')
)
                                                                      

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('doublemuon_HT_SAMPLENAME_sTop_JOBID.root')
)

#============================== configure paths ===============================
process.p1 = cms.Path(process.eventWeight
                    * process.puWeight
                    * process.analysisSeq )
process.outpath = cms.EndPath(process.out)
