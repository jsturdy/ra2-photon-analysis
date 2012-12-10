import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMaker")

#===================== Message Logger =============================
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
            limit = cms.untracked.int32(10),
            reportEvery = cms.untracked.int32(25)
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 250
process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(True)
            )

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

###================= configure poolsource module ===================
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'dcache:/pnfs/cms/WAX/11/store/user/lpcsusyhad/53X_ntuples/DoubleMu_Run2012B-13Jul2012-v4_lpc1/vchetlur/DoubleMu/DoubleMu_Run2012B-13Jul2012-v4_NOCUTS_HLTPFHTInc_12Oct2012V3_lpc1/062a2f358645c26cb65ce944f0ba30c4/susypat_100_1_9LL.root',
        '/store/user/lpcsusyhad/53X_ntuples/DoubleMu_Run2012B-13Jul2012-v4_lpc1/vchetlur/DoubleMu/DoubleMu_Run2012B-13Jul2012-v4_NOCUTS_HLTPFHTInc_12Oct2012V3_lpc1/062a2f358645c26cb65ce944f0ba30c4/susypat_100_1_9LL.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.source.skipEvents = cms.untracked.uint32(0)
process.GlobalTag.globaltag = "FT_53_V6C_AN2::All"

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
from ZInvisibleBkgds.Photons.treemaker_cfi import dimuonTree
process.analysis = dimuonTree.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(True),
    ScaleFactor     = cms.double(1.0),
    DoPUReweight    = cms.bool(False),

    metSource       = cms.InputTag("pfType1MetNoMuon","pfcand"),

    runTopTagger           = cms.bool(True),
    looseTopTaggerSource   = cms.string("diMuonTopTagger5Loose"),
    nominalTopTaggerSource = cms.string("diMuonTopTagger5M"),

    storeExtraVetos    = cms.bool(True),
    muonVetoSource     = cms.InputTag("sTopPFMuonVetoDiMuon"),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoDiMuon"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVetoDiLeptons"),
)

#================ analysis sequence =======================#

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.load('ZInvisibleBkgds.Photons.ZinvMuonJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')

process.load('ZInvisibleBkgds.Photons.zCandFilter_cff')
process.load('ZInvisibleBkgds.Photons.specialMuonCollection_cff')
process.specialMuonCollection.candidateLabel = cms.InputTag("zToMuMu")
process.patMuonsPFID.MuonSource    = cms.InputTag("patMuonsPF")
process.patMuonsPFIDIso.MuonSource = cms.InputTag("patMuonsPF")

from SandBox.Skims.RA2Objects_cff import countPFMuonsIDIso
process.countPFMuonsIDIsoForZ  = countPFMuonsIDIso.clone(minNumber = cms.uint32(2))

##top tagger

process.load('ZInvisibleBkgds.Photons.ZinvMETProducers_cff')
process.load('ZInvisibleBkgds.Photons.ZinvVetos_cff')
process.load('ZInvisibleBkgds.Photons.ZinvTopTaggers_cff')

process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
                                   * process.htPFchs
                                   * process.mhtPFchs
                                   * process.zinvBJetsPF
                                   * process.patMuonsPFIDIso
                                   * process.countPFMuonsIDIsoForZ
                                   * process.zToMuMu
                                   * process.zCandFilter
                                   * process.specialMuonCollection
                                   * process.muonObjectsPF
                                   * process.zmumuMETCollections
                                   * process.zmumuVetos
                                   * process.zmumuTopTaggers
                                   * process.zinvBJetsPFNoMuon
                                   * process.analysis
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('doublemuDataTree.root')
)

#============================== configure paths ===============================
process.p1 = cms.Path(process.eventWeight
#                    * process.puWeight
                    * process.analysisSeq )
