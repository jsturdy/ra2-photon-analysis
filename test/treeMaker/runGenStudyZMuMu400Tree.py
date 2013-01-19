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
        '/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_400ToInf_8TeV_madgraph_Summer12_v1_ext/seema/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_v1_ext_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_9_1_5Kn.root',
        '/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_400ToInf_8TeV_madgraph_Summer12_v1_ext/seema/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_v1_ext_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_99_1_Tvr.root',
        '/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_400ToInf_8TeV_madgraph_Summer12_v1_ext/seema/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_v1_ext_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_98_1_VPV.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================

scaleF = 2.862*10*1000/2727789.
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
    genSrc          = cms.InputTag("zinvBkgdst1ZBosons"),
    debug           = cms.bool(False),
    debugString     = cms.string("status 1 z's"),
    ScaleFactor     = cms.double(scaleF),
    studyAcceptance = cms.bool(False),
    studyRecoIso    = cms.bool(False),
    recoJetSrc      = cms.InputTag("patJetsPFchsPt30"),
    htJetSrc        = cms.InputTag("patJetsPFchsPt50Eta25"),
    bJetSrc         = cms.InputTag("patCSVTJetsPFPt30Eta24"),
    htNoBosonSource  = cms.InputTag("htPFchs"),
    mhtNoBosonSource = cms.InputTag("mhtPFchs"),
    storeExtraVetos    = cms.bool(True),
    muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoZInv"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVeto"),
)
process.st3ZMuBosons = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZMuMuBosons"),
    debugString = cms.string("status 3 z's from di-muons"),
    muonVetoSource     = cms.InputTag("sTopPFMuonVetoDiMuon"),
    tauVetoSource      = cms.InputTag("sTopTauVetoDiMuon"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVetoDiLeptons"),
)
process.allGens = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("genParticles"),
    debugString = cms.string("all gen's"),
)
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')

from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import countPhotonsIDPFIso
process.countGenMuMu   = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZMuMuBosons"),
    minNumber = cms.uint32(1))


process.load('ZInvisibleBkgds.Photons.ZinvMETProducers_cff')
process.load('ZInvisibleBkgds.Photons.ZinvVetos_cff')

process.preAnalysisSeq = cms.Sequence(  process.ra2PFchsJets
                                      * process.zinvBkgdGenZBosons
                                      * process.zinvBkgdGenZMuMuBosons
                                      * process.zinvBJetsPF
                                      * process.muonObjectsPF
                                      * process.zmumuMETCollections
                                      * process.zmumuVetos

)
process.muAnalysisSeq = cms.Sequence(process.countGenMuMu
                                   * process.st3ZMuBosons
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zmumuht400toinf_gen_tree.root')
)

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
#process.p1 = cms.Path(process.puWeight
#                    * process.eventWeight
#                    * process.preAnalysisSeq
#)
process.pmu = cms.Path(process.puWeight
                     * process.eventWeight
                     * process.preAnalysisSeq
                     * process.muAnalysisSeq)

##file = open('wtf_gentree.py','w')
##file.write(str(process.dumpPython()))
##file.close()
