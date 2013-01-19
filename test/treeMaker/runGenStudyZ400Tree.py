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
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_400_HT_inf_ext_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_9_1_gtL.root',
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_400_HT_inf_ext_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_99_1_MVz.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================

scaleF = 6.26*10*1000/5105710.
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
process.st3ZBosons = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZBosons"),
    debugString = cms.string("status 3 z's"),
)
process.st3ZMuBosons = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZMuMuBosons"),
    debugString = cms.string("status 3 z's from di-muons"),
    muonVetoSource     = cms.InputTag("sTopPFMuonVetoDiMuon"),
    tauVetoSource      = cms.InputTag("sTopTauVetoDiMuon"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVetoDiLeptons"),
)
process.st3ZElBosons = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZElElBosons"),
    debugString = cms.string("status 3 z's from di-electrons"),
    electronVetoSource = cms.InputTag("sTopPFElectronVetoDiElectron"),
    tauVetoSource      = cms.InputTag("sTopTauVetoDiElectron"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVetoDiLeptons"),
)
process.st3ZTauBosons = process.st1ZBosons.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZTauTauBosons"),
    debugString = cms.string("status 3 z's from di-taus"),
    muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoDiTau"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVetoDiLeptons"),
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

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
#from ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff import *
#from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import countPhotonsIDPFIso

process.countGenBosons = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZBosons"))
process.countGenNuNu   = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZNuNuBosons"),
    minNumber = cms.uint32(1))
process.countGenMuMu   = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZMuMuBosons"),
    minNumber = cms.uint32(1))
process.countGenElEl   = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZElElBosons"),
    minNumber = cms.uint32(1))
process.countGenTauTau = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZTauTauBosons"),
    minNumber = cms.uint32(1))

process.load('ZInvisibleBkgds.Photons.ZinvVetos_cff')

process.preAnalysisSeq = cms.Sequence(  process.ra2PFchsJets
                                      * process.zinvVetos
                                      * process.zinvBkgdGenZBosons
                                      * process.zinvBkgdGenZElElBosons
                                      * process.zinvBkgdGenZMuMuBosons
                                      * process.zinvBkgdGenZTauTauBosons
                                      * process.zinvBkgdGenZNuNuBosons
                                      * process.zinvBJetsPF
                                      #* process.countGenBosons

)
process.muAnalysisSeq = cms.Sequence(process.countGenMuMu
                                   * process.st3ZMuBosons
)
process.elAnalysisSeq = cms.Sequence(process.countGenElEl
                                   * process.st3ZElBosons
)
process.tauAnalysisSeq = cms.Sequence(process.countGenTauTau
                                    * process.st3ZTauBosons
)
process.nuAnalysisSeq = cms.Sequence(#process.countGenBosons
                                     process.countGenNuNu
                                   * process.st3ZNuBosons
)


#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zinvht400toinf_gen_tree.root')
)

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
#process.p1 = cms.Path(process.puWeight
#                    * process.eventWeight
#                    * process.preAnalysisSeq
#)
#process.pmu = cms.Path(process.puWeight
#                     * process.eventWeight
#                     * process.preAnalysisSeq
#                     * process.muAnalysisSeq)
#process.pel = cms.Path(process.puWeight
#                     * process.eventWeight
#                     * process.preAnalysisSeq
#                     * process.elAnalysisSeq)
#process.ptau = cms.Path(process.puWeight
#                      * process.eventWeight
#                      * process.preAnalysisSeq
#                      * process.tauAnalysisSeq)
process.pnu = cms.Path(process.puWeight
                     * process.eventWeight
                     * process.preAnalysisSeq
                     * process.nuAnalysisSeq
)

##file = open('wtf_gentree.py','w')
##file.write(str(process.dumpPython()))
##file.close()
