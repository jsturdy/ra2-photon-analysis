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
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_400_HT_inf_ext_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_999_1_zY9.root',
        '/store/user/lpcsusyhad/53X_ntuples/ZJetsToNuNu_400_HT_inf_ext_TuneZ2Star_8TeV_madgraph_Summer12/dhare/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_ext/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_985_1_8cs.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

process.source.skipEvents = cms.untracked.uint32(0)

##process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))
##
##process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )
##
##process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)  

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
from ZInvisibleBkgds.Photons.genstudytree_cfi import genzinvtree
process.st3ZBosons = genzinvtree.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZBosons"),
    debugString = cms.string("status 3 z's"),
)
process.st3ZNuBosons = genzinvtree.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZNuNuBosons"),
    debugString = cms.string("status 3 z's from di-nus"),
)

#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvMETProducers_cff')
process.load('ZInvisibleBkgds.Photons.ZinvVetos_cff')

from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import countPhotonsIDPFIso

process.countGenBosons = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZBosons"))
process.countGenNuNu   = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZNuNuBosons"))

process.analysisSeq = cms.Sequence(process.ra2PFchsJets
                                 * process.htPFchs
                                 * process.mhtPFchs
                                 * process.zinvBkgdGenZBosons
#                                 * process.zinvBkgdGenZNuNuBosons
                                 * process.zinvBJetsPF
                                 * process.countGenBosons
#                                 * process.countGenNuNu
                                 * process.zinvVetos
                                 * process.zinvBJetsPF
                                 * process.st3ZBosons
#                                 * process.st3ZNuBosons
)


#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zinvht400toinf_gen_tree.root')
)

#============================== configure paths ===============================
process.p1 = cms.Path(process.puWeight
                    * process.eventWeight
                    * process.analysisSeq
)
#process.pmu = cms.Path(process.muAnalysisSeq)
#process.pel = cms.Path(process.elAnalysisSeq)
#process.ptau = cms.Path(process.tauAnalysisSeq)
#process.pnu = cms.Path(process.nuAnalysisSeq)

##file = open('wtf_gentree.py','w')
##file.write(str(process.dumpPython()))
##file.close()
