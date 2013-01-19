import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMaker")

#===================== Message Logger =============================
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
            limit = cms.untracked.int32(10),
            reportEvery = cms.untracked.int32(100)
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(True)
            )

#================= configure poolsource module ===================
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:///eos/uscms/store/user/bellan/PAT/Summer12_DR53X-PU_S10_START53_V7A-v1/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/v1/patuple_100_1_TSk.root',
        'file:/tmp/sturdy/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph.root'
        #'file:/tmp/sturdy/DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball.root'
        #'dcache:////pnfs/cms/WAX/resilient/bellan/PAT/Summer12_DR53X-PU_S10_START53_V7A-v2/DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball/v0/patuple_100_1_mrv.root',
        #'/store/user/lpcsusyhad/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12/samantha/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_100_1_yAC.root',
        #'/store/user/lpcsusyhad/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12/samantha/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_101_1_iUW.root',
        #'/store/user/lpcsusyhad/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12/samantha/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_102_1_Q7x.root',
        #'/store/user/lpcsusyhad/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12/samantha/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_107_1_flB.root',
        #'/store/user/lpcsusyhad/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph_Summer12/samantha/ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph/Summer12_PU_S7_START52_V9_RA2_NoCuts_09Aug2012V1/f42b6bbcc7ea08f6809cc21efa861dac/susypat_103_1_ErO.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================

scaleF = 5.274*10*1000/5095710.
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
process.analysis4MNoVeto = process.analysis4M.clone()
process.analysis5MNoVeto = process.analysis5M.clone()
process.analysis6MNoVeto = process.analysis6M.clone()

process.analysis4T = process.analysis.clone(
    bJetSrc         = cms.InputTag("patCSVTJetsPFPt30Eta24"),
    topTaggerSource = cms.string("myTopTagger4T"),
)
process.analysis5T = process.analysis4T.clone(
    topTaggerSource = cms.string("myTopTagger5T"),
)
process.analysis6T = process.analysis4T.clone(
    topTaggerSource = cms.string("myTopTagger6T"),
)
process.analysis4TNoVeto = process.analysis4T.clone()
process.analysis5TNoVeto = process.analysis5T.clone()
process.analysis6TNoVeto = process.analysis6T.clone()

#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *

##top tagger
from UserCode.TopTagger.topTagger_cfi import *
process.load("UserCode.TopTagger.topTagger_cfi")
process.myTopTagger = topTagger.clone(
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
process.myTopTagger4T = topTagger4T.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsAK5PFchs"),
)
process.myTopTagger5T = topTagger5T.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsAK5PFchs"),
)
process.myTopTagger6T = topTagger6T.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsAK5PFchs"),
)

####
from SandBox.Skims.basicJetSelector_cfi import selectedPatJets
process.patJetsForIndirectTauVeto = selectedPatJets.clone()
process.patJetsForIndirectTauVeto.src = cms.InputTag("patJetsAK5PFchs")
process.patJetsForIndirectTauVeto.cut = cms.string('pt > 15 && abs(eta) < 2.4 && bDiscriminator("combinedSecondaryVertexBJetTags") <= 0.898')

from SandBox.HadSTopSkims.muonVetoFilter_cfi import muonVetoSTop
process.sTopPFMuonVeto = muonVetoSTop.clone()
from SandBox.HadSTopSkims.electronVetoFilter_cfi import electronVetoSTop
process.sTopPFElectronVeto = electronVetoSTop.clone(
    ElectronSource = cms.InputTag('gsfElectrons'))
from SandBox.HadSTopSkims.indirectTauVeto_cfi import indirectTauVeto
process.sTopTauVeto = indirectTauVeto.clone()

process.patJetsPFchsPt30.src = cms.InputTag("patJetsAK5PFchs")
process.patCSVMJetsPF.src = cms.InputTag('patJetsAK5PFchs')
process.patCSVTJetsPF.src = cms.InputTag('patJetsAK5PFchs')
######
process.analysisSeq = cms.Sequence(process.ra2PFchsJets
                                 * process.htPFchs
                                 * process.mhtPFchs
#ra2Objects
                                 * process.patJetsForIndirectTauVeto
                                 * process.zinvBkgdGenZBosons
                                 * process.zinvBJetsPF
#ra2Objects
                                 * process.myTopTagger4M
                                 * process.myTopTagger5M
                                 * process.myTopTagger6M
                                 * process.myTopTagger4T
                                 * process.myTopTagger5T
                                 * process.myTopTagger6T
)
process.vetoSeq = cms.Sequence(process.sTopPFMuonVeto
                             * process.sTopPFElectronVeto
                             * process.sTopTauVeto
                             * process.analysis4M
                             * process.analysis5M
                             * process.analysis6M
                             * process.analysis4T
                             * process.analysis5T
                             * process.analysis6T
)

process.NoVetoSeq = cms.Sequence(process.analysis4MNoVeto
                               * process.analysis5MNoVeto
                               * process.analysis6MNoVeto
                               * process.analysis4TNoVeto
                               * process.analysis5TNoVeto
                               * process.analysis6TNoVeto
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zinvisible_HT_400toInf_sTop.root')
)

#============================== configure paths ===============================
process.p1 = cms.Path(process.eventWeight
                    * process.puWeight
                    * process.analysisSeq )
process.veto = cms.Path( process.vetoSeq )
process.noVeto = cms.Path(process.noVetoSeq )
#file = open('zinv400tree.py','w')
#file.write(str(process.dumpPython()))
#file.close()
