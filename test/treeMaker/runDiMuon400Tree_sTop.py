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

#================= configure poolsource module ===================
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/tmp/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph.root'
        #'/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_101_1_fEs.root',
        #'/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_102_1_kr1.root',
        #'/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_103_1_3Ho.root',
        #'/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_104_1_LRu.root',
        #'/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_105_1_ceN.root',
        #'/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_107_1_qGd.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================
process.zToMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
                                 decay = cms.string('patMuonsPFIDIso@+ patMuonsPFIDIso@-'),
                                 cut = cms.string('50.0 < mass < 120.0'),
                                 name = cms.string('zToMuMu'),
                                 roles = cms.vstring('muon1', 'muon2')
                                 )

scaleF =  2.826*10*1000/1703863.
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.treemaker_cfi import photonTree
from ZInvisibleBkgds.Photons.treemaker_cfi import dimuonTree
process.analysisNoMuon = dimuonTree.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),
    metSource       = cms.InputTag("pfMetType1NoMuon","selected"),
    topTaggerSource = cms.string("myTopTaggerNoMuon"),
)
process.analysisNoMuon4M = process.analysisNoMuon.clone(
    topTaggerSource = cms.string("myTopTaggerNoMuon4M"),
)
process.analysisNoMuon5M = process.analysisNoMuon4M.clone(
    topTaggerSource = cms.string("myTopTaggerNoMuon5M"),
)
process.analysisNoMuon6M = process.analysisNoMuon4M.clone(
    topTaggerSource = cms.string("myTopTaggerNoMuon6M"),
)
process.analysisNoMuon4T = process.analysisNoMuon.clone(
    bJetSrc         = cms.InputTag("patCSVTJetsPFNoMuonPt30Eta24"),
    topTaggerSource = cms.string("myTopTaggerNoMuon4T"),
)
process.analysisNoMuon5T = process.analysisNoMuon4T.clone(
    topTaggerSource = cms.string("myTopTaggerNoMuon5T"),
)
process.analysisNoMuon6T = process.analysisNoMuon4T.clone(
    topTaggerSource = cms.string("myTopTaggerNoMuon6T"),
)

process.pfMetType1NoMuon = cms.EDProducer('SpecialPATMuonCleanedMETProducer',
    inputMET     = cms.InputTag('pfMetType1'),
    inputObjects = cms.InputTag('specialMuonCollection'),
                                          )
process.analysisWithMuon = process.analysisNoMuon.clone()
process.analysisWithMuon.JetSrc    = cms.InputTag("patJetsPFchsPt30")
process.analysisWithMuon.htJetSrc  = cms.InputTag("patJetsPFchsPt50Eta25")
process.analysisWithMuon.bJetSrc   = cms.InputTag("patCSVMJetsPFPt30Eta24")
process.analysisWithMuon.htSource  = cms.InputTag("htPFchs")
process.analysisWithMuon.mhtSource = cms.InputTag("mhtPFchs")
process.analysisWithMuon.metSource = cms.InputTag("pfMetType1")
process.analysisWithMuon.topTaggerSource = cms.string("myTopTaggerWithMuon")
process.analysisWithMuon4M = process.analysisWithMuon.clone(
    topTaggerSource = cms.string("myTopTaggerWithMuon4M"),
)
process.analysisWithMuon5M = process.analysisWithMuon4M.clone(
    topTaggerSource = cms.string("myTopTaggerWithMuon5M"),
)
process.analysisWithMuon6M = process.analysisWithMuon4M.clone(
    topTaggerSource = cms.string("myTopTaggerWithMuon6M"),
)
process.analysisWithMuon4T = process.analysisWithMuon.clone(
    bJetSrc         = cms.InputTag("patCSVTJetsPFPt30Eta24"),
    topTaggerSource = cms.string("myTopTaggerWithMuon4T"),
)
process.analysisWithMuon5T = process.analysisWithMuon4T.clone(
    topTaggerSource = cms.string("myTopTaggerWithMuon5T"),
)
process.analysisWithMuon6T = process.analysisWithMuon4T.clone(
    topTaggerSource = cms.string("myTopTaggerWithMuon6T"),
)
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvMuonJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.zCandFilter_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *
process.load('ZInvisibleBkgds.Photons.specialMuonCollection_cff')
from ZInvisibleBkgds.Photons.specialMuonCollection_cff import specialMuonCollection
process.specialMuonCollection.candidateLabel = cms.InputTag("zToMuMu")
process.load('ZInvisibleBkgds.Photons.MuonHT_cff')
process.load('ZInvisibleBkgds.Photons.MuonMHT_cff')

from SandBox.Skims.RA2Objects_cff import countPFMuonsIDIso
process.countPFMuonsIDIsoForZ  = countPFMuonsIDIso.clone(minNumber = cms.uint32(2))

##top tagger
from UserCode.TopTagger.topTagger_cfi import *
process.load("UserCode.TopTagger.topTagger_cfi")
process.myTopTaggerNoMuon = topTagger.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)

process.myTopTaggerNoMuon4M = topTagger4M.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTaggerNoMuon5M = topTagger5M.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTaggerNoMuon6M = topTagger6M.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTaggerNoMuon4T = topTagger4T.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTaggerNoMuon5T = topTagger5T.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTaggerNoMuon6T = topTagger6T.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)

process.myTopTaggerWithMuon4M = topTagger4M.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsPFchsPt30"),
)
process.myTopTaggerWithMuon5M = topTagger5M.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsPFchsPt30"),
)
process.myTopTaggerWithMuon6M = topTagger6M.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsPFchsPt30"),
)
process.myTopTaggerWithMuon4T = topTagger4T.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsPFchsPt30"),
)
process.myTopTaggerWithMuon5T = topTagger5T.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsPFchsPt30"),
)
process.myTopTaggerWithMuon6T = topTagger6T.clone(
    metSrc = cms.InputTag("pfMetType1"),
    jetSrc = cms.InputTag("patJetsPFchsPt30"),
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
    ElectronSource = cms.InputTag('patElectronsPFchs'))
from SandBox.HadSTopSkims.indirectTauVeto_cfi import indirectTauVeto
process.sTopTauVeto = indirectTauVeto.clone()

process.patJetsPFchsPt30.src = cms.InputTag("patJetsAK5PFchs")
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

process.patMuonsPFID.MuonSource    = cms.InputTag("patMuonsPFchs")
process.patMuonsPFIDIso.MuonSource = cms.InputTag("patMuonsPFchs")
      
process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
                                   * process.htPFchs
                                   * process.mhtPFchs
                                   * process.patJetsForIndirectTauVeto
                                   * process.zinvBJetsPF
                                   * process.patMuonsPFIDIso
                                   * process.countPFMuonsIDIsoForZ
                                   * process.zToMuMu
                                   * process.zCandFilter
                                   * process.specialMuonCollection
                                   * process.pfMetType1NoMuon
                                   * process.muonCleanedPFJetsPF
                                   * process.zinvBJetsPFNoMuon
                                   * process.zinvBJetsPF
                                   * process.htPFchsNoMuon
                                   * process.mhtPFchsNoMuon
                                   * process.myTopTaggerNoMuon4M
                                   * process.analysisNoMuon4M
                                   * process.myTopTaggerWithMuon4M
                                   * process.analysisWithMuon4M
                                   * process.myTopTaggerNoMuon5M
                                   * process.analysisNoMuon5M
                                   * process.myTopTaggerWithMuon5M
                                   * process.analysisWithMuon5M
                                   * process.myTopTaggerNoMuon6M
                                   * process.analysisNoMuon6M
                                   * process.myTopTaggerWithMuon6M
                                   * process.analysisWithMuon6M
                                   * process.myTopTaggerNoMuon4T
                                   * process.analysisNoMuon4T
                                   * process.myTopTaggerWithMuon4T
                                   * process.analysisWithMuon4T
                                   * process.myTopTaggerNoMuon5T
                                   * process.analysisNoMuon5T
                                   * process.myTopTaggerWithMuon5T
                                   * process.analysisWithMuon5T
                                   * process.myTopTaggerNoMuon6T
                                   * process.analysisNoMuon6T
                                   * process.myTopTaggerWithMuon6T
                                   * process.analysisWithMuon6T
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('dimuonMC400Tree.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')

#process.load('SandBox.Utilities.puWeightProducer_cfi')
##process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/May10_Prompt167151_pudist.root"
#process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/Cert_160404-177515_JSON.pileup.root"

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.p1 = cms.Path(process.eventWeight
                    * process.puWeight
                    * process.analysisSeq )
