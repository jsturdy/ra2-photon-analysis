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
        'file:/tmp/sturdy/DoubleMu_198934-202504_Run2012C-PromptReco-v2_v0.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2500) )
process.source.skipEvents = cms.untracked.uint32(0)

##process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))
##
##process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )
##
##process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)

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
process.analysisNoVeto = dimuonTree.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(True),
    ScaleFactor     = cms.double(1.0),
    DoPUReweight    = cms.bool(False),
    metSource       = cms.InputTag("pfMetType1NoMuon","selected"),
    topTaggerSource = cms.string("myTopTaggerNoMuon"),
)
process.analysisNoVeto4Loose = process.analysisNoVeto.clone(
    topTaggerSource = cms.string("myTopTagger4Loose"),
)
process.analysisNoVeto5Loose = process.analysisNoVeto4Loose.clone(
    topTaggerSource = cms.string("myTopTagger5Loose"),
)
process.analysisNoVeto6Loose = process.analysisNoVeto4Loose.clone(
    topTaggerSource = cms.string("myTopTagger6Loose"),
)
process.analysisNoVeto4M = process.analysisNoVeto.clone(
    topTaggerSource = cms.string("myTopTagger4M"),
)
process.analysisNoVeto5M = process.analysisNoVeto4M.clone(
    topTaggerSource = cms.string("myTopTagger5M"),
)
process.analysisNoVeto6M = process.analysisNoVeto4M.clone(
    topTaggerSource = cms.string("myTopTagger6M"),
)

process.pfMetType1NoMuon = cms.EDProducer('SpecialPATMuonCleanedMETProducer',
    inputMET     = cms.InputTag('pfMetType1'),
    inputObjects = cms.InputTag('specialMuonCollection'),
                                          )

process.analysisVeto4Loose = process.analysisNoVeto4Loose.clone()
process.analysisVeto5Loose = process.analysisNoVeto5Loose.clone()
process.analysisVeto6Loose = process.analysisNoVeto6Loose.clone()
process.analysisVeto4M = process.analysisNoVeto4M.clone()
process.analysisVeto5M = process.analysisNoVeto5M.clone()
process.analysisVeto6M = process.analysisNoVeto6M.clone()
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvMuonJets_cff')
process.load('ZInvisibleBkgds.Photons.zCandFilter_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *
process.load('ZInvisibleBkgds.Photons.specialMuonCollection_cff')
from ZInvisibleBkgds.Photons.specialMuonCollection_cff import specialMuonCollection
process.specialMuonCollection.candidateLabel = cms.InputTag("zToMuMu")
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')

from SandBox.Skims.RA2Objects_cff import countPFMuonsIDIso
process.countPFMuonsIDIsoForZ  = countPFMuonsIDIso.clone(minNumber = cms.uint32(2))

##top tagger
from UserCode.TopTagger.topTagger_cfi import *
process.load("UserCode.TopTagger.topTagger_cfi")
process.myTopTagger = topTagger.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)

process.myTopTagger4Loose = topTagger4Loose.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTagger5Loose = topTagger5Loose.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTagger6Loose = topTagger6Loose.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTagger4M = topTagger4M.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTagger5M = topTagger5M.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)
process.myTopTagger6M = topTagger6M.clone(
    metSrc = cms.InputTag("pfMetType1NoMuon","selected"),
    jetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
)

####
from SandBox.Skims.basicJetSelector_cfi import selectedPatJets
process.patJetsForIndirectTauVeto = selectedPatJets.clone()
process.patJetsForIndirectTauVeto.src = cms.InputTag("patJetsPFNoMuon")
process.patJetsForIndirectTauVeto.cut = cms.string('pt > 15 && abs(eta) < 2.4 && bDiscriminator("combinedSecondaryVertexBJetTags") <= 0.898')

from SandBox.HadSTopSkims.electronVetoFilter_cfi import electronVetoSTop
process.sTopPFElectronVeto = electronVetoSTop.clone(
    DoElectronVeto = cms.bool(True),
    ElectronSource = cms.InputTag('patElectronsPFchs'),
    #ElectronSource = cms.InputTag('gsfElectrons')
)
from SandBox.HadSTopSkims.indirectTauVeto_cfi import indirectTauVeto
process.sTopTauVeto = indirectTauVeto.clone(
    JetSource     = cms.InputTag("patJetsForIndirectTauVeto"),
    METSource     = cms.InputTag("pfMetType1NoMuon","selected"),
    DoTauVeto     = cms.bool(True),
)

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
                                   * process.muonObjectsPF
                                   * process.pfMetType1NoMuon
                                   * process.zinvBJetsPFNoMuon
                                   ####Compute the top-tagger variables
                                   * process.myTopTagger4Loose
                                   * process.myTopTagger5Loose
                                   * process.myTopTagger6Loose
                                   * process.myTopTagger4M
                                   * process.myTopTagger5M
                                   * process.myTopTagger6M

                                   ####Run the analysis with muon-cleaned collections and no lepton vetos
                                   * process.analysisNoVeto4Loose
                                   * process.analysisNoVeto5Loose
                                   * process.analysisNoVeto6Loose
                                   * process.analysisNoVeto4M
                                   * process.analysisNoVeto5M
                                   * process.analysisNoVeto6M

                                   ####Run the analysis with muon-cleaned collections and applying vetos
                                   * process.sTopPFElectronVeto
                                   * process.sTopTauVeto
                                   * process.analysisVeto4Loose
                                   * process.analysisVeto5Loose
                                   * process.analysisVeto6Loose
                                   * process.analysisVeto4M
                                   * process.analysisVeto5M
                                   * process.analysisVeto6M
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zdimuon_DATA_SAMPLENAME_sTop_JOBID.root')
)

#============================== configure paths ===============================
process.p1 = cms.Path(process.eventWeight
                    * process.puWeight
                    * process.analysisSeq )
