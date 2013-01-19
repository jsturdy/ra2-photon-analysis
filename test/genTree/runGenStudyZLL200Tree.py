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
        '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_1238_1_7hy.root',
        '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_1237_1_EgX.root',
        '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_123_1_YW1.root',
        '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_1239_1_HDA.root',
        '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_1001_1_uzn.root',
        '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_1000_1_WbE.root',
        '/store/user/lpcsusyhad/sturdy/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_1003_1_55n.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source.skipEvents = cms.untracked.uint32(0)

##process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))
##
##process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )
##
##process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)  

#========================= analysis module =====================================

scaleF = 19.73*10*1000/3789889.
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.genstudytree_cfi import gendimuontree

process.zToMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('patMuonsPFIDIso@+ patMuonsPFIDIso@-'),
    cut   = cms.string('50.0 < mass < 120.0'),
    name  = cms.string('zToMuMu'),
    roles = cms.vstring('muon1', 'muon2')
)

process.st3ZMuBosons = gendimuontree.clone(
    genSrc      = cms.InputTag("zinvBkgdst3ZMuMuBosons"),
    debugString = cms.string("status 3 z's from di-muons"),
)
#================ configure filters and analysis sequence=======================

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
process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvMETProducers_cff')
process.load('ZInvisibleBkgds.Photons.ZinvVetos_cff')

process.patJetsPFNoMuon.checkOverlaps.muons.src = cms.InputTag('patMuonsPFIDIso')
process.patJetsPFNoMuonPt30.checkOverlaps.muons.src = cms.InputTag('patMuonsPFIDIso')
process.patJetsPFNoMuonPt50Eta25.checkOverlaps.muons.src = cms.InputTag('patMuonsPFIDIso')

from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import countPhotonsIDPFIso
process.countGenMuMu   = countPhotonsIDPFIso.clone(
    src = cms.InputTag("zinvBkgdst3ZMuMuBosons"),
    minNumber = cms.uint32(1))

process.analysisSeq = cms.Sequence(process.ra2PFchsJets
                                 * process.htPFchs
                                 * process.mhtPFchs
                                 #* process.zinvBkgdst3ZBosons
                                 * process.zinvBkgdst3ZMuMuBosons
                                 * process.zinvBJetsPF
                                 * process.countGenMuMu
                                 * process.patMuonsPFIDIso
                                 #* process.zToMuMu
                                 #* process.zCandFilter
                                 #* process.specialMuonCollection
                                 * process.muonObjectsPF
                                 * process.zmumuMETCollections
                                 * process.zmumuVetos
                                 * process.zinvBJetsPFNoMuon
                                 * process.st3ZMuBosons
)
#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zllht200to400_gen_tree.root')
)
#============================== configure paths ===============================
process.p1 = cms.Path(process.puWeight
                    * process.eventWeight
                    * process.analysisSeq
)

##file = open('wtf_gentree.py','w')
##file.write(str(process.dumpPython()))
##file.close()
