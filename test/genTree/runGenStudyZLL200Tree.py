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
        '/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_200To400_8TeV_madgraph_Summer12_v1/seema/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_v1_NoHepTopTagger_NOCUTS_12Oct2012V3/67e73f866d34072a7aee49050064a888/susypat_96_1_l5H.root',
        '/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_200To400_8TeV_madgraph_Summer12_v1/seema/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_v1_NoHepTopTagger_NOCUTS_12Oct2012V3/67e73f866d34072a7aee49050064a888/susypat_98_1_qVn.root',
        '/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_200To400_8TeV_madgraph_Summer12_v1/seema/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_v1_NoHepTopTagger_NOCUTS_12Oct2012V3/67e73f866d34072a7aee49050064a888/susypat_99_1_eKZ.root',
        '/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_200To400_8TeV_madgraph_Summer12_v1/seema/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph/DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph_v1_NoHepTopTagger_NOCUTS_12Oct2012V3/67e73f866d34072a7aee49050064a888/susypat_9_1_tqK.root',
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

process.analysisSeq = cms.Sequence(process.zinvBkgdst3ZMuMuBosons
                                 #* process.zinvBkgdst3ZBosons
                                 * process.countGenMuMu
                                 * process.ra2PFchsJets
                                 * process.htPFchs
                                 * process.mhtPFchs
                                 * process.zinvBJetsPF
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
