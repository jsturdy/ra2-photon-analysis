import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMaker")

#===================== Message Logger =============================
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
            limit = cms.untracked.int32(10),
            reportEvery = cms.untracked.int32(1)
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 250
process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(False)
            )

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#================= configure poolsource module ===================
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    'file:/uscms_data/d2/sturdy07/SUSY/RA2/newRelease/CMSSW_5_3_5/src/SandBox/Skims/test/susypat_mc_zmumu.root'
        #'/store/user/lpcsusyhad/53X_ntuples/DYJetsToLL_HT_400ToInf_8TeV_madgraph_Summer12_v1_ext/seema/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_ext/DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph_v1_ext_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_101_1_Skk.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source.skipEvents = cms.untracked.uint32(0)
process.GlobalTag.globaltag = "START53_V7G::All"
scaleF =  2.826*10*1000/2727789.

#========================= analysis module =====================================
process.zToMuMu = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('patMuonsPFIDIso@+ patMuonsPFIDIso@-'),
#    decay = cms.string('zinvMuonCandidates@+ zinvMuonCandidates@-'),
    #cut   = cms.string('50.0 < mass < 120.0'),
    cut = cms.string('abs(mass-91.2) < 20.0'), 
    name  = cms.string('zToMuMu'),
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
process.analysisNoRem = dimuonTree.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),

    MuonSrc         = cms.InputTag("patMuonsPFIDIso"),

    JetSrc          = cms.InputTag("patJetsPFchsPt30"),
    htJetSrc        = cms.InputTag("patJetsPFchsPt50Eta25"),
    bJetSrc         = cms.InputTag("patCSVMJetsPFPt30Eta24"),
    htSource        = cms.InputTag("htPFchs"),
    mhtSource       = cms.InputTag("mhtPFchs"),
    metSource       = cms.InputTag("newMETwPhiCorr"),

)
process.analysisIDIso = dimuonTree.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),
)
process.analysisGEN = process.analysisNoRem.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),
    runGenStudy     = cms.bool(True),
    
)
#================ configure filters and analysis sequence=======================

process.load("SandBox.Skims.RA2Leptons_cff")
process.load("SandBox.Skims.jesChange_cfi")
process.newJetsMET.JECLevel = cms.string('ak5PFchsL1FastL2L3')
process.patMETPhiCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc
                        
process.load('SandBox.Skims.RA2Objects_cff')
process.patJetsPFchsPt30.src      = cms.InputTag('newJetsMET')

#process.zinvMuonCandidates = process.patMuonsPFIDIso.clone()
#process.zinvMuonCandidates.MuonSource = cms.InputTag('patMuonsPFIDIso')
#process.zinvMuonCandidates.MinMuPt    = cms.double(20)
#process.zinvMuonCandidates.MinMuEta   = cms.double(2.1)

process.load('SandBox.Skims.RA2Selection_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenLeptons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvMuonJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')

process.load('ZInvisibleBkgds.Photons.zCandFilter_cff')
process.zCandFilter.minPt = cms.double(10.0)
process.zCandPtFilter = process.zCandFilter.clone(minPt = cms.double(100.0))

process.load('ZInvisibleBkgds.Photons.specialMuonCollection_cff')
process.specialMuonCollection.candidateLabel = cms.InputTag("zToMuMu")
process.patMuonsPFID.MuonSource    = cms.InputTag("patMuonsPF")
process.patMuonsPFIDIso.MuonSource = cms.InputTag("patMuonsPF")

from SandBox.Skims.RA2Objects_cff import countJetsPFchsPt50Eta25
process.countJetsPFchsPt50Eta25.src = cms.InputTag('patJetsPFNoMuonPt50Eta25')
process.countJetsPFchsPt50Eta25.minNumber = cms.uint32(2)

from SandBox.Skims.RA2Objects_cff import countPFMuonsIDIso
process.countPFMuonsIDIsoForZ      = countPFMuonsIDIso.clone(
#process.countPFMuonsIDIsoForZ      = countPFMuonsIDIso.clone(src = cms.InputTag("zinvMuonCandidates"),
                                                             minNumber = cms.uint32(2))
process.countPFMuonsIDIsoForZVeto  = countPFMuonsIDIso.clone(minNumber = cms.uint32(2))
process.countPFMuonsIDIsoForZVeto  = process.countPFMuonsIDIsoForZ.clone(maxNumber = cms.uint32(2))

##top tagger

process.load('ZInvisibleBkgds.Photons.ZinvMETProducers_cff')
process.load('ZInvisibleBkgds.Photons.ZinvVetos_cff')

process.load('SandBox.Skims.RA2CleaningFilterResults_cfg')
process.load('SandBox.Skims.RA2CaloVsPFMHTFilterSequence_cff')
process.RA2CaloVsPFMHTFilter.TaggingMode = cms.bool(False)
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')

process.countZsGEN = process.countPhotonsID.clone()
process.countZsGEN.src = cms.InputTag("zinvBkgdst3ZMuMuBosons")
process.countZsGEN.minNumber = cms.uint32(1)

process.countMaxZsGEN = process.countZsGEN.clone()
process.countMaxZsGEN.maxNumber = cms.uint32(1)

from SandBox.Skims.jetMHTDPhiFilter_cfi  import *
process.zmumuDPhiFilter   = jetMHTDPhiFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoMuon"),
                                                   JetSource = cms.InputTag("patJetsPFNoMuonPt30"))
from SandBox.Skims.htFilter_cfi  import *
process.zmumuHTPreFilter   = htFilter.clone(HTSource = cms.InputTag("htPFchsNoMuon"),MinHT = cms.double(200.0))
process.zmumuHTFilter      = htFilter.clone(HTSource = cms.InputTag("htPFchsNoMuon"),MinHT = cms.double(500.0))
from SandBox.Skims.mhtFilter_cfi import *
process.zmumuMHTPreFilter   = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoMuon"),MinMHT = cms.double(100))
process.zmumuMHTFilter      = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoMuon"),MinMHT = cms.double(200))

process.setupSeq = cms.Sequence(  process.newra2PFchsJets
                                * process.ra2Electrons
                                * process.ra2PFchsJets
                                * process.htPFchs
                                * process.mhtPFchs
                                * process.zmumuVetos
                                * process.zinvBJetsPF
                                )
process.cleaningSeq = cms.Sequence(  process.ecalLaserCorrFilter
                                   * process.cleaningOnFilterResults
                                   #* process.RA2CaloVsPFMHTFilterSequence
                                  )
process.analysisSeqIDIso = cms.Sequence(process.countPFMuonsIDIsoForZ
                                   * process.zToMuMu
                                   * process.zCandFilter
                                   * process.specialMuonCollection
                                   * process.muonObjectsPF
                                   * process.zmumuMETCollections
                                   * process.zinvBJetsPFNoMuon
#                                   * process.zmumuHTPreFilter
#                                   * process.zmumuMHTPreFilter
#                                   * process.countJetsPFchsPt50Eta25
                                   * process.analysisIDIso
)

process.analysisSeqNoRem = cms.Sequence(process.countPFMuonsIDIsoForZ
#                                   * process.zToMuMu
#                                   * process.zCandFilter
#                                   * process.specialMuonCollection
#                                   * process.muonObjectsPF
#                                   * process.zmumuMETCollections
#                                   * process.zinvBJetsPFNoMuon
#                                   * process.zmumuHTPreFilter
#                                   * process.zmumuMHTPreFilter
#                                   * process.countJetsPFchsPt50Eta25
                                   * process.analysisNoRem
)

process.analysisSeqGEN = cms.Sequence(  process.zinvBkgdGenZMuMuBosons
                                   * process.zinvBkgdGenMuonSeq
                                   * process.countZsGEN
#                                   * process.zToMuMu
##                                   * process.zCandFilter
#                                   * process.specialMuonCollection
#                                   * process.muonObjectsPF
#                                   * process.zmumuMETCollections
#                                   * process.zinvBJetsPFNoMuon
                                   * process.analysisGEN
)


#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('zmumuMC400Tree.root')
)

#============================== configure paths ===============================
#process.raw = cms.Path(process.eventWeight
#                     * process.puWeight
#                     * process.cleaningSeq
#                     * process.setupSeq
#                     * process.analysisSeqNoRem
#                       )
#process.idiso = cms.Path(process.eventWeight
#                       * process.puWeight
#                       * process.cleaningSeq
#                       * process.setupSeq
#                       * process.analysisSeqIDIso
#                         )
process.gen = cms.Path(process.eventWeight
                     * process.puWeight
                     * process.cleaningSeq
                     * process.setupSeq
                     * process.analysisSeqGEN
                       )
#process.withCaloPFMHT = cms.Path(process.eventWeight
#                               * process.puWeight
#                               * process.cleaningSeq
#                               * process.RA2CaloVsPFMHTFilterSequence
#                               * process.setupSeq
#                               * process.analysisSeq
#                               )
file = open('dimuon400tree.py','w')
file.write(str(process.dumpPython()))
file.close()
