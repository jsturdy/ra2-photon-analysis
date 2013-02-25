import FWCore.ParameterSet.Config as cms

process = cms.Process("TemplateMaker")

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

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "FT_P_V42_AN3::All"

#================= configure poolsource module ===================
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012D_PromptReco_v1_lpc1/seema/JetHT/JetHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_1_1_EjF.root',
        '/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012D_PromptReco_v1_lpc1/seema/JetHT/JetHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_2_1_PYU.root',
        '/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012D_PromptReco_v1_lpc1/seema/JetHT/JetHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_3_1_rtj.root',
        '/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012D_PromptReco_v1_lpc1/seema/JetHT/JetHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_4_1_TIa.root',
        '/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012D_PromptReco_v1_lpc1/seema/JetHT/JetHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_5_1_ZL3.root',
        '/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012D_PromptReco_v1_lpc1/seema/JetHT/JetHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_6_1_vfX.root',
        '/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012D_PromptReco_v1_lpc1/seema/JetHT/JetHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_7_1_Kmc.root',
        '/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012D_PromptReco_v1_lpc1/seema/JetHT/JetHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_8_1_p0M.root',
        '/store/user/lpcsusyhad/53X_ntuples/JetHT_Run2012D_PromptReco_v1_lpc1/seema/JetHT/JetHT_Run2012D-PromptReco-v1_NoHepTopTagger_NOCUTS_HLTPFHTInc_12Oct2012V3/b187b7cb1a0ab43ecf500fbb58f6aeb9/susypat_9_1_FY3.root',
        #
        #        '/store/user/lpcsusyhad/lacroix/2012NOV29/MET_D/lacroix/MET/MET_Run2012D-PromptReco-v1_NOCUTS_RA253X_12Oct2012V3/a50d9b1f0a6f357bceb66deb66d50138/susypat_9_1_nGk.root',
        #        '/store/user/lpcsusyhad/lacroix/2012NOV29/MET_D/lacroix/MET/MET_Run2012D-PromptReco-v1_NOCUTS_RA253X_12Oct2012V3/a50d9b1f0a6f357bceb66deb66d50138/susypat_101_1_5sh.root',
        #        '/store/user/lpcsusyhad/lacroix/2012NOV29/MET_D/lacroix/MET/MET_Run2012D-PromptReco-v1_NOCUTS_RA253X_12Oct2012V3/a50d9b1f0a6f357bceb66deb66d50138/susypat_102_1_YvX.root',
        #        '/store/user/lpcsusyhad/lacroix/2012NOV29/MET_D/lacroix/MET/MET_Run2012D-PromptReco-v1_NOCUTS_RA253X_12Oct2012V3/a50d9b1f0a6f357bceb66deb66d50138/susypat_103_2_cwn.root',
        #        '/store/user/lpcsusyhad/lacroix/2012NOV29/MET_D/lacroix/MET/MET_Run2012D-PromptReco-v1_NOCUTS_RA253X_12Oct2012V3/a50d9b1f0a6f357bceb66deb66d50138/susypat_105_1_S87.root',
        #        '/store/user/lpcsusyhad/lacroix/2012NOV29/MET_D/lacroix/MET/MET_Run2012D-PromptReco-v1_NOCUTS_RA253X_12Oct2012V3/a50d9b1f0a6f357bceb66deb66d50138/susypat_106_1_5oF.root',
        #        '/store/user/lpcsusyhad/lacroix/2012NOV29/MET_D/lacroix/MET/MET_Run2012D-PromptReco-v1_NOCUTS_RA253X_12Oct2012V3/a50d9b1f0a6f357bceb66deb66d50138/susypat_107_1_eUe.root',
        #        '/store/user/lpcsusyhad/lacroix/2012NOV29/MET_D/lacroix/MET/MET_Run2012D-PromptReco-v1_NOCUTS_RA253X_12Oct2012V3/a50d9b1f0a6f357bceb66deb66d50138/susypat_100_1_vfD.root',
        #        '/store/user/lpcsusyhad/lacroix/2012NOV29/MET_D/lacroix/MET/MET_Run2012D-PromptReco-v1_NOCUTS_RA253X_12Oct2012V3/a50d9b1f0a6f357bceb66deb66d50138/susypat_110_1_NXG.root',
        #        '/store/user/lpcsusyhad/lacroix/2012NOV29/MET_D/lacroix/MET/MET_Run2012D-PromptReco-v1_NOCUTS_RA253X_12Oct2012V3/a50d9b1f0a6f357bceb66deb66d50138/susypat_111_1_45O.root',
        #        '/store/user/lpcsusyhad/lacroix/2012NOV29/MET_D/lacroix/MET/MET_Run2012D-PromptReco-v1_NOCUTS_RA253X_12Oct2012V3/a50d9b1f0a6f357bceb66deb66d50138/susypat_104_1_E0h.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.source.skipEvents = cms.untracked.uint32(0)
#========================= analysis module =====================================

scaleF = 1.
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(1.0),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.templatemaker_cfi import photonTemplate
process.analysisIDPFIso = photonTemplate.clone(
    ##Debug           = cms.bool(True),
    Data            = cms.bool(True),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(False),
    PhotonSrc       = cms.InputTag("patPhotonsIDPFIsoTight"),
    TightPhotonSrc  = cms.InputTag("patPhotonsIDPFIsoTight"),

    JetSrc          = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30"),
    htJetSrc        = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25"),
    bJetSrc         = cms.InputTag("patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    htSource        = cms.InputTag("htPFchsNoPhotIDPFIso"),
    mhtSource       = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
    metSource       = cms.InputTag("pfType1MetNoPhotonIDPFIso","pfcand"),

    MakeTemplates      = cms.bool(True),
    storeExtraVetos    = cms.bool(True),
    muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoPhotonIDPFIso"),
    isoTrkVetoSource   = cms.InputTag("sTopTrkIsolationMaker","trkIsoVeto"),
)

process.analysisFakes = process.analysisIDPFIso.clone(
    ##Debug           = cms.bool(True),
    DebugString     = cms.string("photonFakes"),
    PhotonSrc       = cms.InputTag("patJetFakePhotons"),
)
#================ configure filters and analysis sequence=======================

process.load("SandBox.Skims.RA2Leptons_cff")
process.load("SandBox.Skims.jesChange_cfi")
process.newJetsMET.JECLevel = cms.string('ak5PFchsL1FastL2L3')
process.patMETPhiCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc
                        
process.load('SandBox.Skims.RA2Objects_cff')
process.patJetsPFchsPt30.src      = cms.InputTag('newJetsMET')

process.load('SandBox.Skims.RA2Selection_cff')
from SandBox.Skims.RA2Objects_cff import countJetsPFchsPt50Eta25
process.countJetsPFchsPt50Eta25.src = cms.InputTag('patJetsPFchsPt50Eta25')
process.countJetsPFchsPt50Eta25.minNumber = cms.uint32(2)

process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvPhotonJets_cff')
from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')

process.load('ZInvisibleBkgds.Photons.addphotonuserdata_cfi')

from ZInvisibleBkgds.Photons.photonmap_cfi import *
process.rhoToPhotonMap = photonmap.clone(photonSrc = cms.InputTag("patPhotonsRA2"))
from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import *
process.patPhotonsUser1 = addphotonuserdata1.clone()
process.patPhotonsUser1.photonLabel = cms.InputTag("patPhotonsRA2")
process.patPhotonsUser1.userData.userFloats = cms.PSet(
    src = cms.VInputTag(
        cms.InputTag("rhoToPhotonMap")
    )
)
process.patPhotonsUserData = addphotonuserdata2.clone()
process.patPhotonsUserData.photonLabel = cms.InputTag("patPhotonsUser1")

process.load('ZInvisibleBkgds.Photons.ZinvMETProducers_cff')
process.load('ZInvisibleBkgds.Photons.ZinvVetos_cff')

process.load('SandBox.Skims.RA2CleaningFilterResults_cfg')
process.load('SandBox.Skims.RA2CaloVsPFMHTFilterSequence_cff')
process.RA2CaloVsPFMHTFilter.TaggingMode = cms.bool(False)
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')

from SandBox.Skims.jetMHTDPhiFilter_cfi  import *
process.photonDPhiFilter   = jetMHTDPhiFilter.clone(MHTSource = cms.InputTag("mhtPFchs"),
                                                  JetSource = cms.InputTag("patJetsPFchsPt30"))
from SandBox.Skims.htFilter_cfi  import *
process.photonIDPFIsoHTFilter      = htFilter.clone(HTSource = cms.InputTag("htPFchsNoPhotIDPFIso"),
                                                    MinHT = cms.double(250))
from SandBox.Skims.mhtFilter_cfi import *
process.photonIDPFIsoMHTFilter      = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
                                                      MinMHT = cms.double(100))

####
process.analysisSeq = cms.Sequence(  process.ecalLaserCorrFilter
                                   * process.cleaningOnFilterResults
                                   * process.newra2PFchsJets
                                   * process.ra2Electrons
                                   * process.countRA2ElectronsIDIso
                                   * process.countRA2MuonsPFIDIso
                                   * process.ra2PFchsJets
                                   * process.htPFchs
                                   * process.mhtPFchs
                                   * process.zinvBJetsPF
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjectsPF
                                   * process.photonMETCollections
                                   #* process.photonIDPFIsoHTFilter
                                   #* process.photonIDPFIsoMHTFilter
                                   * process.photonVetos
                                   * process.zinvBJetsPFNoPhotonIDPFIsoSpecial
)
process.truePhotons = cms.Sequence(  process.countPhotonsIDPFIso
                                   * process.countMaxPhotonsIDPFIso
                                   * process.analysisIDPFIso
)
process.fakePhotons = cms.Sequence(  process.countJetFakePhotons
                                   * process.countMaxJetFakePhotons
                                   * process.analysisFakes
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonDataTemplatePromptJetHT.root')
)

#============================== configure paths ===============================
#process.p1 = cms.Path(process.eventWeight
#                    * process.analysisSeq )
process.true = cms.Path(process.eventWeight
                      * process.analysisSeq
                      * process.truePhotons )
process.fake = cms.Path(process.eventWeight
                      * process.analysisSeq
                      * process.fakePhotons )

process.mySched = cms.Schedule(#process.p1,
                               process.true,
                               process.fake)
