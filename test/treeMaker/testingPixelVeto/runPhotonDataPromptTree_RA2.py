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

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "GR_P_V41_AN2::All"

#================= configure poolsource module ===================
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/lpcsusyhad/53X_ntuples/PhotonHad_Run2012C_PromptReco_v2_lpc1/seema/PhotonHad/PhotonHad_Run2012C-PromptReco-v2_NoHepTopTagger_NOCUTS_HLTPhoton70Inc_12Oct2012V3/0cc7c0df13c0d8758c7e8d5139d63072/susypat_777_1_PFV.root',
        '/store/user/lpcsusyhad/53X_ntuples/PhotonHad_Run2012C_PromptReco_v2_lpc1/seema/PhotonHad/PhotonHad_Run2012C-PromptReco-v2_NoHepTopTagger_NOCUTS_HLTPhoton70Inc_12Oct2012V3/0cc7c0df13c0d8758c7e8d5139d63072/susypat_778_1_jUd.root',
        '/store/user/lpcsusyhad/53X_ntuples/PhotonHad_Run2012C_PromptReco_v2_lpc1/seema/PhotonHad/PhotonHad_Run2012C-PromptReco-v2_NoHepTopTagger_NOCUTS_HLTPhoton70Inc_12Oct2012V3/0cc7c0df13c0d8758c7e8d5139d63072/susypat_773_1_lDn.root',
        '/store/user/lpcsusyhad/53X_ntuples/PhotonHad_Run2012C_PromptReco_v2_lpc1/seema/PhotonHad/PhotonHad_Run2012C-PromptReco-v2_NoHepTopTagger_NOCUTS_HLTPhoton70Inc_12Oct2012V3/0cc7c0df13c0d8758c7e8d5139d63072/susypat_774_1_7PV.root',
        '/store/user/lpcsusyhad/53X_ntuples/PhotonHad_Run2012C_PromptReco_v2_lpc1/seema/PhotonHad/PhotonHad_Run2012C-PromptReco-v2_NoHepTopTagger_NOCUTS_HLTPhoton70Inc_12Oct2012V3/0cc7c0df13c0d8758c7e8d5139d63072/susypat_775_1_ZMU.root',
        '/store/user/lpcsusyhad/53X_ntuples/PhotonHad_Run2012C_PromptReco_v2_lpc1/seema/PhotonHad/PhotonHad_Run2012C-PromptReco-v2_NoHepTopTagger_NOCUTS_HLTPhoton70Inc_12Oct2012V3/0cc7c0df13c0d8758c7e8d5139d63072/susypat_781_1_zBt.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
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
from ZInvisibleBkgds.Photons.treemaker_cfi import photonTree
process.analysisID = photonTree.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(True),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(False),

    metSource          = cms.InputTag("pfType1MetNoPhotonID","pfcand"),

    runTopTagger           = cms.bool(True),
    looseTopTaggerSource   = cms.string("photonTopTaggerID5Loose"),
    nominalTopTaggerSource = cms.string("photonTopTaggerID5M"),

    storeExtraVetos    = cms.bool(True),
    muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoPhotonID"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVeto"),
)

process.analysisIDPFIso = process.analysisID.clone(
    PhotonSrc       = cms.InputTag("patPhotonsIDPFIso"),
    JetSrc          = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30"),
    htJetSrc        = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25"),
    bJetSrc         = cms.InputTag("patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    htSource        = cms.InputTag("htPFchsNoPhotIDPFIso"),
    mhtSource       = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
    metSource       = cms.InputTag("pfType1MetNoPhotonIDPFIso","pfcand"),

    looseTopTaggerSource   = cms.string("photonTopTaggerIDPFIso5Loose"),
    nominalTopTaggerSource = cms.string("photonTopTaggerIDPFIso5M"),

    tauVetoSource      = cms.InputTag("sTopTauVetoPhotonIDPFIso"),
)
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')

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
process.load('ZInvisibleBkgds.Photons.ZinvTopTaggers_cff')

process.load('SandBox.Skims.RA2CleaningFilterResults_cfg')
process.load('RecoMET.METFilters.ecalLaserCorrFilter_cfi')
from SandBox.Skims.htFilter_cfi  import *
process.photonIDHTFilter      = htFilter.clone(HTSource = cms.InputTag("htPFchsNoPhotID"),MinHT = cms.double(250))
from SandBox.Skims.mhtFilter_cfi import *
process.photonIDMHTFilter      = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoPhotID"),MinMHT = cms.double(100))

####
process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
                                   * process.htPFchs
                                   * process.mhtPFchs
                                   * process.zinvBJetsPF
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjectsPF
                                   * process.countPhotonsID
                                   * process.photonIDHTFilter
                                   #* process.photonIDMHTFilter
                                   * process.photonMETCollections
                                   * process.photonVetos
                                   * process.photonTopTaggers
                                   * process.ecalLaserCorrFilter
                                   * process.cleaningOnFilterResults
                                   * process.zinvBJetsPFNoPhotonIDSpecial
                                   * process.analysisID
)
#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonDataTreePrompt.root')
)

#============================== configure paths ===============================
process.p1 = cms.Path(process.eventWeight
                    * process.analysisSeq )
