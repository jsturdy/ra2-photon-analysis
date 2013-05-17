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
process.GlobalTag.globaltag = "FT_P_V42C_AN3::All"

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
process.analysisNoRem = photonTree.clone(
    #Debug           = cms.bool(True),
    Data            = cms.bool(True),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(False),
    PhotonSrc       = cms.InputTag("patPhotonsIDPFIso"),
    TightPhotonSrc  = cms.InputTag("patPhotonsIDPFIso"),

    JetSrc          = cms.InputTag("patJetsPFchsPt30"),
    htJetSrc        = cms.InputTag("patJetsPFchsPt50Eta25"),
    bJetSrc         = cms.InputTag("patCSVMJetsPFPt30Eta24"),
    htSource        = cms.InputTag("htPFchs"),
    mhtSource       = cms.InputTag("mhtPFchs"),
    metSource       = cms.InputTag("newMETwPhiCorr"),
    ra2HTSource     = cms.InputTag("htPFchsNoPhotIDPFIso"),
    ra2MHTSource    = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
)

process.analysisIDPFIso = process.analysisNoRem.clone(
    #Debug           = cms.bool(True),
    DebugString     = cms.string("photonIDPFIso"),
    PhotonSrc       = cms.InputTag("patPhotonsIDPFIso"),
    JetSrc          = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30"),
    htJetSrc        = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25"),
    bJetSrc         = cms.InputTag("patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    htSource        = cms.InputTag("htPFchsNoPhotIDPFIso"),
    mhtSource       = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
    metSource       = cms.InputTag("pfType1MetNoPhotonIDPFIso","pfcand"),
    ra2HTSource     = cms.InputTag("htPFchs"),
    ra2MHTSource    = cms.InputTag("mhtPFchs"),
)
from ZInvisibleBkgds.Photons.templatemaker_cfi import photonTemplate
process.analysisFitTemplate = photonTemplate.clone(
#    ##Debug           = cms.bool(True),
    Data            = cms.bool(True),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(False),
    PhotonSrc       = cms.InputTag("patFitTemplatePhotons"),
    TightPhotonSrc  = cms.InputTag("patPhotonsIDPFIso"),

    JetSrc          = cms.InputTag("patJetsPFNoPhotonFitTemplateSpecialPt30"),
    htJetSrc        = cms.InputTag("patJetsPFNoPhotonFitTemplateSpecialPt50Eta25"),
    htSource        = cms.InputTag("htPFchsNoPhotFitTemplate"),
    mhtSource       = cms.InputTag("mhtPFchsNoPhotFitTemplate"),
    metSource       = cms.InputTag("pfType1MetNoPhotonFitTemplate","pfcand"),
)

process.analysisFakes = process.analysisFitTemplate.clone(
    ##Debug           = cms.bool(True),
    DebugString     = cms.string("photonFakes"),
    PhotonSrc       = cms.InputTag("patJetFakePhotons"),
    JetSrc          = cms.InputTag("patJetsPFNoPhotonJetFakeSpecialPt30"),
    htJetSrc        = cms.InputTag("patJetsPFNoPhotonJetFakeSpecialPt50Eta25"),
    htSource        = cms.InputTag("htPFchsNoPhotJetFake"),
    mhtSource       = cms.InputTag("mhtPFchsNoPhotJetFake"),
    metSource       = cms.InputTag("pfType1MetNoPhotonJetFake","pfcand"),
)
#================ configure filters and analysis sequence=======================

process.load("SandBox.Skims.RA2Leptons_cff")
process.load("SandBox.Skims.jesChange_cfi")
process.newJetsMET.JECLevel = cms.string('ak5PFchsL1FastL2L3Residual')
process.patMETPhiCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_data
                        
process.load('SandBox.Skims.RA2Objects_cff')
process.patJetsPFchsPt30.src      = cms.InputTag('newJetsMET')

process.load('SandBox.Skims.RA2Selection_cff')
from SandBox.Skims.RA2Objects_cff import countJetsPFchsPt50Eta25
process.countJetsPFchsPt50Eta25.src = cms.InputTag('patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25')
process.countJetsPFchsPt50Eta25.minNumber = cms.uint32(2)
process.countJetsPFchsPt50Eta25NoRem = process.countJetsPFchsPt50Eta25.clone()
process.countJetsPFchsPt50Eta25NoRem.src = cms.InputTag('patJetsPFchsPt50Eta25')
process.countJetsPFchsPt50Eta25NoRem.minNumber = cms.uint32(2)

process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvPhotonJets_cff')
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
from RecoMET.METFilters.multiEventFilter_cfi import multiEventFilter
process.hcalLaserEventsPhotonHad = multiEventFilter.clone(
    file = cms.FileInPath('RA2Classic/AdditionalInputFiles/data/HCALLaserEventList_20Nov2012-v2_PhotonHad.txt'))

from SandBox.Skims.jetMHTDPhiFilter_cfi  import *
process.photonDPhiFilter   = jetMHTDPhiFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoPhotFitTemplate"),
                                                    JetSource = cms.InputTag("patJetsPFNoPhotonFitTemplateSpecialPt30"))
from SandBox.Skims.htFilter_cfi  import *
process.photonNoRemHTFilter      = htFilter.clone(HTSource = cms.InputTag("htPFchs"),
                                               MinHT = cms.double(200))
process.photonIDPFIsoHTFilter     = process.photonNoRemHTFilter.clone(HTSource = cms.InputTag("htPFchsNoPhotIDPFIso"))
process.photonFitTemplateHTFilter = process.photonNoRemHTFilter.clone(HTSource = cms.InputTag("htPFchsNoPhotFitTemplate"))
process.photonJetFakeHTFilter     = process.photonNoRemHTFilter.clone(HTSource = cms.InputTag("htPFchsNoPhotJetFake"))

from SandBox.Skims.mhtFilter_cfi import *
process.photonNoRemMHTFilter      = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchs"),
                                                 MinMHT = cms.double(100))
process.photonIDPFIsoMHTFilter     = process.photonNoRemMHTFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoPhotIDPFIso"))
process.photonFitTemplateMHTFilter = process.photonNoRemMHTFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoPhotFitTemplate"))
process.photonJetFakeMHTFilter     = process.photonNoRemMHTFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoPhotJetFake"))

####
process.analysisSeq = cms.Sequence(
    process.ecalLaserCorrFilter
    * process.cleaningOnFilterResults
    * process.hcalLaserEventsPhotonHad
    * process.newra2PFchsJets
    * process.ra2Electrons
    * process.ra2PFchsJets
    * process.htPFchs
    * process.mhtPFchs
    * process.zinvBJetsPF
    * process.rhoToPhotonMap
    * process.patPhotonsUser1
    * process.patPhotonsUserData
    * process.photonObjectsPF
    * process.photonMETCollections
    * process.photonVetos
)
process.rawPhotons = cms.Sequence(
      process.countJetsPFchsPt50Eta25NoRem
    * process.countPhotonsIDPFIso
    * process.photonNoRemHTFilter
    * process.analysisNoRem
    * process.countMaxPhotonsIDPFIso
    * process.photonNoRemMHTFilter
    )
process.idisoPhotons = cms.Sequence(
      process.countJetsPFchsPt50Eta25
    * process.countPhotonsIDPFIso
    * process.photonIDPFIsoHTFilter
    * process.photonIDPFIsoMHTFilter
    * process.zinvBJetsPFNoPhotonIDPFIsoSpecial
    * process.analysisIDPFIso
    * process.countMaxPhotonsIDPFIso
)
process.truePhotons = cms.Sequence(
    process.countFitTemplatePhotons
    * process.photonTemplateObjectsPF
    * process.pfType1MetNoPhotonFitTemplate
    #* process.photonFitTemplateHTFilter
    #* process.photonFitTemplateMHTFilter
    * process.analysisFitTemplate
    * process.countMaxFitTemplatePhotons
)
process.fakePhotons = cms.Sequence(
    process.countJetFakePhotons
    * process.photonJetFakeObjectsPF
    * process.pfType1MetNoPhotonJetFake
    #* process.photonJetFakeHTFilter
    #* process.photonJetFakeMHTFilter
    * process.analysisFakes
    * process.countMaxJetFakePhotons
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonDataTrees.root')
)

#============================== configure paths ===============================
process.raw = cms.Path(  process.eventWeight
                       * process.analysisSeq
                       * process.rawPhotons )
process.idiso = cms.Path(process.eventWeight
                       * process.analysisSeq
                       * process.idisoPhotons )
process.true = cms.Path(process.eventWeight
                      * process.analysisSeq
                      * process.truePhotons )
process.fake = cms.Path(process.eventWeight
                      * process.analysisSeq
                      * process.fakePhotons )

process.mySched = cms.Schedule(process.raw,
                               process.idiso,
                               process.true,
                               process.fake)
