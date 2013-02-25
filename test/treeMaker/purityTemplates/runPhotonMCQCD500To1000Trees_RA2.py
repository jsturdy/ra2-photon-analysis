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

process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

#================= configure poolsource module ===================
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/lpcsusyhad/53X_ntuples/QCD_HT_500To1000_MGPythia_v1_lpc1/seema/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia_v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_982_1_okC.root',
        '/store/user/lpcsusyhad/53X_ntuples/QCD_HT_500To1000_MGPythia_v1_lpc1/seema/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia_v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_984_1_SSH.root',
        '/store/user/lpcsusyhad/53X_ntuples/QCD_HT_500To1000_MGPythia_v1_lpc1/seema/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia6/QCD_HT-500To1000_TuneZ2star_8TeV-madgraph-pythia_v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_983_1_TXN.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )
process.source.skipEvents = cms.untracked.uint32(0)
process.GlobalTag.globaltag = "START53_V7G::All"
###========================= analysis module =====================================

scaleF = 8426.0*10.*1000./30599239.
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.treemaker_cfi import photonTree
process.analysisID = photonTree.clone(
    #Debug           = cms.bool(True),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),
    PhotonSrc       = cms.InputTag("patPhotonsID"),
    TightPhotonSrc  = cms.InputTag("patPhotonsIDPFIso"),

    JetSrc          = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30"),
    htJetSrc        = cms.InputTag("patJetsPFNoPhotonIDSpecialPt50Eta25"),
    htSource        = cms.InputTag("htPFchsNoPhotID"),
    mhtSource       = cms.InputTag("mhtPFchsNoPhotID"),
    metSource       = cms.InputTag("pfType1MetNoPhotonID","pfcand"),
)

process.analysisIDPFIso = process.analysisID.clone(
    #Debug           = cms.bool(True),
    DebugString     = cms.string("photonIDPFIso"),
    PhotonSrc       = cms.InputTag("patPhotonsIDPFIso"),
    JetSrc          = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30"),
    htJetSrc        = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25"),
    bJetSrc         = cms.InputTag("patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    htSource        = cms.InputTag("htPFchsNoPhotIDPFIso"),
    mhtSource       = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
    metSource       = cms.InputTag("pfType1MetNoPhotonIDPFIso","pfcand"),
)
process.analysisGEN = process.analysisID.clone(
    #Debug           = cms.bool(True),
    DebugString     = cms.string("photonGEN"),
    PhotonSrc       = cms.InputTag("patPhotonsID"),
    JetSrc          = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30"),
    htJetSrc        = cms.InputTag("patJetsPFNoPhotonIDSpecialPt50Eta25"),
    htSource        = cms.InputTag("htPFchsNoPhotID"),
    mhtSource       = cms.InputTag("mhtPFchsNoPhotID"),
    metSource       = cms.InputTag("pfType1MetNoPhotonID","pfcand"),
    runGenStudy     = cms.bool(True),
    genSrc          = cms.InputTag("zinvBkgdDirectPhotons")
)
from ZInvisibleBkgds.Photons.templatemaker_cfi import photonTemplate
process.analysisFitTemplate = photonTemplate.clone(
    #Debug           = cms.bool(True),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),
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
process.newJetsMET.JECLevel = cms.string('ak5PFchsL1FastL2L3')
process.patMETPhiCorr.parameter = process.pfMEtSysShiftCorrParameters_2012runABCvsNvtx_mc
                        
process.load('SandBox.Skims.RA2Objects_cff')
process.patJetsPFchsPt30.src      = cms.InputTag('newJetsMET')

process.load('SandBox.Skims.RA2Selection_cff')
from SandBox.Skims.RA2Objects_cff import countJetsPFchsPt50Eta25
process.countJetsPFchsPt50Eta25.src = cms.InputTag('patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25')
process.countJetsPFchsPt50Eta25.minNumber = cms.uint32(2)

process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvPhotonJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')

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
process.photonDPhiFilter   = jetMHTDPhiFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoPhotFitTemplate"),
                                                    JetSource = cms.InputTag("patJetsPFNoPhotonFitTemplateSpecialPt30"))
from SandBox.Skims.htFilter_cfi  import *
process.photonIDHTFilter      = htFilter.clone(HTSource = cms.InputTag("htPFchsNoPhotID"),
                                               MinHT = cms.double(300))
process.photonIDPFIsoHTFilter     = htFilter.clone(HTSource = cms.InputTag("htPFchsNoPhotIDPFIso"))
process.photonFitTemplateHTFilter = htFilter.clone(HTSource = cms.InputTag("htPFchsNoPhotFitTemplate"))
process.photonJetFakeHTFilter     = htFilter.clone(HTSource = cms.InputTag("htPFchsNoPhotJetFake"))

from SandBox.Skims.mhtFilter_cfi import *
process.photonIDMHTFilter      = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoPhotID"),
                                                 MinMHT = cms.double(100))
process.photonIDPFIsoMHTFilter     = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoPhotIDPFIso"))
process.photonFitTemplateMHTFilter = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoPhotFitTemplate"))
process.photonJetFakeMHTFilter     = mhtFilter.clone(MHTSource = cms.InputTag("mhtPFchsNoPhotJetFake"))

####
process.analysisSeq = cms.Sequence(
    process.ecalLaserCorrFilter
    * process.cleaningOnFilterResults
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
process.idisoPhotons = cms.Sequence(
    #  process.countPhotonsID
    #* process.photonIDHTFilter
    #* process.photonIDMHTFilter
    #* process.zinvBJetsPFNoPhotonIDSpecial
    #* process.analysisID
    #* process.countMaxPhotonsID
    #)
    #process.idisoPhotons = cms.Sequence(
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
    #                                   * process.photonFitTemplateHTFilter
    #                                   * process.photonFitTemplateMHTFilter
    * process.analysisFitTemplate
    * process.countMaxFitTemplatePhotons
)
process.fakePhotons = cms.Sequence(
    process.countJetFakePhotons
    * process.photonJetFakeObjectsPF
    * process.pfType1MetNoPhotonJetFake
    #                                   * process.photonJetFakeHTFilter
    #                                   * process.photonJetFakeMHTFilter
    * process.analysisFakes
    * process.countMaxJetFakePhotons
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonMCQCD500Trees.root')
)

process.idiso = cms.Path(process.puWeight
                       * process.eventWeight
                       * process.analysisSeq
                       * process.idisoPhotons )
process.true = cms.Path(process.puWeight
                      * process.eventWeight
                      * process.analysisSeq
                      * process.truePhotons )
process.fake = cms.Path(process.puWeight
                      * process.eventWeight
                      * process.analysisSeq
                      * process.fakePhotons )

process.mySched = cms.Schedule(process.idiso,
                               process.true,
                               process.fake)
