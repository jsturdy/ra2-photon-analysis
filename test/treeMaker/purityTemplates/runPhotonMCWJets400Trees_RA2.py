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
        '/store/user/lpcsusyhad/53X_ntuples/WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12/dhare/WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_976_1_DMI.root',
        '/store/user/lpcsusyhad/53X_ntuples/WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12/dhare/WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_975_1_MKr.root',
        '/store/user/lpcsusyhad/53X_ntuples/WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12/dhare/WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_974_1_zri.root',
        '/store/user/lpcsusyhad/53X_ntuples/WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12/dhare/WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_973_1_e82.root',
        '/store/user/lpcsusyhad/53X_ntuples/WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12/dhare/WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_972_1_LE3.root',
        '/store/user/lpcsusyhad/53X_ntuples/WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12/dhare/WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_971_1_8W0.root',
        '/store/user/lpcsusyhad/53X_ntuples/WJetsToLNu_HT-400ToInf_8TeV-madgraph_Summer12/dhare/WJetsToLNu_HT-400ToInf_8TeV-madgraph_v2/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/b9d339f81100b66394e7e5c0a998fe80/susypat_970_1_535.root',
        )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
process.source.skipEvents = cms.untracked.uint32(0)
process.GlobalTag.globaltag = "START53_V7G::All"
###========================= analysis module =====================================

scaleF = 30.08*10*1000/6619651.
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.treemaker_cfi import photonTree
process.analysisNoRem = photonTree.clone(
    #Debug           = cms.bool(True),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),
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
process.countJetsPFchsPt50Eta25NoRem = process.countJetsPFchsPt50Eta25.clone()
process.countJetsPFchsPt50Eta25NoRem.src = cms.InputTag('patJetsPFchsPt50Eta25')
process.countJetsPFchsPt50Eta25NoRem.minNumber = cms.uint32(2)

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

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonMCWJetsTrees.root')
)

process.raw = cms.Path(process.puWeight
                       * process.eventWeight
                       * process.analysisSeq
                       * process.rawPhotons )
process.idiso = cms.Path(process.puWeight
                       * process.eventWeight
                       * process.analysisSeq
                       * process.idisoPhotons
)
