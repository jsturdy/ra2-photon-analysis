import FWCore.ParameterSet.Config as cms

process = cms.Process("TreeMaker")

#===================== Message Logger =============================
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
            limit = cms.untracked.int32(10),
            reportEvery = cms.untracked.int32(100)
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 25
process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(True)
            )

process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.GlobalTag.globaltag = "START52_V5::All"
#if runningOnMC == False:
process.GlobalTag.globaltag = "GR_P_V39_AN1::All"

#================= configure poolsource module ===================

###process.load('SusyAnalysis.PhotonAnalysis.PhotonRun2011AMay10ReReco_160404to163869_cfi');
###process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))
##
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/lpcsusyhad/sturdy/PhotonHad/Run2012B-PromptReco-v1/6dd1668cbb6d1245210a6dfb141f845c/susypat_938_1_IeO.root',
        '/store/user/lpcsusyhad/sturdy/PhotonHad/Run2012B-PromptReco-v1/6dd1668cbb6d1245210a6dfb141f845c/susypat_939_1_AvP.root',
        '/store/user/lpcsusyhad/sturdy/PhotonHad/Run2012B-PromptReco-v1/6dd1668cbb6d1245210a6dfb141f845c/susypat_93_1_qwX.root',
        '/store/user/lpcsusyhad/sturdy/PhotonHad/Run2012B-PromptReco-v1/6dd1668cbb6d1245210a6dfb141f845c/susypat_940_1_ORd.root',
        '/store/user/lpcsusyhad/sturdy/PhotonHad/Run2012B-PromptReco-v1/6dd1668cbb6d1245210a6dfb141f845c/susypat_1000_1_l3P.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source.skipEvents = cms.untracked.uint32(0)
####FILELIST = ['file:/eos/uscms/store/user/seema/SusyRA2Analysis2012/25July2012_HTMHT_Run2012B_PromptRecoV3/susypat_1329_1_jK5.root']
####MAXEVENTS = 1000
####SKIPEVENTS = 0
##process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))
##
##process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )
##
##process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)  

#========================= analysis module =====================================

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
    ScaleFactor     = cms.double(1.),
    topTaggerSource = cms.string("myTopTaggerID"),
)
process.analysisID4M = process.analysisID.clone(
    topTaggerSource = cms.string("myTopTaggerID4M"),
)
process.analysisID5M = process.analysisID4M.clone(
    topTaggerSource = cms.string("myTopTaggerID5M"),
)
process.analysisID6M = process.analysisID4M.clone(
    topTaggerSource = cms.string("myTopTaggerID6M"),
)
process.analysisID4T = process.analysisID.clone(
    bJetSrc         = cms.InputTag("patCSVTJetsPFNoPhotonIDSpecialPt30Eta24"),
    topTaggerSource = cms.string("myTopTaggerID4T"),
)
process.analysisID5T = process.analysisID4T.clone(
    topTaggerSource = cms.string("myTopTaggerID5T"),
)
process.analysisID6T = process.analysisID4T.clone(
    topTaggerSource = cms.string("myTopTaggerID6T"),
)
process.analysisIDPFIso = process.analysisID.clone(
    PhotonSrc       = cms.InputTag("patPhotonsIDPFIso"),
    JetSrc          = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30"),
    htJetSrc        = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25"),
    bJetSrc         = cms.InputTag("patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    htSource        = cms.InputTag("htPFchsNoPhotIDPFIso"),
    mhtSource       = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
    topTaggerSource = cms.string("myTopTaggerIDPFIso"),
)
process.analysisIDPFIso4M = process.analysisIDPFIso.clone(
    topTaggerSource = cms.string("myTopTaggerIDPFIso4M"),
)
process.analysisIDPFIso5M = process.analysisIDPFIso4M.clone(
    topTaggerSource = cms.string("myTopTaggerIDPFIso5M"),
)
process.analysisIDPFIso6M = process.analysisIDPFIso4M.clone(
    topTaggerSource = cms.string("myTopTaggerIDPFIso6M"),
)
process.analysisIDPFIso4T = process.analysisIDPFIso.clone(
    bJetSrc         = cms.InputTag("patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    topTaggerSource = cms.string("myTopTaggerIDPFIso4T"),
)
process.analysisIDPFIso5T = process.analysisIDPFIso4T.clone(
    topTaggerSource = cms.string("myTopTaggerIDPFIso5T"),
)
process.analysisIDPFIso6T = process.analysisIDPFIso4T.clone(
    topTaggerSource = cms.string("myTopTaggerIDPFIso6T"),
)
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvPhotonJets_cff')
from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')

process.load('ZInvisibleBkgds.Photons.adduserdata_cfi')

from ZInvisibleBkgds.Photons.photonmap_cfi import *
process.rhoToPhotonMap = photonmap.clone()
from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import *
process.patPhotonsUser1 = addphotonuserdata1.clone()
process.patPhotonsUser1.photonLabel = cms.InputTag("patPhotonsAlt")
process.patPhotonsUser1.userData.userFloats = cms.PSet(
    src = cms.VInputTag(
        cms.InputTag("rhoToPhotonMap")
    )
)
process.patPhotonsUserData = addphotonuserdata2.clone()
process.patPhotonsUserData.photonLabel = cms.InputTag("patPhotonsUser1")

##process.countPhotonsID     = countPhotonsIDPFIso.clone()
##process.countPhotonsIDISO  = countPhotonsIDPFIso.clone()

##top tagger
from UserCode.TopTagger.topTagger_cfi import *
process.load("UserCode.TopTagger.topTagger_cfi")
process.myTopTaggerID        = topTagger.clone()
process.myTopTaggerID.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID4M = topTagger4M.clone()
process.myTopTaggerID4M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID5M = topTagger5M.clone()
process.myTopTaggerID5M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID6M = topTagger6M.clone()
process.myTopTaggerID6M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID4T = topTagger4T.clone()
process.myTopTaggerID4T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID5T = topTagger5T.clone()
process.myTopTaggerID5T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID6T = topTagger6T.clone()
process.myTopTaggerID6T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")

process.myTopTaggerIDPFIso        = topTagger.clone()
process.myTopTaggerIDPFIso.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso4M = topTagger4M.clone()
process.myTopTaggerIDPFIso4M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso5M = topTagger5M.clone()
process.myTopTaggerIDPFIso5M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso6M = topTagger6M.clone()
process.myTopTaggerIDPFIso6M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso4T = topTagger4T.clone()
process.myTopTaggerIDPFIso4T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso5T = topTagger5T.clone()
process.myTopTaggerIDPFIso5T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso6T = topTagger6T.clone()
process.myTopTaggerIDPFIso6T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
      
process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                     process.ra2MuonVeto
                                   * process.ra2ElectronVeto
                                   * process.ra2PFchsJets
                                   * process.zinvBJetsPF
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjectsPF
                                   * process.countPhotonsID
                                   * process.zinvBJetsPFNoPhotonIDSpecial
                                   * process.myTopTaggerID4M
                                   * process.analysisID4M
                                   * process.myTopTaggerID5M
                                   * process.analysisID5M
                                   * process.myTopTaggerID6M
                                   * process.analysisID6M
                                   * process.myTopTaggerID4T
                                   * process.analysisID4T
                                   * process.myTopTaggerID5T
                                   * process.analysisID5T
                                   * process.myTopTaggerID6T
                                   * process.analysisID6T
                                   * process.countPhotonsIDPFIso
                                   * process.zinvBJetsPFNoPhotonIDPFIsoSpecial
                                   * process.myTopTaggerIDPFIso4M
                                   * process.analysisIDPFIso4M
                                   * process.myTopTaggerIDPFIso5M
                                   * process.analysisIDPFIso5M
                                   * process.myTopTaggerIDPFIso6M
                                   * process.analysisIDPFIso6M
                                   * process.myTopTaggerIDPFIso4T
                                   * process.analysisIDPFIso4T
                                   * process.myTopTaggerIDPFIso5T
                                   * process.analysisIDPFIso5T
                                   * process.myTopTaggerIDPFIso6T
                                   * process.analysisIDPFIso6T
)

#======================= output module configuration ===========================
##process.out = cms.OutputModule("PoolOutputModule",
##    fileName = cms.untracked.string('data_userdata.root'),
##    SelectEvents = cms.untracked.PSet(
##        SelectEvents = cms.vstring('p1')
##    ),
##    outputCommands = cms.untracked.vstring('keep *'),
##    dropMetaData = cms.untracked.string('DROPPED')
##)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonDataTree.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')
#process.load('SandBox.Utilities.puWeightProducer_cfi')
#============================== configure paths ===============================
process.p1 = cms.Path(process.puWeight
                    * process.eventWeight
                    * process.analysisSeq )
##process.outpath = cms.EndPath(process.out)
##process.p1 = cms.Path( process.runRangeFilter1 * process.analysisSeq )  #160431 - 161176
##file = open('test_output.py','w')
##file.write(str(process.dumpPython()))
##file.close()
