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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

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

from ZInvisibleBkgds.Photons.treemaker_cfi import photonTree
process.analysis = photonTree.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(True),
    ScaleFactor     = cms.double(1.),
)

#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvPhotonJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
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

process.countPhotons  = countPhotonsIDPFIso.clone()


process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                     process.ra2PFchsJets
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjectsPF
                                   * process.zinvBJetsPFNoPhotonSpecial
                                   * process.zinvBJetsPF
                                   * process.countPhotons
                                   * process.ra2MuonVeto
                                   * process.ra2ElectronVeto
                                   * process.analysis
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
process.load('SandBox.Utilities.puWeightProducer_cfi')
#============================== configure paths ===============================
process.p1 = cms.Path(process.puWeight * process.analysisSeq )
##process.outpath = cms.EndPath(process.out)
##process.p1 = cms.Path( process.runRangeFilter1 * process.analysisSeq )  #160431 - 161176
##file = open('test_output.py','w')
##file.write(str(process.dumpPython()))
##file.close()
