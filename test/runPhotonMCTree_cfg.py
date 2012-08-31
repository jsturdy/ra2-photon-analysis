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

#process.load("Configuration.StandardSequences.Geometry_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.MagneticField_cff")
##process.GlobalTag.globaltag = "START52_V5::All"
##if runningOnMC == False:
#process.GlobalTag.globaltag = "GR_R_52_V9D::All"

#================= configure poolsource module ===================

#process.load('SusyAnalysis.PhotonAnalysis.PhotonRun2011AMay10ReReco_160404to163869_cfi');
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_98_1_YIN.root',
        '/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_99_1_Gav.root',
        '/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_96_1_KJU.root',
        '/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_94_1_9ZP.root',
        '/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_95_1_wIR.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================
scaleF = 107.5*10*1000/1611963.
process.analysis = cms.EDAnalyzer('RA2ZInvPhotonTreeMaker',
                                  Debug           = cms.bool(False),
                                  Data            = cms.bool(False),
                                  ScaleFactor     = cms.double(scaleF),
                                  PhotonSrc       = cms.InputTag("patPhotonsUserData"),
                                  VertexSrc       = cms.InputTag("goodVertices"),
                                  JetSrc          = cms.InputTag("patJetsPFchsPt30"),
                                  bJetSrc         = cms.InputTag("patCSVJetsPFPt30Eta24"),
                                  JetHTSource     = cms.InputTag("patJetsPFchsPt50Eta25"),
                                  DoPUReweight    = cms.bool(True),
                                  PUWeightSource  = cms.InputTag("puWeight"),
)

#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')

from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
process.patPhotonsIDPFIsoLoose  = patPhotonsIDIso.clone(cut = photonISOCutLoose)
process.patPhotonsIDPFIsoMedium = patPhotonsIDIso.clone(cut = photonISOCutMedium)
process.countPhotonsIDPFIsoLoose  = countPhotonsIDPFIso.clone(src = cms.InputTag("patPhotonsIDPFIsoLoose"))
process.countPhotonsIDPFIsoMedium = countPhotonsIDPFIso.clone(src = cms.InputTag("patPhotonsIDPFIsoMedium"))

#process.load('ZInvisibleBkgds.Photons.ZinvPhotonJets_cff')
from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')

from ZInvisibleBkgds.Photons.specialObjectCleaner_cff import specialPhotonCleanedJets
process.patJetsPFSpecial = specialPhotonCleanedJets.clone()
process.patJetsPFSpecial.debug = cms.bool(True)
process.patJetsPFSpecial.objectLabel = cms.InputTag("patPhotonsIDPFIso")

process.load('ZInvisibleBkgds.Photons.addphotonuserdata_cfi')

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

process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                     process.ra2PFchsJets
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjectsPF
                                   * process.patPhotonsIDPFIsoLoose
                                   * process.patPhotonsIDPFIsoMedium
                                   * process.patPhotonsIDPFIso
                                   * process.patJetsPFSpecial
                                   * process.countPhotonsIDPFIsoLoose
                                   * process.countPhotonsIDPFIsoMedium
                                   * process.countPhotonsIDPFIso
                                   * process.zinvBJetsPFNoPhoton
                                   * process.zinvBJetsPF
                                   * process.ra2MuonVeto
                                   * process.ra2ElectronVeto
                                   * process.analysis
)

#======================= output module configuration ===========================
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('mc_userdata.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p1')
    ),
    outputCommands = cms.untracked.vstring('keep *'),
    dropMetaData = cms.untracked.string('DROPPED')
)

                                                                      
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonMCTree.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')
process.load('SandBox.Utilities.puWeightProducer_cfi')
#============================== configure paths ===============================
process.p1 = cms.Path(process.puWeight * process.analysisSeq )
##process.outpath = cms.EndPath(process.out)
##process.p1 = cms.Path( process.runRangeFilter1 * process.analysisSeq )  #160431 - 161176
file = open('wtf_mc.py','w')
file.write(str(process.dumpPython()))
file.close()
