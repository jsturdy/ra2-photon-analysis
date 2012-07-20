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
        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_9_1_F8i.root',
        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_99_1_dpV.root',
        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_98_1_pqg.root',
        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_97_1_7I1.root',
        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_96_1_Jmz.root',
        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_95_1_JiH.root',
        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_94_1_50a.root',
        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_93_1_tdg.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2500) )

process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================

process.analysis = cms.EDAnalyzer('RA2ZInvPhotonAnalyzer',
                                  Debug           = cms.bool(False),
                                  Data            = cms.bool(True),
                                  PhotonSrc       = cms.InputTag("patPhotonsUserData"),
                                  PhotonIsoSrc    = cms.InputTag("patPhotonsIDPFIso"),
                                  JetSrc          = cms.InputTag("patJetsAK5PFPt30"),
                                  JetNoPhotSrc    = cms.InputTag("patJetsAK5PFPt30NoPhotonIDIso"),
                                  #JetNoPhotSrc    = cms.InputTag("patJetsAK5PFPt50Eta25NoPhotonIDIso"),
                                  JetHTSource     = cms.InputTag("patJetsAK5PFPt50Eta25"),
                                  RA2NJets        = cms.uint32(3),
                                  RA2HT           = cms.double(350.0),
                                  RA2MHT          = cms.double(200.0),
                                  RA2ApplyDphiCuts= cms.bool(True),
                                  DoPUReweight    = cms.bool(False),
                                  PUWeightSource  = cms.InputTag("puWeight"),
                                  PFRhoSource     = cms.InputTag("kt6PFJetsForIsolation","rho"),
                                  DoOptimizePlots = cms.bool(True),
                                  GenParticles    = cms.InputTag("genParticles"),
                                  DoGenAnalysis   = cms.bool(False),
)

#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
#from SusyAnalysis.MyAnalysis.filterBoolean_cfi import *
#process.load("SusyAnalysis.MyAnalysis.filterBoolean_cfi")

process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.addphotonuserdata_cfi')
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

from ZInvisibleBkgds.Photons.photonmap_cfi import *
process.rhoToPhotonMap = photonmap.clone()
from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import *
process.patPhotonsUser1 = addphotonuserdata1.clone()
process.patPhotonsUser1.photonLabel = cms.InputTag("patPhotons")
process.patPhotonsUser1.userData.userFloats = cms.PSet(
    src = cms.VInputTag(
        cms.InputTag("rhoToPhotonMap")
    )
)
process.patPhotonsUserData = addphotonuserdata2.clone()
process.patPhotonsUserData.photonLabel = cms.InputTag("patPhotonsUser1")

process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                   process.ra2Objects
                                   * process.kt6PFJetsForIsolation
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjects
                                   * process.ra2PFMuonVeto
                                   * process.ra2PFElectronVeto
                                   * process.analysis
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
   #fileName = cms.string('histPhotonRun2011AMay10ReReco_160404to163869.root')
    fileName = cms.string('photonData.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')

process.load('SandBox.Utilities.puWeightProducer_cfi')
##process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/May10_Prompt167151_pudist.root"
#process.puWeight.DataPileUpHistFile = "SandBox/Utilities/data/Cert_160404-177515_JSON.pileup.root"

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.p1 = cms.Path(process.puWeight * process.analysisSeq )
##process.p1 = cms.Path( process.runRangeFilter1 * process.analysisSeq )  #160431 - 161176
##process.p2 = cms.Path( process.runRangeFilter2 * process.analysisSeq )  #161217 - 163261
##process.p3 = cms.Path( process.runRangeFilter3 * process.analysisSeq )  #163270 - 163869
##process.p4 = cms.Path( process.runRangeFilter4 * process.analysisSeq )  #165088 - 165633
###process.p5 = cms.Path( process.runRangeFilter5 * process.analysisSeq )  #165088 - 165633 
##process.p6 = cms.Path( process.runRangeFilter6 * process.analysisSeq )  #165970 - 166967
##process.p7 = cms.Path( process.runRangeFilter7 * process.analysisSeq )  #167039 - 167784
##process.p8 = cms.Path( process.runRangeFilter8 * process.analysisSeq )  #170249 - 173198
##process.p9 = cms.Path( process.runRangeFilter9 * process.analysisSeq )  #173236 - 177515

