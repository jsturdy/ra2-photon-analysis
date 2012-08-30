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
process.GlobalTag.globaltag = "GR_R_52_V9D::All"

#================= configure poolsource module ===================

###process.load('SusyAnalysis.PhotonAnalysis.PhotonRun2011AMay10ReReco_160404to163869_cfi');
###process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))
##
##process.source = cms.Source("PoolSource",
##    fileNames = cms.untracked.vstring(
##        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_9_1_F8i.root',
##        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_99_1_dpV.root',
##        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_98_1_pqg.root',
##        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_97_1_7I1.root',
##        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_96_1_Jmz.root',
##        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_95_1_JiH.root',
##        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_94_1_50a.root',
##        '/store/user/sturdy07/RA2_525_Skims/PhotonHad_2012B_PromptV1_June1_cmslpc/sturdy/PhotonHad/RA2_525_Skims_PhotonHad_2012B_PromptV1_June1_cmslpc/3826414aa9c6dc9311aa462f179e0de0/susypat_data_93_1_tdg.root'
##    )
##)
##
##process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
##
##process.source.skipEvents = cms.untracked.uint32(0)
##FILELIST = ['file:/eos/uscms/store/user/seema/SusyRA2Analysis2012/25July2012_HTMHT_Run2012B_PromptRecoV3/susypat_1329_1_jK5.root']
##MAXEVENTS = 1000
##SKIPEVENTS = 0
process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(MAXEVENTS) )

process.source.skipEvents = cms.untracked.uint32(SKIPEVENTS)  

#========================= analysis module =====================================

process.analysis = cms.EDAnalyzer('RA2ZInvPhotonTreeMaker',
                                  Debug           = cms.bool(False),
                                  Data            = cms.bool(True),
                                  ScaleFactor     = cms.double(1.),
                                  PhotonSrc       = cms.InputTag("patPhotonsUserData"),
                                  VertexSrc       = cms.InputTag("goodVertices"),
                                  JetSrc          = cms.InputTag("patJetsPFPt30"),
                                  bJetSrc         = cms.InputTag("patCSVJetsPFPt30Eta24"),
                                  JetHTSource     = cms.InputTag("patJetsPFPt50Eta25"),
#                                  RA2NJets        = cms.uint32(3),
#                                  RA2HT           = cms.double(350.0),
#                                  RA2MHT          = cms.double(200.0),
#                                  RA2ApplyDphiCuts= cms.bool(True),
                                  DoPUReweight    = cms.bool(False),
                                  PUWeightSource  = cms.InputTag("puWeight"),
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
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.addphotonuserdata_cfi')
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

from ZInvisibleBkgds.Photons.photonmap_cfi import *
process.rhoToPhotonMap = photonmap.clone()
from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import *

##################
#from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFPhotonIso
#process.phoIsoSequence = setupPFPhotonIso(process, 'photons')
#process.phoIsoSequence = setupPFPhotonIso(process, 'patPhotons')
#process.photonPFIsolationDepositsSequencePFIso = cms.Sequence(process.phPFIsoDepositChargedPFIso+
#                                                              process.phPFIsoDepositChargedAllPFIso+
#                                                              process.phPFIsoDepositGammaPFIso+
#                                                              process.phPFIsoDepositNeutralPFIso+
#                                                              process.phPFIsoDepositPUPFIso)
##
##process.phPFIsoValueCharged03PFIdPF    = process.phPFIsoValueCharged03PFId   .clone(src = cms.InputTag("phPFIsoDepositChargedPF"))
##process.phPFIsoValueChargedAll03PFIdPF = process.phPFIsoValueChargedAll03PFId.clone(src = cms.InputTag("phPFIsoDepositChargedAllPF"))
##process.phPFIsoValueGamma03PFIdPF      = process.phPFIsoValueGamma03PFId     .clone(src = cms.InputTag("phPFIsoDepositGammaPF"))
##process.phPFIsoValueNeutral03PFIdPF    = process.phPFIsoValueNeutral03PFId   .clone(src = cms.InputTag("phPFIsoDepositNeutralPF"))
##process.phPFIsoValuePU03PFIdPF         = process.phPFIsoValuePU03PFId        .clone(src = cms.InputTag("phPFIsoDepositPUPF"))
##
##process.phPFIsoValueCharged03PFId.src    = cms.InputTag("phPFIsoDepositCharged")
##process.phPFIsoValueChargedAll03PFId.src = cms.InputTag("phPFIsoDepositChargedAll")
##process.phPFIsoValueGamma03PFId.src      = cms.InputTag("phPFIsoDepositGamma")
##process.phPFIsoValueNeutral03PFId.src    = cms.InputTag("phPFIsoDepositNeutral")
##process.phPFIsoValuePU03PFId.src         = cms.InputTag("phPFIsoDepositPU")
##
##process.pfPhotonIsolationSequence = cms.Sequence(
##    process.phPFIsoValueCharged03PFId+
##    process.phPFIsoValueChargedAll03PFId+
##    process.phPFIsoValueGamma03PFId+
##    process.phPFIsoValueNeutral03PFId+
##    process.phPFIsoValuePU03PFId
##)
##process.pfPhotonIsolationSequencePF = cms.Sequence(
##    process.phPFIsoValueCharged03PFIdPF+
##    process.phPFIsoValueChargedAll03PFIdPF+
##    process.phPFIsoValueGamma03PFIdPF+
##    process.phPFIsoValueNeutral03PFIdPF+
##    process.phPFIsoValuePU03PFIdPF
##)
##
###process.load('PhysicsTools.PatAlgos.producersLayer1.photonProducer_cfi')
###process.patPhotons.userData.userFloats = cms.PSet(
###    src = cms.VInputTag(cms.InputTag("kt6PFJetsForIsolation","rho"))
###)
##process.phoIsoSequence = cms.Sequence(process.pfPhotonIsolationSequence
##                                      +process.pfPhotonIsolationSequencePF)
###process.makePatPhotons = cms.Sequence(process.phoIsoSequence+
###                                      process.patPhotons)
################

process.patPhotonsUser1 = addphotonuserdata1.clone()
process.patPhotonsUser1.photonLabel = cms.InputTag("patPhotonsAlt")
process.patPhotonsUser1.userData.userFloats = cms.PSet(
    src = cms.VInputTag(
        cms.InputTag("rhoToPhotonMap")
    )
)
#process.patPhotonsUser1.isoDeposits = cms.PSet(
#    user = cms.VInputTag(
#        cms.InputTag("phPFIsoDepositCharged"),
#        cms.InputTag("phPFIsoDepositChargedAll"),
#        cms.InputTag("phPFIsoDepositNeutral"),
#        cms.InputTag("phPFIsoDepositGamma"),
#        cms.InputTag("phPFIsoDepositPU")
#    ),
#)
#process.patPhotonsUser1.userIsolation = cms.PSet(
#    user = cms.VPSet(
#        cms.PSet( src = cms.InputTag("phPFIsoValueCharged03PFId")),
#        cms.PSet( src = cms.InputTag("phPFIsoValueChargedAll03PFId")),
#        cms.PSet( src = cms.InputTag("phPFIsoValueNeutral03PFId")),
#        cms.PSet( src = cms.InputTag("phPFIsoValueGamma03PFId")),
#        cms.PSet( src = cms.InputTag("phPFIsoValuePU03PFId"))
#    )
#)
process.patPhotonsUserData = addphotonuserdata2.clone()
process.patPhotonsUserData.photonLabel = cms.InputTag("patPhotonsUser1")

process.analysisSeq = cms.Sequence(#process.ra2PostCleaning   *
                                   process.ra2ObjectsPF
                                   #* process.phoIsoSequence
                                   * process.kt6PFJetsForIsolation
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjectsPF
                                   * process.zinvBJetsPF
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
    fileName = cms.string('photonDataTree_JOBID.root')
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
