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
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/tmp/sturdy/GJets_HT-400ToInf_8TeV-madgraph.root'
        #'/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_98_1_YIN.root',
        #'/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_99_1_Gav.root',
        #'/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_96_1_KJU.root',
        #'/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_94_1_9ZP.root',
        #'/store/user/lpcsusyhad/sturdy/GJets_HT-400ToInf_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_95_1_wIR.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source.skipEvents = cms.untracked.uint32(0)

###========================= analysis module =====================================

scaleF = 107.5*10*1000/9539562.
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
    Debug           = cms.bool(False),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),
    metSource       = cms.InputTag("pfType1MetNoPhotonID","selected"),
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
    metSource       = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected"),
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

process.load('ZInvisibleBkgds.Photons.addphotonuserdata_cfi')

from ZInvisibleBkgds.Photons.photonmap_cfi import *
process.rhoToPhotonMap = photonmap.clone(photonSrc = cms.InputTag("patPhotons"))
from ZInvisibleBkgds.Photons.addphotonuserdata_cfi import *
process.patPhotonsUser1 = addphotonuserdata1.clone()
process.patPhotonsUser1.photonLabel = cms.InputTag("patPhotons")
process.patPhotonsUser1.candidateLabel   = cms.InputTag("pfNoPileUpPFchs")
process.patPhotonsUser1.userData.userFloats = cms.PSet(
    src = cms.VInputTag(
        cms.InputTag("rhoToPhotonMap")
    )
)
process.patPhotonsUserData = addphotonuserdata2.clone()
process.patPhotonsUserData.photonLabel = cms.InputTag("patPhotonsUser1")
from ZInvisibleBkgds.Photons.specialMETCleaner_cff import specialPhotonCleanedMET
process.pfType1MetNoPhotonIDPFIso = specialPhotonCleanedMET.clone()
process.pfType1MetNoPhotonIDPFIso.inputObjects = cms.InputTag("patPhotonsIDPFIso")
process.pfType1MetNoPhotonID = specialPhotonCleanedMET.clone()
process.pfType1MetNoPhotonID.inputObjects = cms.InputTag("patPhotonsID")
#process.countPhotonsID     = countPhotonsIDPFIso.clone()
#process.countPhotonsIDISO  = countPhotonsIDPFIso.clone()

##top tagger
from UserCode.TopTagger.topTagger_cfi import *
process.load("UserCode.TopTagger.topTagger_cfi")

process.myTopTaggerID4M = topTagger4M.clone()
process.myTopTaggerID4M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID4M.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
process.myTopTaggerID5M = topTagger5M.clone()
process.myTopTaggerID5M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID5M.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
process.myTopTaggerID6M = topTagger6M.clone()
process.myTopTaggerID6M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID6M.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
process.myTopTaggerID4T = topTagger4T.clone()
process.myTopTaggerID4T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID4T.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
process.myTopTaggerID5T = topTagger5T.clone()
process.myTopTaggerID5T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID5T.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
process.myTopTaggerID6T = topTagger6T.clone()
process.myTopTaggerID6T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30")
process.myTopTaggerID6T.metSrc = cms.InputTag("pfType1MetNoPhotonID","selected")
      
process.myTopTaggerIDPFIso4M = topTagger4M.clone()
process.myTopTaggerIDPFIso4M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso4M.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")
process.myTopTaggerIDPFIso5M = topTagger5M.clone()
process.myTopTaggerIDPFIso5M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso5M.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")
process.myTopTaggerIDPFIso6M = topTagger6M.clone()
process.myTopTaggerIDPFIso6M.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso6M.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")
process.myTopTaggerIDPFIso4T = topTagger4T.clone()
process.myTopTaggerIDPFIso4T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso4T.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")
process.myTopTaggerIDPFIso5T = topTagger5T.clone()
process.myTopTaggerIDPFIso5T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso5T.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")
process.myTopTaggerIDPFIso6T = topTagger6T.clone()
process.myTopTaggerIDPFIso6T.jetSrc = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30")
process.myTopTaggerIDPFIso6T.metSrc = cms.InputTag("pfType1MetNoPhotonIDPFIso","selected")

####
from SandBox.Skims.basicJetSelector_cfi import selectedPatJets
process.patJetsForIndirectTauVeto = selectedPatJets.clone()
process.patJetsForIndirectTauVeto.src = cms.InputTag("patJetsAK5PFchs")
process.patJetsForIndirectTauVeto.cut = cms.string('pt > 15 && abs(eta) < 2.4 && bDiscriminator("combinedSecondaryVertexBJetTags") <= 0.898')

from SandBox.HadSTopSkims.muonVetoFilter_cfi import muonVetoSTop
process.sTopPFMuonVeto = muonVetoSTop.clone()
from SandBox.HadSTopSkims.electronVetoFilter_cfi import electronVetoSTop
process.sTopPFElectronVeto = electronVetoSTop.clone(
    ElectronSource = cms.InputTag('patElectronsPFchs'))
from SandBox.HadSTopSkims.indirectTauVeto_cfi import indirectTauVeto
process.sTopTauVeto = indirectTauVeto.clone()

process.patJetsPFchsPt30.src = cms.InputTag("patJetsAK5PFchs")
process.patCSVMJetsPF.src = cms.InputTag('patJetsAK5PFchs')
process.patCSVTJetsPF.src = cms.InputTag('patJetsAK5PFchs')

process.patJetsPFNoPhotonID.src = cms.InputTag("patJetsAK5PFchs")
process.patJetsPFNoPhotonID.checkOverlaps.taus.src                = cms.InputTag('patTausPFchs')
process.patJetsPFNoPhotonID.checkOverlaps.electrons.src           = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoPhotonID.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoPhotonID.checkOverlaps.muons.src               = cms.InputTag('patMuonsPFchs')

process.patJetsPFNoPhotonIDPFIso.src = cms.InputTag("patJetsAK5PFchs")
process.patJetsPFNoPhotonIDPFIso.checkOverlaps.taus.src                = cms.InputTag('patTausPFchs')
process.patJetsPFNoPhotonIDPFIso.checkOverlaps.electrons.src           = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoPhotonIDPFIso.checkOverlaps.tkIsoElectrons.src      = cms.InputTag('patElectronsPFchs')
process.patJetsPFNoPhotonIDPFIso.checkOverlaps.muons.src               = cms.InputTag('patMuonsPFchs')

process.patJetsPFNoPhotonIDSpecial.jetLabel      = cms.InputTag("patJetsAK5PFchs")
process.patJetsPFNoPhotonIDPFIsoSpecial.jetLabel = cms.InputTag("patJetsAK5PFchs")
####
process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
                                   * process.htPFchs
                                   * process.mhtPFchs
                                   * process.patJetsForIndirectTauVeto
                                   * process.zinvBJetsPF
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjectsPF
                                   * process.pfType1MetNoPhotonID
                                   * process.pfType1MetNoPhotonIDPFIso
)

process.phoIDSeq = cms.Sequence(  process.sTopPFMuonVeto
                                * process.sTopPFElectronVeto
                                * process.sTopTauVeto
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
)

process.phoIDPFIsoSeq = cms.Sequence( process.sTopPFMuonVeto
                                    * process.sTopPFElectronVeto
                                    * process.sTopTauVeto
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
##    fileName = cms.untracked.string('mc_userdata.root'),
##    SelectEvents = cms.untracked.PSet(
##        SelectEvents = cms.vstring('p1')
##    ),
##    outputCommands = cms.untracked.vstring('keep *'),
##    dropMetaData = cms.untracked.string('DROPPED')
##)

                                                                      
process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonMC400Tree.root')
)

#=================== run range & HLT filters ===============================
#process.load('SusyAnalysis.PhotonAnalysis.Photon_RunRangeHLTSeq_cfi')
#process.load('SandBox.Utilities.puWeightProducer_cfi')
#============================== configure paths ===============================
process.p1 = cms.Path(process.eventWeight
                    * process.puWeight
                    * process.analysisSeq )
process.id = cms.Path(process.phoIDSeq)
process.idiso = cms.Path(process.phoIDPFIsoSeq)
##process.outpath = cms.EndPath(process.out)
##process.p1 = cms.Path( process.runRangeFilter1 * process.analysisSeq )  #160431 - 161176
##file = open('wtf_mc.py','w')
##file.write(str(process.dumpPython()))
##file.close()
