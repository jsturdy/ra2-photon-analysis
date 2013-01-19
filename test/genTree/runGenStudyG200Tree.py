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
        'file:/uscms_data/d2/sturdy07/SUSY/RA2/CMSSW_5_3_5/src/ZInvisibleBkgds/Photons/test/genTree/susypat_gjets.root'
        ##'/store/user/lpcsusyhad/sturdy/GJets_HT-200To400_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_784_1_1Ml.root',
        ##'/store/user/lpcsusyhad/sturdy/GJets_HT-200To400_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_785_1_H28.root',
        ##'/store/user/lpcsusyhad/sturdy/GJets_HT-200To400_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_786_1_5oX.root',
        ##'/store/user/lpcsusyhad/sturdy/GJets_HT-200To400_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_787_1_m5s.root',
        ##'/store/user/lpcsusyhad/sturdy/GJets_HT-200To400_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_788_1_8Iw.root',
        ##'/store/user/lpcsusyhad/sturdy/GJets_HT-200To400_8TeV-madgraph/RA2_525_Skims/b5ca2c28b0caa65e44d094fff6785510/susypat_789_1_rUM.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================

scaleF = 960.5*10*1000/10457117.
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.genstudytree_cfi import genphotontree
process.directPhotons = genphotontree.clone(
    debug = cms.bool(False),
    debugString = cms.string("direct photons"),
    ScaleFactor     = cms.double(scaleF),
)
process.secondaryPhotons = process.directPhotons.clone(
    genSrc           = cms.InputTag("zinvBkgdSecondaryPhotons"),
    debugString      = cms.string("secondary photons"))
process.fragmentationPhotons = process.directPhotons.clone(
    genSrc           = cms.InputTag("zinvBkgdFragmentationPhotons"),
    debugString      = cms.string("fragmentation photons"))
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
#from ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff import *
#from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
#from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
#from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *

process.load('ZInvisibleBkgds.Photons.adduserdata_cfi')

from ZInvisibleBkgds.Photons.photonmap_cfi import *
process.rhoToPhotonMap = photonmap.clone()
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

from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import countPhotonsIDPFIso
process.countDirectPhotons        = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdDirectPhotons"))
process.countSecondaryPhotons     = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdSecondaryPhotons"))
process.countFragmentationPhotons = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdFragmentationPhotons"))

process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
                                   * process.htPFchs
                                   * process.mhtPFchs
                                   * process.zinvBJetsPF
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjectsPF
                                   * process.photonMETCollections
                                   * process.photonVetos
                                   * process.zinvBJetsPFNoPhotonIDSpecial
                                   * process.zinvBkgdGenPhotons
)

process.directAnalysisSeq = cms.Sequence(process.countDirectPhotons
                                       * process.directPhotons
)

process.secondaryAnalysisSeq = cms.Sequence(process.countSecondaryPhotons
                                          * process.secondaryPhotons
)

process.fragmentationAnalysisSeq = cms.Sequence(process.countFragmentationPhotons
                                              * process.fragmentationPhotons
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('gjetsht200to400_gen_tree.root')
)

#============================== configure paths ===============================
#process.p1 = cms.Path(process.puWeight
#                    * process.eventWeight
#                    * process.analysisSeq )
process.pdirect = cms.Path(process.puWeight
                         * process.eventWeight
                         * process.analysisSeq 
                         * process.directAnalysisSeq)
process.psecondary = cms.Path(process.puWeight
                            * process.eventWeight
                            * process.analysisSeq 
                            * process.secondaryAnalysisSeq)
process.pfragmentation = cms.Path(process.puWeight
                                * process.eventWeight
                                * process.analysisSeq 
                                * process.fragmentationAnalysisSeq)


##file = open('wtf_gentree.py','w')
##file.write(str(process.dumpPython()))
##file.close()
