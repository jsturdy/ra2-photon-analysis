import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

#===================== Message Logger =============================
process.load("FWCore.MessageService.MessageLogger_cfi")

process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
            limit = cms.untracked.int32(10),
            reportEvery = cms.untracked.int32(100)
            )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options = cms.untracked.PSet(
            wantSummary = cms.untracked.bool(True)
            )

#================= configure poolsource module ===================

#process.load('SusyAnalysis.PhotonAnalysis.PhotonRun2011AMay10ReReco_160404to163869_cfi');
#process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(FILELIST ))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/lpcsusyhad/53X_ntuples/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12/gdujany/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_9_1_DrL.root',
        '/store/user/lpcsusyhad/53X_ntuples/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12/gdujany/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_99_1_COa.root',
        '/store/user/lpcsusyhad/53X_ntuples/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12/gdujany/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_999_1_41G.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.source.skipEvents = cms.untracked.uint32(0)

#========================= analysis module =====================================

scaleF = -1.0
from RA2Classic.WeightProducer.weightProducer_cfi import weightProducer
process.eventWeight = weightProducer.clone(
    weight = cms.double(scaleF),
    Method = cms.string("PtHat"),
    XS = cms.double(2.99913994E10),
    NumberEvts = cms.double(9991647),
    Lumi = cms.double(10000.),
    Exponent = cms.double(-4.5),
    LumiScale = cms.double( 1.0 ),
)
from RA2Classic.WeightProducer.puWeightProducer_cfi import puWeightProducer
process.puWeight = puWeightProducer.clone(
    weight = cms.double(1.0),
)
from ZInvisibleBkgds.Photons.genstudytree_cfi import *
process.directPhotonsID = genstudytree.clone(
    debug            = cms.bool(False),
    genSrc           = cms.InputTag("zinvBkgdDirectPhotons"),
    debugString      = cms.string("direct photons"),
    ScaleFactor      = cms.double(scaleF),
    recoPhotonSrc    = cms.InputTag("patPhotonsID"),
    recoJetSrc       = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30"),
    htJetSrc         = cms.InputTag("patJetsPFNoPhotonIDSpecialPt50Eta25"),
    bJetSrc          = cms.InputTag("patCSVTJetsPFNoPhotonIDSpecialPt30Eta24"),
    htNoBosonSource  = cms.InputTag("htPFchsNoPhotID"),
    mhtNoBosonSource = cms.InputTag("mhtPFchsNoPhotID"),
)
process.directPhotonsIDPFIso = process.directPhotonsID.clone(
    recoPhotonSrc    = cms.InputTag("patPhotonsIDPFIso"),
    recoJetSrc       = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30"),
    htJetSrc         = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25"),
    bJetSrc          = cms.InputTag("patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    htNoBosonSource  = cms.InputTag("htPFchsNoPhotIDPFIso"),
    mhtNoBosonSource = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
)
process.directPhotonsIDNoVeto = process.directPhotonsID.clone()
process.directPhotonsIDPFIsoNoVeto = process.directPhotonsIDPFIso.clone()
process.secondaryPhotonsID = process.directPhotonsID.clone(
    genSrc           = cms.InputTag("zinvBkgdSecondaryPhotons"),
    debugString      = cms.string("secondary photons"))
process.secondaryPhotonsIDPFIso = process.directPhotonsIDPFIso.clone(
    genSrc           = cms.InputTag("zinvBkgdSecondaryPhotons"),
    debugString      = cms.string("secondary photons"))
process.secondaryPhotonsIDNoVeto = process.secondaryPhotonsID.clone()
process.secondaryPhotonsIDPFIsoNoVeto = process.secondaryPhotonsIDPFIso.clone()
process.fragmentationPhotonsID = process.directPhotonsID.clone(
    genSrc           = cms.InputTag("zinvBkgdFragmentationPhotons"),
    debugString      = cms.string("fragmentation photons"))
process.fragmentationPhotonsIDPFIso = process.directPhotonsIDPFIso.clone(
    genSrc           = cms.InputTag("zinvBkgdFragmentationPhotons"),
    debugString      = cms.string("fragmentation photons"))
process.fragmentationPhotonsIDNoVeto = process.fragmentationPhotonsID.clone()
process.fragmentationPhotonsIDPFIsoNoVeto = process.fragmentationPhotonsIDPFIso.clone()
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')

process.load('ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdGenPhotons_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
from ZInvisibleBkgds.Photons.ZinvBkgdJets_cff import *

process.load('ZInvisibleBkgds.Photons.adduserdata_cfi')

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

process.countDirectPhotons        = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdDirectPhotons"))
process.countSecondaryPhotons     = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdSecondaryPhotons"))
process.countFragmentationPhotons = countPhotonsIDPFIso.clone(src = cms.InputTag("zinvBkgdFragmentationPhotons"))

process.patJetsPFNoPhotonID.checkOverlaps.electrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonID.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonID.checkOverlaps.muons.src = cms.InputTag('patMuonsPF')
process.patJetsPFNoPhotonID.checkOverlaps.taus.src = cms.InputTag('selectedPatTausPF')

process.patJetsPFNoPhotonIDPt30.checkOverlaps.electrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonIDPt30.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonIDPt30.checkOverlaps.muons.src = cms.InputTag('patMuonsPF')
process.patJetsPFNoPhotonIDPt30.checkOverlaps.taus.src = cms.InputTag('selectedPatTausPF')

process.patJetsPFNoPhotonIDPt50Eta25.checkOverlaps.electrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonIDPt50Eta25.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonIDPt50Eta25.checkOverlaps.muons.src = cms.InputTag('patMuonsPF')
process.patJetsPFNoPhotonIDPt50Eta25.checkOverlaps.taus.src = cms.InputTag('selectedPatTausPF')

process.patJetsPFNoPhotonIDPFIso.checkOverlaps.electrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonIDPFIso.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonIDPFIso.checkOverlaps.muons.src = cms.InputTag('patMuonsPF')
process.patJetsPFNoPhotonIDPFIso.checkOverlaps.taus.src = cms.InputTag('selectedPatTausPF')

process.patJetsPFNoPhotonIDPFIsoPt30.checkOverlaps.electrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonIDPFIsoPt30.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonIDPFIsoPt30.checkOverlaps.muons.src = cms.InputTag('patMuonsPF')
process.patJetsPFNoPhotonIDPFIsoPt30.checkOverlaps.taus.src = cms.InputTag('selectedPatTausPF')

process.patJetsPFNoPhotonIDPFIsoPt50Eta25.checkOverlaps.electrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonIDPFIsoPt50Eta25.checkOverlaps.tkIsoElectrons.src = cms.InputTag('patElectrons')
process.patJetsPFNoPhotonIDPFIsoPt50Eta25.checkOverlaps.muons.src = cms.InputTag('patMuonsPF')
process.patJetsPFNoPhotonIDPFIsoPt50Eta25.checkOverlaps.taus.src = cms.InputTag('selectedPatTausPF')

process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.zinvBkgdGenPhotons
#                                   * process.zinvBkgdGenZBosons
                                   * process.photonObjectsPF
                                   * process.zinvBJetsPFNoPhotonIDSpecial
                                   * process.zinvBJetsPFNoPhotonIDPFIsoSpecial
                                   * process.zinvBJetsPF
#                                   * process.countGenBosons
#                                   * process.directPhotonsIDNoVeto
#                                   * process.directPhotonsIDPFIsoNoVeto
#                                   * process.ra2MuonVeto
#                                   * process.ra2ElectronVeto
#                                   * process.directPhotonsID
#                                   * process.directPhotonsIDPFIso
)

process.directAnalysisSeq = cms.Sequence(process.countDirectPhotons
                                       * process.directPhotonsIDNoVeto
                                       * process.directPhotonsIDPFIsoNoVeto
                                       * process.ra2MuonVeto
                                       * process.ra2ElectronVeto
                                       * process.directPhotonsID
                                       * process.directPhotonsIDPFIso
)

process.secondaryAnalysisSeq = cms.Sequence(process.countSecondaryPhotons
                                          * process.secondaryPhotonsIDNoVeto
                                          * process.secondaryPhotonsIDPFIsoNoVeto
                                          * process.ra2MuonVeto
                                          * process.ra2ElectronVeto
                                          * process.secondaryPhotonsID
                                          * process.secondaryPhotonsIDPFIso
)

process.fragmentationAnalysisSeq = cms.Sequence(process.countFragmentationPhotons
                                              * process.fragmentationPhotonsIDNoVeto
                                              * process.fragmentationPhotonsIDPFIsoNoVeto
                                              * process.ra2MuonVeto
                                              * process.ra2ElectronVeto
                                              * process.fragmentationPhotonsID
                                              * process.fragmentationPhotonsIDPFIso
)

#======================= output module configuration ===========================

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('qcd_flat_pt_15to3000_gen_tree.root')
)

#============================== configure paths ===============================
#process.p1 = cms.Path( process.analysisSeq )
process.p1 = cms.Path(process.puWeight
                    * process.eventWeight
                    * process.analysisSeq )
process.pdirect        = cms.Path(process.directAnalysisSeq)
process.psecondary     = cms.Path(process.secondaryAnalysisSeq)
process.pfragmentation = cms.Path(process.fragmentationAnalysisSeq)
