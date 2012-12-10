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

#================= configure poolsource module ===================
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/user/lpcsusyhad/53X_ntuples/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12/gdujany/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_9_1_DrL.root',
        '/store/user/lpcsusyhad/53X_ntuples/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12/gdujany/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_99_1_COa.root',
        '/store/user/lpcsusyhad/53X_ntuples/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6_Summer12/gdujany/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/Summer12_DR53X-PU_S10_V7A-v1_NOCUTS_12Oct2012V3/c31a0db70b73f0b9a355af58227b92dd/susypat_999_1_41G.root',
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
process.source.skipEvents = cms.untracked.uint32(0)
###========================= analysis module =====================================

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
from ZInvisibleBkgds.Photons.treemaker_cfi import photonTree
process.analysisID = photonTree.clone(
    Debug           = cms.bool(False),
    Data            = cms.bool(False),
    ScaleFactor     = cms.double(scaleF),
    DoPUReweight    = cms.bool(True),
    muonVetoSource = cms.InputTag("ra2MuonVeto"),
    electronVetoSource = cms.InputTag("ra2ElectronVeto"),
    metSource       = cms.InputTag("pfType1MetNoPhotonID","pfcandjec"),
)

process.analysisIDPFIso = process.analysisID.clone(
    PhotonSrc       = cms.InputTag("patPhotonsIDPFIso"),
    JetSrc          = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30"),
    htJetSrc        = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25"),
    bJetSrc         = cms.InputTag("patCSVMJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    htSource        = cms.InputTag("htPFchsNoPhotIDPFIso"),
    mhtSource       = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
    metSource       = cms.InputTag("pfType1MetNoPhotonIDPFIso","pfcandjec"),
)
####no vetos applied
process.analysisIDNoVeto      = process.analysisID.clone()
process.analysisIDPFIsoNoVeto = process.analysisIDPFIso.clone()
#================ configure filters and analysis sequence=======================

process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')

process.load('ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff')
from ZInvisibleBkgds.Photons.ZinvBkgdPhotons_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvPhotonJets_cff')
from ZInvisibleBkgds.Photons.ZinvPhotonJets_cff import *
process.load('ZInvisibleBkgds.Photons.ZinvBkgdJets_cff')
process.load('ZInvisibleBkgds.Photons.ZinvBkgdObjects_cff')

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

from ZInvisibleBkgds.Photons.specialMETCleaner_cff import specialPhotonCleanedMET
process.pfType1MetNoPhotonIDPFIso              = specialPhotonCleanedMET.clone(
    inputObjects = cms.InputTag("patPhotonsIDPFIso"),
    inputMET     = cms.InputTag("patMETsPF"),
    inputJets    = cms.InputTag("patJetsPF"),
    inputPFCands = cms.InputTag("particleFlow"),
)
process.pfType1MetNoPhotonID = process.pfType1MetNoPhotonIDPFIso.clone(
    inputObjects = cms.InputTag("patPhotonsID"),
)

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
      
####
process.analysisSeq = cms.Sequence(  process.ra2PFchsJets
                                   * process.htPFchs
                                   * process.mhtPFchs
                                   * process.zinvBJetsPF
                                   * process.rhoToPhotonMap
                                   * process.patPhotonsUser1
                                   * process.patPhotonsUserData
                                   * process.photonObjectsPF
)

process.phoIDSeq = cms.Sequence(  process.countPhotonsID
                                * process.pfType1MetNoPhotonID
                                * process.zinvBJetsPFNoPhotonIDSpecial
                                * process.analysisIDNoVeto
                                * process.ra2MuonVeto
                                * process.ra2ElectronVeto
                                * process.analysisID
)
process.phoIDPFIsoSeq = cms.Sequence( process.countPhotonsIDPFIso
                                    * process.pfType1MetNoPhotonIDPFIso
                                    * process.zinvBJetsPFNoPhotonIDPFIsoSpecial
                                    * process.analysisIDPFIsoNoVeto
                                    * process.ra2MuonVeto
                                    * process.ra2ElectronVeto
                                    * process.analysisIDPFIso
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('photonQCDMC15to3000Tree.root')
)

process.p1 = cms.Path(process.puWeight
                    * process.eventWeight
                    * process.analysisSeq )
process.id    = cms.Path(process.phoIDSeq)
process.idiso = cms.Path(process.phoIDPFIsoSeq)
