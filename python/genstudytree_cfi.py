import FWCore.ParameterSet.Config as cms

genstudytree = cms.EDProducer("GenStudyTree",
    debug          = cms.bool(False),
    debugString    = cms.string("genstudytree"),
    genSrc       = cms.InputTag("zinvBkgdDirectPhotons"),
    genJetSrc    = cms.InputTag("ak5GenJets"),

    ScaleFactor   = cms.double(1.),
    VertexSrc     = cms.InputTag("goodVertices"),

    recoPhotonSrc = cms.InputTag("patPhotonsIDPFIso"),

    recoJetSrc    = cms.InputTag("patJetsPFNoPhotonSpecialPt30"),
    htJetSrc      = cms.InputTag("patJetsPFNoPhotonSpecialPt50Eta25"),
    bJetSrc       = cms.InputTag("patCSVJetsPFNoPhotonSpecialPt30Eta24"),
    htSource      = cms.InputTag("htPFchs"),
    mhtSource     = cms.InputTag("mhtPFchs"),
    htNoBosonSource  = cms.InputTag("htPFchsNoPhot"),
    mhtNoBosonSource = cms.InputTag("mhtPFchsNoPhot"),


    doPUReweight = cms.bool(True),
    puWeight     = cms.InputTag("puWeight"),

    studyAcceptance = cms.bool(True),
    studyRecoIso    = cms.bool(True),
    bosonMinPt    = cms.double(50.),
    bosonEBMaxEta = cms.double(1.4442),
    bosonEEMinEta = cms.double(1.566),
    bosonEEMaxEta = cms.double(2.5)
)
