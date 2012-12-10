import FWCore.ParameterSet.Config as cms

genstudytree = cms.EDProducer("GenStudyTree",
    debug          = cms.bool(False),
    debugString    = cms.string("genstudytree"),
    genSrc       = cms.InputTag("zinvBkgdDirectPhotons"),
    genJetSrc    = cms.InputTag("ak5GenJets"),
    genMETSrc    = cms.InputTag("genMetTrue"),

    ScaleFactor   = cms.double(1.),
    VertexSrc     = cms.InputTag("goodVertices"),

    recoPhotonSrc = cms.InputTag("patPhotonsIDPFIso"),

    recoJetSrc    = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt30"),
    htJetSrc      = cms.InputTag("patJetsPFNoPhotonIDPFIsoSpecialPt50Eta25"),
    bJetSrc       = cms.InputTag("patCSVTJetsPFNoPhotonIDPFIsoSpecialPt30Eta24"),
    htSource      = cms.InputTag("htPFchs"),
    mhtSource     = cms.InputTag("mhtPFchs"),
    htNoBosonSource  = cms.InputTag("htPFchsNoPhotIDPFIso"),
    mhtNoBosonSource = cms.InputTag("mhtPFchsNoPhotIDPFIso"),
    metSource     = cms.InputTag("patMETsPF"),


    doPUReweight = cms.bool(True),
    puWeights    = cms.InputTag("puWeight"   ,"weight"),
    eventWeights = cms.InputTag("eventWeight","weight"),

    storeExtraVetos    = cms.bool(False),
    electronVetoSource = cms.InputTag("electronVeto"),
    muonVetoSource     = cms.InputTag("muonVeto"),
    tauVetoSource      = cms.InputTag("tauVeto"),
    isoTrkVetoSource   = cms.InputTag("isoTrkVeto"),

    studyAcceptance = cms.bool(True),
    studyRecoIso    = cms.bool(True),
    bosonMinPt    = cms.double(50.),
    bosonEBMaxEta = cms.double(1.4442),
    bosonEEMinEta = cms.double(1.566),
    bosonEEMaxEta = cms.double(2.5)
)
