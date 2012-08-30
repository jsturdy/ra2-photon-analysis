import FWCore.ParameterSet.Config as cms

genstudy = cms.EDProducer("GenStudy",
    debug          = cms.bool(False),
    debugString    = cms.string("genstudy"),
    genLabel       = cms.InputTag("genParticles"),
    pdgId          = cms.int32(22),
    genStatus      = cms.int32(3),
    mompdgId       = cms.int32(22),
    momgenStatus   = cms.int32(3),
    genJetLabel    = cms.InputTag("ak5GenJets"),

    recoJetLabel    = cms.InputTag("patJetsAK5PFPt30"),
    recoPhotonLabel = cms.InputTag("patPhotonsIDPFIso"),


    doPUReweight = cms.bool(False),
    puWeight     = cms.InputTag("puWeight"),

    htBins      = cms.vdouble(0.,500., 900., 1300., 2500.),
    mhtBins     = cms.vdouble(0.,200., 350., 500., 1500.),
    nJetBins    = cms.vdouble(0.,2,5,7,9,25),
    minHT       = cms.double(500.),
    minMHT      = cms.double(200.),

    studyAcceptance = cms.bool(False),
    studyRecoIso = cms.bool(False),
    removePhoton    = cms.bool(True),
    bosonPtBins   = cms.vdouble(50., 80., 100., 120., 200., 300., 500., 1000.),
    bosonEtaBins  = cms.vdouble(-5., -2.5, -2.1, -1.566, -1.4442, -0.9, 0.0, 0.9, 1.4442, 1.566, 2.1, 2.5, 5.),
    bosonMinPt    = cms.double(50.),
    bosonEBMaxEta = cms.double(1.4442),
    bosonEEMinEta = cms.double(1.566),
    bosonEEMaxEta = cms.double(2.5)
)
