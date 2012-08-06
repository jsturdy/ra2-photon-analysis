import FWCore.ParameterSet.Config as cms

genstudy = cms.EDProducer("GenStudy",
    debug          = cms.bool(False),
    debugString    = cms.string("genstudy"),
    genLabel       = cms.InputTag("genParticles"),
    #genDirectLabel = cms.InputTag("zinvBkgdGenDirectPhotons"),
    #genSecLabel    = cms.InputTag("zinvBkgdGenSecondaryPhotons"),
    #genFragLabel   = cms.InputTag("zinvBkgdGenFragmentationPhotons"),
    #genMisLabel    = cms.InputTag("zinvBkgdGenMistagPhotons"),
    pdgId          = cms.int32(22),
    genStatus      = cms.int32(3),
    mompdgId       = cms.int32(22),
    momgenStatus   = cms.int32(3),
    genJetLabel    = cms.InputTag("ak5GenJets"),
    #genJetLabel    = cms.InputTag("ak5GenJetsNoNu"),


    doPUReweight = cms.bool(False),
    puWeight     = cms.InputTag("puWeight"),

    htBins      = cms.vdouble(300., 500., 800., 1000., 1200., 1400., 1600.),
    mhtBins     = cms.vdouble(100., 200., 350., 500.,  600.,  800.),
    minHT       = cms.double(500.),
    minMHT      = cms.double(200.),

    studyAcceptance = cms.bool(False),
    removePhoton    = cms.bool(True),
    bosonPtBins   = cms.vdouble(50., 80., 100., 120., 200., 300., 500., 1000.),
    bosonEtaBins  = cms.vdouble(-5., -2.5, -2.1, -1.566, -1.4442, -0.9, 0.0, 0.9, 1.4442, 1.566, 2.1, 2.5, 5.),
    bosonMinPt    = cms.double(50.),
    bosonEBMaxEta = cms.double(1.4442),
    bosonEEMinEta = cms.double(1.566),
    bosonEEMaxEta = cms.double(2.5)
)
