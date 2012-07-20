import FWCore.ParameterSet.Config as cms

acceptanceG = cms.EDProducer("AcceptanceG",
    debug          = cms.bool(False),
    debugString    = cms.string("acceptanceG"),
    genLabel       = cms.InputTag("genParticles"),
    genStatus      = cms.int32(3),

    doPUReweight = cms.bool(False),
    puWeight     = cms.InputTag("puWeight"),

    jetHTLabel  = cms.InputTag("patJetsAK5PFPt50Eta25"),
    jetMHTLabel = cms.InputTag("patJetsAK5PFPt30"),
    htLabel     = cms.InputTag("htPF"),
    mhtLabel    = cms.InputTag("mhtPF"),
    htBins      = cms.vdouble(-1., 300., 500., 800., 1000., 1200., 1400., 1600.),
    mhtBins     = cms.vdouble(-1., 100., 200., 350., 500.,  600.,  800.),

    cleanJetHTLabel  = cms.InputTag("patJetsAK5PFPt50Eta25NoPhotonIDIso"),
    cleanJetMHTLabel = cms.InputTag("patJetsAK5PFPt30NoPhotonIDIso"),
    cleanHTLabel     = cms.InputTag("htPFNoPhotIDIso"),
    cleanMHTLabel    = cms.InputTag("mhtPFNoPhotIDIso"),

    photonLabel    = cms.InputTag("patPhotonsIDIso"),
    #isoPhotonLabel    = cms.InputTag("patPhotonsIDIso"),
    photonPtBins   = cms.vdouble(50., 80., 100., 120., 200., 300., 500., 1000.),
    photonEtaBins  = cms.vdouble(-5., -2.5, -2.1, -1.566, -1.4442, -0.9, 0.0, 0.9, 1.4442, 1.566, 2.1, 2.5, 5.),
    photonMinPt    = cms.double(80.),
    photonEBMaxEta = cms.double(1.4442),
    photonEEMinEta = cms.double(1.566),
    photonEEMaxEta = cms.double(2.5)
)
