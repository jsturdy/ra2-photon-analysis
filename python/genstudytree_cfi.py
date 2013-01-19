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

    recoJetSrc    = cms.InputTag("patJetsPFchsPt30"),
    htJetSrc      = cms.InputTag("patJetsPFchsPt50Eta25"),
    bJetSrc       = cms.InputTag("patCSVTJetsPFPt30Eta24"),
    htSource      = cms.InputTag("htPFchs"),
    mhtSource     = cms.InputTag("mhtPFchs"),
    metSource     = cms.InputTag("patMETsPF"),
    #htNoBosonSource  = cms.InputTag("htPFchs"),
    #mhtNoBosonSource = cms.InputTag("mhtPFchs"),
    #metNoBosonSource = cms.InputTag("patMETsPF"),


    doPUReweight = cms.bool(True),
    puWeights    = cms.InputTag("puWeight"   ,"weight"),
    eventWeights = cms.InputTag("eventWeight","weight"),

    storeExtraVetos    = cms.bool(False),
    electronVetoSource = cms.InputTag("electronVeto"),
    muonVetoSource     = cms.InputTag("muonVeto"),
    tauVetoSource      = cms.InputTag("tauVeto"),
    isoTrkVetoSource   = cms.InputTag("isoTrkVeto"),

    studyAcceptance = cms.bool(False),
    studyRecoIso    = cms.bool(False),
    nParticles      = cms.int32(2),
    bosonMinPt    = cms.double(50.),
    bosonEBMaxEta = cms.double(1.4442),
    bosonEEMinEta = cms.double(1.566),
    bosonEEMaxEta = cms.double(2.5)
)

genphotontree = genstudytree.clone(
    genSrc       = cms.InputTag("zinvBkgdDirectPhotons"),
    recoPhotonSrc = cms.InputTag("patPhotonsID"),

    recoJetSrc = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30"),
    htJetSrc   = cms.InputTag("patJetsPFNoPhotonIDSpecialPt50Eta25"),
    bJetSrc    = cms.InputTag("patCSVTJetsPFNoPhotonIDSpecialPt30Eta24"),
    htSource   = cms.InputTag("htPFchsNoPhotID"),
    mhtSource  = cms.InputTag("mhtPFchsNoPhotID"),
    metSource  = cms.InputTag("pfType1MetNoPhotonID","pfcand"),
    
    studyAcceptance = cms.bool(True),
    studyRecoIso    = cms.bool(True),
    bosonMinPt    = cms.double(50.),
    nParticles      = cms.int32(1),

    storeExtraVetos    = cms.bool(True),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoPhotonID"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVeto"),
    )

gendimuontree = genstudytree.clone(
    genSrc       = cms.InputTag("zinvBkgdst3ZMuMuBosons"),
    recoPhotonSrc = cms.InputTag("specialMuonCollection"),

    recoJetSrc = cms.InputTag("patJetsPFNoMuonPt30"),
    htJetSrc   = cms.InputTag("patJetsPFNoMuonPt50Eta25"),
    bJetSrc    = cms.InputTag("patCSVMJetsPFNoMuonPt30Eta24"),
    htSource   = cms.InputTag("htPFchsNoMuon"),
    mhtSource  = cms.InputTag("mhtPFchsNoMuon"),
    metSource  = cms.InputTag("pfType1MetNoMuon","pfcand"),
    
    studyAcceptance = cms.bool(False),
    studyRecoIso    = cms.bool(False),
    bosonMinPt    = cms.double(50.),
    nParticles      = cms.int32(2),

    storeExtraVetos    = cms.bool(True),
    muonVetoSource     = cms.InputTag("sTopPFMuonVetoDiMuon"),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoDiMuon"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVetoDiLeptons"),
    )

genzinvtree = genstudytree.clone(
    #genSrc       = cms.InputTag("zinvBkgdst3ZBosons"),
    genSrc       = cms.InputTag("zinvBkgdst3ZNuNuBosons"),
    recoPhotonSrc = cms.InputTag(""),

    storeExtraVetos    = cms.bool(True),
    electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
    muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
    tauVetoSource      = cms.InputTag("sTopTauVetoZInv"),
    isoTrkVetoSource   = cms.InputTag("sTopIsoTrkVeto"),
    )
