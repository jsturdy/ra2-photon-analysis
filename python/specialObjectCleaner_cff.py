import FWCore.ParameterSet.Config as cms

specialMuonCleanedJets = cms.EDProducer("SpecialPATMuonCleanedJetCollection",
    debug          = cms.bool(False),
    debugString    = cms.string("specialMuonCleanedJets"),
    objectLabel    = cms.InputTag("patMuons"),
    jetLabel       = cms.InputTag("patJetsPF"),
    maxDR          = cms.double(0.2),
    arbitration    = cms.bool(False)
)

specialElectronCleanedJets = cms.EDProducer("SpecialPATElectronCleanedJetCollection",
    debug          = cms.bool(False),
    debugString    = cms.string("specialElectronCleanedJets"),
    objectLabel    = cms.InputTag("patElectrons"),
    jetLabel       = cms.InputTag("patJetsPF"),
    maxDR          = cms.double(0.2),
    arbitration    = cms.bool(False)
)

specialPhotonCleanedJets = cms.EDProducer("SpecialPATPhotonCleanedJetCollection",
    debug          = cms.bool(False),
    debugString    = cms.string("specialPhotonCleanedJets"),
    objectLabel    = cms.InputTag("patPhotonsAlt"),
    jetLabel       = cms.InputTag("patJetsPF"),
    maxDR          = cms.double(0.2),
    arbitration    = cms.bool(False)
)

specialTauCleanedJets = cms.EDProducer("SpecialPATTauCleanedJetCollection",
    debug          = cms.bool(False),
    debugString    = cms.string("specialTauCleanedJets"),
    objectLabel    = cms.InputTag("patTaus"),
    jetLabel       = cms.InputTag("patJetsPF"),
    maxDR          = cms.double(0.2),
    arbitration    = cms.bool(False)
)
