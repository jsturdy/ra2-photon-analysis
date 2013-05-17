import FWCore.ParameterSet.Config as cms

#photon cuts
#from piet
zbosonst1mu    = 'pt>70.0 && abs(pdgId) == 23 && status == 1 && abs(daughter(0).pdgId) == 13'
zbosonst3mu    = 'pt>70.0 && abs(pdgId) == 23 && status == 3 && abs(daughter(0).pdgId) == 13'
zbosonst1el    = 'pt>70.0 && abs(pdgId) == 23 && status == 1 && abs(daughter(0).pdgId) == 11'
zbosonst3el    = 'pt>70.0 && abs(pdgId) == 23 && status == 3 && abs(daughter(0).pdgId) == 11'
zbosonst1tau   = 'pt>70.0 && abs(pdgId) == 23 && status == 1 && abs(daughter(0).pdgId) == 15'
zbosonst3tau   = 'pt>70.0 && abs(pdgId) == 23 && status == 3 && abs(daughter(0).pdgId) == 15'
zbosonst1nu    = 'pt>70.0 && (abs(daughter(0).pdgId) == 12 || abs(daughter(0).pdgId) == 14 || abs(daughter(0).pdgId) == 16) && status == 3 && abs(pdgId) == 23 && status == 1'
zbosonst3nu    = 'pt>70.0 && (abs(daughter(0).pdgId) == 12 || abs(daughter(0).pdgId) == 14 || abs(daughter(0).pdgId) == 16) && status == 3 && abs(pdgId) == 23 && status == 3'
########

zinvBkgdGenElectrons = cms.EDFilter("GenParticleSelector",
                                    src = cms.InputTag("genParticles"),
                                    cut = cms.string('status == 3 && pt>10 && abs(pdgId) == 11')
                                    )

zinvBkgdGenElectronRefs = cms.EDFilter("GenParticleRefSelector",
                                       src = cms.InputTag("genParticles"),
                                       cut = cms.string('status == 3 && pt>10 && abs(pdgId) == 11')
                                       )
zinvBkgdGenElectronDaughters = cms.EDFilter("GenParticleSelector",
                                            src = cms.InputTag("genParticles"),
                                            cut = cms.string('status == 3 && pt>10 && abs(pdgId) == 11')
                                            )

zinvBkgdGenElectronDaughterRefs = cms.EDFilter("GenParticleRefSelector",
                                               src = cms.InputTag("genParticles"),
                                               cut = cms.string('status == 3 && pt>10 && abs(pdgId) == 11')
                                               )

zinvBkgdGenElectronSeq = cms.Sequence(
    zinvBkgdGenElectrons   *
    zinvBkgdGenElectronRefs *
    zinvBkgdGenElectronDaughters   *
    zinvBkgdGenElectronDaughterRefs
)

########

zinvBkgdGenMuons = cms.EDFilter("GenParticleSelector",
                                src = cms.InputTag("genParticles"),
                                cut = cms.string('status == 3 && pt>10 && abs(pdgId) == 13')
                                )

zinvBkgdGenMuonRefs = cms.EDFilter("GenParticleRefSelector",
                                   src = cms.InputTag("genParticles"),
                                   cut = cms.string('status == 3 && pt>10 && abs(pdgId) == 13')
                                   )

zinvBkgdGenMuonDaughters = cms.EDFilter("GenParticleSelector",
                                src = cms.InputTag("genParticles"),
                                cut = cms.string('status == 3 && pt>10 && abs(pdgId) == 13')
                                )

zinvBkgdGenMuonDaughterRefs = cms.EDFilter("GenParticleRefSelector",
                                   src = cms.InputTag("genParticles"),
                                   cut = cms.string('status == 3 && pt>10 && abs(pdgId) == 13')
                                   )

zinvBkgdGenMuonSeq = cms.Sequence(
    zinvBkgdGenMuons   *
    zinvBkgdGenMuonRefs *
    zinvBkgdGenMuonDaughters   *
    zinvBkgdGenMuonDaughterRefs
)
