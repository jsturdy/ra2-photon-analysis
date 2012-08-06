import FWCore.ParameterSet.Config as cms

#photon cuts
#from piet
directpcut       = 'pt>10 && abs(pdgId) == 22 && (abs(mother.pdgId) == 22  && mother.status == 3)'
secondarypcut    = 'pt>10 && abs(pdgId) == 22 && (abs(mother.pdgId) > 100 && mother.status == 2)'
#seema's definitions
directcut        = 'pt>10 && abs(pdgId) == 22 && status == 3 && (abs(mother.pdgId) < 25 || mother.pdgId != 22)'
secondarycut     = 'pt>10 && abs(pdgId) == 22 && status == 1 && abs(mother.pdgId) > 100'
fragmentationcut = 'pt>10 && abs(pdgId) == 22 && status == 1 && (abs(mother.pdgId) < 25 || mother.pdgId != 22)'
mistagcut        = 'pt>10 && abs(pdgId) == 11 && status == 1'

zinvBkgdDirectPhotons = cms.EDFilter("GenParticleSelector",
                                     src = cms.InputTag("genParticles"),
                                     cut = cms.string(directcut)
                                     )

zinvBkgdGenDirectPhotonRefs = cms.EDFilter("GenParticleRefSelector",
                                           src = cms.InputTag("genParticles"),
                                           cut = cms.string(directcut)
                                           )

zinvBkgdSecondaryPhotons = cms.EDFilter("GenParticleSelector",
                                     src = cms.InputTag("genParticles"),
                                     cut = cms.string(secondarycut)
                                     )

zinvBkgdGenSecondaryPhotonRefs = cms.EDFilter("GenParticleRefSelector",
                                           src = cms.InputTag("genParticles"),
                                           cut = cms.string(secondarycut)
                                           )

zinvBkgdFragmentationPhotons = cms.EDFilter("GenParticleSelector",
                                     src = cms.InputTag("genParticles"),
                                     cut = cms.string(fragmentationcut)
                                     )

zinvBkgdGenFragmentationPhotonRefs = cms.EDFilter("GenParticleRefSelector",
                                           src = cms.InputTag("genParticles"),
                                           cut = cms.string(fragmentationcut)
                                           )

zinvBkgdMistagPhotons = cms.EDFilter("GenParticleSelector",
                                     src = cms.InputTag("genParticles"),
                                     cut = cms.string(mistagcut)
                                     )

zinvBkgdGenMistagPhotonRefs = cms.EDFilter("GenParticleRefSelector",
                                           src = cms.InputTag("genParticles"),
                                           cut = cms.string(mistagcut)
                                           )

zinvBkgdGenPhotons = cms.Sequence(
      zinvBkgdGenDirectPhotons   
    * zinvBkgdGenDirectPhotonRefs
    * zinvBkgdGenSecondaryPhotons   
    * zinvBkgdGenSecondaryPhotonRefs
    * zinvBkgdGenFragmentationPhotons   
    * zinvBkgdGenFragmentationPhotonRefs
    * zinvBkgdGenMistagPhotons   
    * zinvBkgdGenMistagPhotonRefs
)

########

zinvBkgdGenElectrons = cms.EDFilter("GenParticleSelector",
                                    src = cms.InputTag("genParticles"),
                                    cut = cms.string('pt>10 && abs(pdgId) == 11')
                                    )

zinvBkgdGenElectronRefs = cms.EDFilter("GenParticleRefSelector",
                                       src = cms.InputTag("genParticles"),
                                       cut = cms.string('pt>10 && abs(pdgId) == 11')
                                       )

zinvBkgdGenElectrons = cms.Sequence(
    zinvBkgdGenElectrons   *
    zinvBkgdGenElectronRefs
)

########

zinvBkgdGenMuons = cms.EDFilter("GenParticleSelector",
                                src = cms.InputTag("genParticles"),
                                cut = cms.string('pt>10 && abs(pdgId) == 13')
                                )

zinvBkgdGenMuonRefs = cms.EDFilter("GenParticleRefSelector",
                                   src = cms.InputTag("genParticles"),
                                   cut = cms.string('pt>10 && abs(pdgId) == 13')
                                   )

zinvBkgdGenMuons = cms.Sequence(
    zinvBkgdGenMuons   *
    zinvBkgdGenMuonRefs
)
