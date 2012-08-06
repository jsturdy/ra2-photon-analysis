import FWCore.ParameterSet.Config as cms

#photon cuts
#from piet
directpcut       = 'pt>0.0 && abs(pdgId) == 22 && (abs(mother.pdgId) == 22  && mother.status == 3)'
secondarypcut    = 'pt>0.0 && abs(pdgId) == 22 && (abs(mother.pdgId) > 100 && mother.status == 2)'
#seema's definitions
directcut        = 'pt>0.0 && abs(pdgId) == 22 && status == 3 && (abs(mother.pdgId) < 25 || mother.pdgId != 22)'
secondarycut     = 'pt>0.0 && abs(pdgId) == 22 && status == 1 && abs(mother.pdgId) > 100'
fragmentationcut = 'pt>0.0 && abs(pdgId) == 22 && status == 1 && (abs(mother.pdgId) < 25 || mother.pdgId != 22)'
mistagcut        = 'pt>0.0 && abs(pdgId) == 11 && status == 1'
zbosonst1cut     = 'pt>0.0 && abs(pdgId) == 23 && status == 1'
zbosonst3cut     = 'pt>0.0 && abs(pdgId) == 23 && status == 3'

zinvBkgdDirectPhotons = cms.EDFilter("GenParticleSelector",
                                     src = cms.InputTag("genParticles"),
                                     cut = cms.string(directcut)
                                     )

zinvBkgdDirectPhotonRefs = cms.EDFilter("GenParticleRefSelector",
                                        src = cms.InputTag("genParticles"),
                                        cut = cms.string(directcut)
                                        )

zinvBkgdSecondaryPhotons = cms.EDFilter("GenParticleSelector",
                                        src = cms.InputTag("genParticles"),
                                        cut = cms.string(secondarycut)
                                        )

zinvBkgdSecondaryPhotonRefs = cms.EDFilter("GenParticleRefSelector",
                                           src = cms.InputTag("genParticles"),
                                           cut = cms.string(secondarycut)
                                           )

zinvBkgdFragmentationPhotons = cms.EDFilter("GenParticleSelector",
                                            src = cms.InputTag("genParticles"),
                                            cut = cms.string(fragmentationcut)
                                            )

zinvBkgdFragmentationPhotonRefs = cms.EDFilter("GenParticleRefSelector",
                                               src = cms.InputTag("genParticles"),
                                               cut = cms.string(fragmentationcut)
                                               )

zinvBkgdMistagPhotons = cms.EDFilter("GenParticleSelector",
                                     src = cms.InputTag("genParticles"),
                                     cut = cms.string(mistagcut)
                                     )

zinvBkgdMistagPhotonRefs = cms.EDFilter("GenParticleRefSelector",
                                        src = cms.InputTag("genParticles"),
                                        cut = cms.string(mistagcut)
                                        )

zinvBkgdst3ZBosons = cms.EDFilter("GenParticleSelector",
                                  src = cms.InputTag("genParticles"),
                                  cut = cms.string(zbosonst3cut)
                                  )

zinvBkgdst3ZBosonRefs = cms.EDFilter("GenParticleRefSelector",
                                     src = cms.InputTag("genParticles"),
                                     cut = cms.string(zbosonst3cut)
                                     )

zinvBkgdst1ZBosons = cms.EDFilter("GenParticleSelector",
                                  src = cms.InputTag("genParticles"),
                                  cut = cms.string(zbosonst1cut)
                                  )

zinvBkgdst1ZBosonRefs = cms.EDFilter("GenParticleRefSelector",
                                     src = cms.InputTag("genParticles"),
                                     cut = cms.string(zbosonst1cut)
                                     )

zinvBkgdGenPhotons = cms.Sequence(
      zinvBkgdDirectPhotons   
    * zinvBkgdDirectPhotonRefs
    * zinvBkgdSecondaryPhotons   
    * zinvBkgdSecondaryPhotonRefs
    * zinvBkgdFragmentationPhotons   
    * zinvBkgdFragmentationPhotonRefs
    * zinvBkgdMistagPhotons   
    * zinvBkgdMistagPhotonRefs
)
zinvBkgdGenZBosons = cms.Sequence(
      zinvBkgdst3ZBosons   
    * zinvBkgdst3ZBosonRefs
    * zinvBkgdst1ZBosons   
    * zinvBkgdst1ZBosonRefs
)
