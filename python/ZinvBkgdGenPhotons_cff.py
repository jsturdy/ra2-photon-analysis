import FWCore.ParameterSet.Config as cms

#photon cuts
#from piet
directpcut       = 'pt>70.0 && abs(pdgId) == 22 && (abs(mother.pdgId) == 22  && mother.status == 3)'
secondarypcut    = 'pt>70.0 && abs(pdgId) == 22 && (abs(mother.pdgId) > 100 && mother.status == 2)'
#seema's definitions
directcut        = 'pt>70.0 && abs(pdgId) == 22 && status == 3 && (abs(mother.pdgId) < 25 || mother.pdgId != 22)'
secondarycut     = 'pt>70.0 && abs(pdgId) == 22 && status == 1 && abs(mother.pdgId) > 100'
fragmentationcut = 'pt>70.0 && abs(pdgId) == 22 && status == 1 && (abs(mother.pdgId) < 25 || mother.pdgId != 22)'
mistagcut        = 'pt>70.0 && abs(pdgId) == 11 && status == 1'
zbosonst1cut     = 'pt>70.0 && abs(pdgId) == 23 && status == 1'
zbosonst3cut     = 'pt>70.0 && abs(pdgId) == 23 && status == 3'

##zbosonst1mu    = 'pt>70.0 && abs(pdgId) == 13 && status == 3 && abs(mother.pdgId) == 23 && mother.status == 1'
##zbosonst3mu    = 'pt>70.0 && abs(pdgId) == 13 && status == 3 && abs(mother.pdgId) == 23 && mother.status == 3'
##zbosonst1el    = 'pt>70.0 && abs(pdgId) == 11 && status == 3 && abs(mother.pdgId) == 23 && mother.status == 1'
##zbosonst3el    = 'pt>70.0 && abs(pdgId) == 11 && status == 3 && abs(mother.pdgId) == 23 && mother.status == 3'
##zbosonst1tau   = 'pt>70.0 && abs(pdgId) == 15 && status == 3 && abs(mother.pdgId) == 23 && mother.status == 1'
##zbosonst3tau   = 'pt>70.0 && abs(pdgId) == 15 && status == 3 && abs(mother.pdgId) == 23 && mother.status == 3'
##zbosonst1nu    = 'pt>70.0 && (abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16) && status == 3 && abs(mother.pdgId) == 23 && mother.status == 1'
##zbosonst3nu    = 'pt>70.0 && (abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16) && status == 3 && abs(mother.pdgId) == 23 && mother.status == 3'

zbosonst1mu    = 'pt>70.0 && abs(pdgId) == 23 && status == 1 && abs(daughter(0).pdgId) == 13'
zbosonst3mu    = 'pt>70.0 && abs(pdgId) == 23 && status == 3 && abs(daughter(0).pdgId) == 13'
zbosonst1el    = 'pt>70.0 && abs(pdgId) == 23 && status == 1 && abs(daughter(0).pdgId) == 11'
zbosonst3el    = 'pt>70.0 && abs(pdgId) == 23 && status == 3 && abs(daughter(0).pdgId) == 11'
zbosonst1tau   = 'pt>70.0 && abs(pdgId) == 23 && status == 1 && abs(daughter(0).pdgId) == 15'
zbosonst3tau   = 'pt>70.0 && abs(pdgId) == 23 && status == 3 && abs(daughter(0).pdgId) == 15'
zbosonst1nu    = 'pt>70.0 && (abs(daughter(0).pdgId) == 12 || abs(daughter(0).pdgId) == 14 || abs(daughter(0).pdgId) == 16) && status == 3 && abs(pdgId) == 23 && status == 1'
zbosonst3nu    = 'pt>70.0 && (abs(daughter(0).pdgId) == 12 || abs(daughter(0).pdgId) == 14 || abs(daughter(0).pdgId) == 16) && status == 3 && abs(pdgId) == 23 && status == 3'

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

##dimuon
zinvBkgdst3ZMuMuBosons = cms.EDFilter("GenParticleSelector",
                                  src = cms.InputTag("zinvBkgdst3ZBosons"),
                                  cut = cms.string(zbosonst3mu)
                                  )

zinvBkgdst3ZMuMuBosonRefs = cms.EDFilter("GenParticleRefSelector",
                                     src = cms.InputTag("zinvBkgdst3ZBosons"),
                                     cut = cms.string(zbosonst3mu)
                                     )

zinvBkgdst1ZMuMuBosons = cms.EDFilter("GenParticleSelector",
                                  src = cms.InputTag("zinvBkgdst1ZBosons"),
                                  cut = cms.string(zbosonst1mu)
                                  )

zinvBkgdst1ZMuMuBosonRefs = cms.EDFilter("GenParticleRefSelector",
                                     src = cms.InputTag("zinvBkgdst1ZBosons"),
                                     cut = cms.string(zbosonst1mu)
                                     )

##dielectron
zinvBkgdst3ZElElBosons = cms.EDFilter("GenParticleSelector",
                                  src = cms.InputTag("zinvBkgdst3ZBosons"),
                                  cut = cms.string(zbosonst3el)
                                  )

zinvBkgdst3ZElElBosonRefs = cms.EDFilter("GenParticleRefSelector",
                                     src = cms.InputTag("zinvBkgdst3ZBosons"),
                                     cut = cms.string(zbosonst3el)
                                     )

zinvBkgdst1ZElElBosons = cms.EDFilter("GenParticleSelector",
                                  src = cms.InputTag("zinvBkgdst1ZBosons"),
                                  cut = cms.string(zbosonst1el)
                                  )

zinvBkgdst1ZElElBosonRefs = cms.EDFilter("GenParticleRefSelector",
                                     src = cms.InputTag("zinvBkgdst1ZBosons"),
                                     cut = cms.string(zbosonst1el)
                                     )
##ditaus
zinvBkgdst3ZTauTauBosons = cms.EDFilter("GenParticleSelector",
                                  src = cms.InputTag("zinvBkgdst3ZBosons"),
                                  cut = cms.string(zbosonst3tau)
                                  )

zinvBkgdst3ZTauTauBosonRefs = cms.EDFilter("GenParticleRefSelector",
                                     src = cms.InputTag("zinvBkgdst3ZBosons"),
                                     cut = cms.string(zbosonst3tau)
                                     )

zinvBkgdst1ZTauTauBosons = cms.EDFilter("GenParticleSelector",
                                  src = cms.InputTag("zinvBkgdst1ZBosons"),
                                  cut = cms.string(zbosonst1tau)
                                  )

zinvBkgdst1ZTauTauBosonRefs = cms.EDFilter("GenParticleRefSelector",
                                     src = cms.InputTag("zinvBkgdst1ZBosons"),
                                     cut = cms.string(zbosonst1tau)
                                     )
##dinu
zinvBkgdst3ZNuNuBosons = cms.EDFilter("GenParticleSelector",
                                  src = cms.InputTag("zinvBkgdst3ZBosons"),
                                  cut = cms.string(zbosonst3nu)
                                  )

zinvBkgdst3ZNuNuBosonRefs = cms.EDFilter("GenParticleRefSelector",
                                     src = cms.InputTag("zinvBkgdst3ZBosons"),
                                     cut = cms.string(zbosonst3nu)
                                     )

zinvBkgdst1ZNuNuBosons = cms.EDFilter("GenParticleSelector",
                                  src = cms.InputTag("zinvBkgdst1ZBosons"),
                                  cut = cms.string(zbosonst1nu)
                                  )

zinvBkgdst1ZNuNuBosonRefs = cms.EDFilter("GenParticleRefSelector",
                                     src = cms.InputTag("zinvBkgdst1ZBosons"),
                                     cut = cms.string(zbosonst1nu)
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

zinvBkgdGenZMuMuBosons = cms.Sequence(
      zinvBkgdst3ZMuMuBosons   
    * zinvBkgdst3ZMuMuBosonRefs
    * zinvBkgdst1ZMuMuBosons   
    * zinvBkgdst1ZMuMuBosonRefs
)
zinvBkgdGenZElElBosons = cms.Sequence(
      zinvBkgdst3ZElElBosons   
    * zinvBkgdst3ZElElBosonRefs
    * zinvBkgdst1ZElElBosons   
    * zinvBkgdst1ZElElBosonRefs
)
zinvBkgdGenZTauTauBosons = cms.Sequence(
      zinvBkgdst3ZTauTauBosons   
    * zinvBkgdst3ZTauTauBosonRefs
    * zinvBkgdst1ZTauTauBosons   
    * zinvBkgdst1ZTauTauBosonRefs
)
zinvBkgdGenZNuNuBosons = cms.Sequence(
      zinvBkgdst3ZNuNuBosons   
    * zinvBkgdst3ZNuNuBosonRefs
    * zinvBkgdst1ZNuNuBosons   
    * zinvBkgdst1ZNuNuBosonRefs
)
