import FWCore.ParameterSet.Config as cms

adduserdata = cms.EDProducer("AddPhotonUserData",
    debug          = cms.bool(False),
    debugString    = cms.string("adduserdata"),
    photonLabel    = cms.InputTag("patPhotons"),
    floatLabels    = cms.VInputTag(cms.InputTag("kt6PFJetsForIsolation","rho")),
    floatNames     = cms.vstring("rho25"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring('pfChargedEA','pfNeutralEA','pfGammaEA','pfChargedRel','pfNeutralRel','pfGammaRel','pfChargedPURel','pfNeutralPURel','pfGammaPURel','combIso','combIsoPU'),
        userFunctions = cms.vstring(
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.002 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.003 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.004 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.006 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.006 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.004 :
                                                         0.003""",
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.024 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.037 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.037 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.034 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.043 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.047 :
                                                         0.066""",
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.053 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.052 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.037 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.073 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.107 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.123 :
                                                         0.133""",
            'userIsolation("User1Iso")/et',
            'userIsolation("User3Iso")/et',
            'userIsolation("User4Iso")/et',
            'max((userIsolation("User1Iso") - userFloat("pfChargedEA")*userFloat("rho25"))/et,0.)',
            'max((userIsolation("User3Iso") - userFloat("pfNeutralEA")*userFloat("rho25"))/et,0.)',
            'max((userIsolation("User4Iso") - userFloat("pfGammaEA")*userFloat("rho25"))/et,0.)',
            'trkSumPtSolidConeDR04 + ecalRecHitSumEtConeDR04 + hcalTowerSumEtConeDR04',
            'trkSumPtSolidConeDR04 + ecalRecHitSumEtConeDR04 + hcalTowerSumEtConeDR04 - 3.141593*0.4*0.4*userFloat("rho25")'
        )
    )
)
