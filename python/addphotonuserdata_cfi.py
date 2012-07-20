import FWCore.ParameterSet.Config as cms

addphotonuserdata1 = cms.EDProducer("AddPhotonUserData",
    debug          = cms.bool(False),
    debugString    = cms.string("addphotonuserdata"),
    photonLabel    = cms.InputTag("patPhotons"),
    floatLabels    = cms.VInputTag(cms.InputTag("kt6PFJetsForIsolation","rho")),
    floatNames     = cms.vstring("rho25"),
    embedConversionInfo = cms.bool(True),
    gsfElectronLabel = cms.InputTag("gsfElectrons"),
    conversionsLabel = cms.InputTag("conversions"),
    beamspotLabel    = cms.InputTag("offlineBeamSpot"),
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
        userFunctionLabels = cms.vstring('hcalIsoConeDR03_2012','hcalIsoConeDR04_2012',
                                         'pfChargedEA','pfNeutralEA','pfGammaEA',
                                         'hadTowOverEmTightCut' ,'showerShapeTightCut' ,'pfChargedTightCut' ,'pfNeutralTightCut' ,'pfGammaTightCut',
                                         'hadTowOverEmMediumCut','showerShapeMediumCut','pfChargedMediumCut','pfNeutralMediumCut','pfGammaMediumCut',
                                         'hadTowOverEmLooseCut' ,'showerShapeLooseCut' ,'pfChargedLooseCut' ,'pfNeutralLooseCut' ,'pfGammaLooseCut',
                                         'pfChargedRel','pfNeutralRel','pfGammaRel',
                                         'combIsoR03','combIsoR04'),
        userFunctions = cms.vstring(
            """hcalTowerSumEtConeDR03 +
               (hadronicOverEm - hadTowOverEm)*superCluster.energy/cosh(superCluster.eta)""", #hcalIsoConeDR03_2012
            """hcalTowerSumEtConeDR04 +
               (hadronicOverEm - hadTowOverEm)*superCluster.energy/cosh(superCluster.eta)""", #hcalIsoConeDR04_2012
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.002 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.003 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.004 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.006 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.006 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.004 :
                                                         0.003""", #chargedEA
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.024 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.037 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.037 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.034 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.043 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.047 :
                                                         0.066""", #neutralEA
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.053 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.052 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.037 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.073 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.107 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.123 :
                                                         0.133""", #gammaEA

            ####Cut values:singleTowerHOverE,showerShape,pfChargedRel,pfNeutralRel,pfGammaRel
            ###Tight
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05 """,#singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.011:0.030""",#showerShape
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.02 :0.02 """,#pfChargedRel
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05 """,#pfNeutralRel
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.09 :0.09 """,#pfGammaRel
            ###Medium
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05 """,#singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.011:0.031""",#showerShape      
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.03 :0.03 """,#pfChargedRel     
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.07 :0.07 """,#pfNeutralRel     
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.08 :0.09 """,#pfGammaRel       
            ###Loose
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.06 :0.05 """,#singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.011:0.034""",#showerShape      
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.06 :0.05 """,#pfChargedRel     
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.16 :0.10 """,#pfNeutralRel     
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.08 :0.12 """,#pfGammaRel       

            'userIsolation("User1Iso")/et',
            'userIsolation("User3Iso")/et',
            'userIsolation("User4Iso")/et',
            'trkSumPtSolidConeDR03 + ecalRecHitSumEtConeDR03 + userFloat("hcalIsoConeDR03_2012")',
            'trkSumPtSolidConeDR04 + ecalRecHitSumEtConeDR04 + userFloat("hcalIsoConeDR04_2012")'
        )
    )
)
addphotonuserdata2 = cms.EDProducer("AddPhotonUserData",
    debug          = cms.bool(False),
    debugString    = cms.string("addphotonuserdata"),
    photonLabel    = cms.InputTag("patPhotonsUser1"),
    floatLabels    = cms.VInputTag(cms.InputTag("kt6PFJetsForIsolation","rho")),
    floatNames     = cms.vstring("rho25"),
    embedConversionInfo = cms.bool(False),
    gsfElectronLabel = cms.InputTag("gsfElectrons"),
    conversionsLabel = cms.InputTag("conversions"),
    beamspotLabel    = cms.InputTag("offlineBeamSpot"),
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
        userFunctionLabels = cms.vstring('pfChargedPUSub','pfNeutralPUSub','pfGammaPUSub',
                                         'pfChargedPURel','pfNeutralPURel','pfGammaPURel',
                                         'combIsoR03PU','combIsoR04PU'),
        userFunctions = cms.vstring(
            'userFloat("pfChargedEA")*userFloat("rho25")', #chargedSub
            'userFloat("pfNeutralEA")*userFloat("rho25")', #neutralSub
            'userFloat("pfGammaEA")  *userFloat("rho25")', #gammaSub
            'max((userIsolation("User1Iso") - userFloat("pfChargedEA")*userFloat("rho25"))/et,0.)',
            'max((userIsolation("User3Iso") - userFloat("pfNeutralEA")*userFloat("rho25"))/et,0.)',
            'max((userIsolation("User4Iso") - userFloat("pfGammaEA")  *userFloat("rho25"))/et,0.)',
            'trkSumPtSolidConeDR03 + ecalRecHitSumEtConeDR03 + userFloat("hcalIsoConeDR03_2012") - 3.141593*0.3*0.3*userFloat("rho25")',
            'trkSumPtSolidConeDR04 + ecalRecHitSumEtConeDR04 + userFloat("hcalIsoConeDR04_2012") - 3.141593*0.4*0.4*userFloat("rho25")'
        )
    )
)
