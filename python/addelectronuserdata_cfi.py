import FWCore.ParameterSet.Config as cms

addelectronuserdata1 = cms.EDProducer("AddElectronUserData",
    debug          = cms.bool(False),
    debugString    = cms.string("addelectronuserdata"),
    electronLabel    = cms.InputTag("patElectrons"),
    floatLabels    = cms.VInputTag(cms.InputTag("kt6PFJetsForIsolation","rho")),
    floatNames     = cms.vstring("rho25"),
    embedConversionInfo = cms.bool(True),
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
        userFunctionLabels = cms.vstring('pfIsoEA', 'pfNeutralEA', 'pfGammaEA',

                                         'dEtaInTightCut', 'dPhiInTightCut', 'd0TightCut', 'dzTightCut',
                                         'invEminPTightCut', 'vtxFitProbTightCut', 'misHitsTightCut',
                                         'hadTowOverEmTightCut', 'showerShapeTightCut', 'pfIsoRelTightCut',

                                         'dEtaInMediumCut', 'dPhiInMediumCut', 'd0MediumCut', 'dzMediumCut',
                                         'invEminPMediumCut', 'vtxFitProbMediumCut', 'misHitsMediumCut',
                                         'hadTowOverEmMediumCut', 'showerShapeMediumCut', 'pfRelIsoMediumCut',

                                         'dEtaInLooseCut', 'dPhiInLooseCut', 'd0LooseCut', 'dzLooseCut',
                                         'invEminPLooseCut', 'vtxFitProbLooseCut', 'misHitsLooseCut',
                                         'hadTowOverEmLooseCut', 'showerShapeLooseCut', 'pfIsoRelLooseCut',

                                         'pfIsoRel', 'pfNeutralRel', 'pfGammaRel'),
        userFunctions = cms.vstring(
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.135 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.168 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.068 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.116 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.160 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.241 :
                                                         0.240""", #neutral+gammaEA
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.013 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.021 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.013 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.010 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.024 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.020 :
                                                         0.019""", #neutralEA
            """?0.0   <= abs(superCluster.eta) < 1.0   ? 0.122 :
               ?1.0   <= abs(superCluster.eta) < 1.479 ? 0.147 :
               ?1.479 <= abs(superCluster.eta) < 2.0   ? 0.055 :
               ?2.0   <= abs(superCluster.eta) < 2.2   ? 0.106 :
               ?2.2   <= abs(superCluster.eta) < 2.3   ? 0.136 :
               ?2.3   <= abs(superCluster.eta) < 2.4   ? 0.221 :
                                                         0.211""", #gammaEA

            ####Cut values:singleTowerHOverE,showerShape,pfChargedRel,pfNeutralRel,pfGammaRel
            ###Tight
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.004:0.005""",#dEtaIn
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.03 :0.02""", #dPhiIn
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.02 :0.02""", #d0
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.1  :0.1""",  #dZ
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05""", #fabs(1/E-1/p)
            """?0.0   <= abs(superCluster.eta) < 1.4442? 1e-6 :1e-6""", #vtxFitProb
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0    :0""",    #missHits
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.12 :0.10""", #singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.01 :0.03""", #showerShape
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.10 :0.10""", #pfIsoRel
            ###Medium
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.004:0.007""",#dEtaIn
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.06 :0.03""", #dPhiIn
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.02 :0.02""", #d0
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.1  :0.1""",  #dZ
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05""", #fabs(1/E-1/p)
            """?0.0   <= abs(superCluster.eta) < 1.4442? 1e-6 :1e-6""", #vtxFitProb
            """?0.0   <= abs(superCluster.eta) < 1.4442? 1    :1""",    #missHits
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.12:0.10""",  #singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.01:0.03""",  #showerShape      
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.15:0.15""",  #pfIsoRel     
            ###Loose
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.007:0.009""",#dEtaIn
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.15 :0.10""", #dPhiIn
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.02 :0.02""", #d0
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.2  :0.2""",  #dZ
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.05 :0.05""", #fabs(1/E-1/p)
            """?0.0   <= abs(superCluster.eta) < 1.4442? 1e-6 :1e-6""", #vtxFitProb
            """?0.0   <= abs(superCluster.eta) < 1.4442? 1    :1""",    #missHits
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.12:0.20""",  #singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.01:0.03""",  #showerShape      
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.15:0.15""",  #pfIsoRel     

            '(neutralHadronIso+photonIso)/et',
            'neutralHadronIso/et',
            'photonIso/et'
            )
    )
)
addelectronuserdata2 = cms.EDProducer("AddElectronUserData",
    debug          = cms.bool(False),
    debugString    = cms.string("addelectronuserdata"),
    electronLabel    = cms.InputTag("patElectronsUser1"),
    floatLabels    = cms.VInputTag(cms.InputTag("kt6PFJetsForIsolation","rho")),
    floatNames     = cms.vstring("rho25"),
    embedConversionInfo = cms.bool(False),
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
        userFunctionLabels = cms.vstring('pfIsoPUSub','pfNeutralPUSub','pfGammaPUSub',
                                         'pfIsoPURel','pfNeutralPURel','pfGammaPURel'),
        userFunctions = cms.vstring(
            'userFloat("pfIsoEA")    *userFloat("rho25")', #netural+gamma Sub
            'userFloat("pfNeutralEA")*userFloat("rho25")', #neutralSub
            'userFloat("pfGammaEA")  *userFloat("rho25")', #gammaSub

            'max((neutralHadronIso+photonIso - userFloat("pfIsoEA")    *userFloat("rho25"))/et,0.)',
            'max((neutralHadronIso           - userFloat("pfNeutralEA")*userFloat("rho25"))/et,0.)',
            'max((photonIso                  - userFloat("pfGammaEA")  *userFloat("rho25"))/et,0.)'
        )
    )
)
