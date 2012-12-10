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
        userFunctionLabels = cms.vstring('maxDz', 'maxRelIso','chargedIso','otherIso',
                                         'pfIsoEA', 'pfNeutralEA', 'pfGammaEA',
                                         'dEtaInCut', 'dPhiInCut', 'dzCut',
                                         'hadTowOverEmCut', 'showerShapeCut', 'pfIsoRelCut',
                                         'pfIsoRel', 'pfNeutralRel', 'pfGammaRel'
                                         ),
        userFunctions = cms.vstring(
            """?((0.0   <= pt < 5.0  ) and (0 <= abs(superCluster.eta) < 1.5)) ? 0.03 :
               ?((5.0   <= pt < 10.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 0.05 :
               ?((10.0  <= pt < 15.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 0.05 :
               ?((15.0  <= pt < 20.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 0.05 :
               ?((20.0  <= pt < 40.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 0.2 :
               ?((40.0  <= pt < 80.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 1.0 :
               ?((80.0  <= pt < 140.0) and (0 <= abs(superCluster.eta) < 1.5)) ? 1.0 :
               ?((140.0 <= pt < 200.0) and (0 <= abs(superCluster.eta) < 1.5)) ? 1.0 :
               ?((200.0 <= pt)         and (0 <= abs(superCluster.eta) < 1.5)) ? 1.0 :
               ?((0.0   <= pt < 5.0  ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 0.09 :
               ?((5.0   <= pt < 10.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 0.09 :
               ?((10.0  <= pt < 15.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 0.09 :
               ?((15.0  <= pt < 20.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 0.11 :
               ?((20.0  <= pt < 40.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 0.1 :
               ?((40.0  <= pt < 80.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 1.0 :
               ?((80.0  <= pt < 140.0) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 1.0 :
               ?((140.0 <= pt < 200.0) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 1.0 :
               ?((200.0 <= pt)         and (1.5 <= abs(superCluster.eta) < 2.5)) ? 1.0 :
                                                                                   1.0""", #maxDz
            """?((0.0   <= pt < 5.0  ) and (0 <= abs(superCluster.eta) < 1.5)) ? 0.5 :
               ?((5.0   <= pt < 10.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 1.5 :
               ?((10.0  <= pt < 15.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 4.5 :
               ?((15.0  <= pt < 20.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 7.5 :
               ?((20.0  <= pt < 40.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 10.0:
               ?((40.0  <= pt < 80.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 18.5:
               ?((80.0  <= pt < 140.0) and (0 <= abs(superCluster.eta) < 1.5)) ? 44.0:
               ?((140.0 <= pt < 200.0) and (0 <= abs(superCluster.eta) < 1.5)) ? 81.5:
               ?((200.0 <= pt)         and (0 <= abs(superCluster.eta) < 1.5)) ? 81.5:
               ?((0.0   <= pt < 5.0  ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 0.5 :
               ?((5.0   <= pt < 10.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 2.5 :
               ?((10.0  <= pt < 15.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 6.5 :
               ?((15.0  <= pt < 20.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 9.0 :
               ?((20.0  <= pt < 40.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 10.5:
               ?((40.0  <= pt < 80.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 18.5:
               ?((80.0  <= pt < 140.0) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 66.5:
               ?((140.0 <= pt < 200.0) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 70.0:
               ?((200.0 <= pt)         and (1.5 <= abs(superCluster.eta) < 2.5)) ? 70.0:
                                                                                   70.0""", #chargedIso
            """?((0.0   <= pt < 5.0  ) and (0 <= abs(superCluster.eta) < 1.5)) ? 0.5 :
               ?((5.0   <= pt < 10.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 1.5 :
               ?((10.0  <= pt < 15.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 4.5 :
               ?((15.0  <= pt < 20.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 7.5 :
               ?((20.0  <= pt < 40.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 10.0:
               ?((40.0  <= pt < 80.0 ) and (0 <= abs(superCluster.eta) < 1.5)) ? 18.5:
               ?((80.0  <= pt < 140.0) and (0 <= abs(superCluster.eta) < 1.5)) ? 44.0:
               ?((140.0 <= pt < 200.0) and (0 <= abs(superCluster.eta) < 1.5)) ? 81.5:
               ?((200.0 <= pt)         and (0 <= abs(superCluster.eta) < 1.5)) ? 81.5:
               ?((0.0   <= pt < 5.0  ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 0.5 :
               ?((5.0   <= pt < 10.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 2.5 :
               ?((10.0  <= pt < 15.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 6.5 :
               ?((15.0  <= pt < 20.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 9.0 :
               ?((20.0  <= pt < 40.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 10.5:
               ?((40.0  <= pt < 80.0 ) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 18.5:
               ?((80.0  <= pt < 140.0) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 66.5:
               ?((140.0 <= pt < 200.0) and (1.5 <= abs(superCluster.eta) < 2.5)) ? 70.0:
               ?((200.0 <= pt)         and (1.5 <= abs(superCluster.eta) < 2.5)) ? 70.0:
                                                                                   70.0""", #relIso
            '(pfIsolationR03().sumChargedHadronPt+max(0.,pfIsolationR03().sumNeutralHadronEt+pfIsolationR03().sumPhotonEt-0.5*pfIsolationR03().sumPUPt))/pt()'

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
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.007:0.01""",#dEtaIn
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.8 :0.7""", #dPhiIn
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.4 :0.4""", #d0
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.2  :0.2""",  #dZ
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.15 :1e6""", #singleTowerHOverE
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.01 :0.03""", #showerShape
            """?0.0   <= abs(superCluster.eta) < 1.4442? 0.15 :0.15""", #pfIsoRel

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
