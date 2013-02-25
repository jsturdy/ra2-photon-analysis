
from PhysicsTools.PatAlgos.selectionLayer1.photonSelector_cfi import *
from PhysicsTools.PatAlgos.selectionLayer1.photonCountFilter_cfi import *

photonIDCutTight  = cms.string('et > 70.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && '
#                               'userInt("passElectronConvVeto") > 0 && '
                               'hadTowOverEm < userFloat("hadTowOverEmTightCut") && '
                               'sigmaIetaIeta < userFloat("showerShapeTightCut")'
                              )

photonIDCutMedium = cms.string('et > 70.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && '
#                               'userInt("passElectronConvVeto") > 0 && '
                               'hadTowOverEm < userFloat("hadTowOverEmMediumCut") && '
                               'sigmaIetaIeta < userFloat("showerShapeMediumCut")'
                              )

photonIDCutLoose  = cms.string('et > 70.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && '
#                               'userInt("passElectronConvVeto") > 0 && '
                               'hadTowOverEm < userFloat("hadTowOverEmLooseCut") && '
                               'sigmaIetaIeta < userFloat("showerShapeLooseCut")'
                              )

photoncombiso03cut = cms.string('trkSumPtSolidConeDR03 + ecalRecHitSumEtConeDR03 + '
                                'userFloat("hcalIsoConeDR03") - 3.141593*0.3*0.3*userFloat("rho25")< 5.'
                               )
photoncombiso04cut = cms.string('trkSumPtSolidConeDR04 + ecalRecHitSumEtConeDR04 + '
                                'userFloat("hcalIsoConeDR04") - 3.141593*0.4*0.4*userFloat("rho25")< 5.'
                               )

photonISOCutTight  = cms.string('userFloat("pfChargedPU") < userFloat("pfChargedTightCut") && '
                                'userFloat("pfNeutralPU") < userFloat("pfNeutralTightCut") && '
                                'userFloat("pfGammaPU")   < userFloat("pfGammaTightCut")'
                               )

photonISOCutMedium = cms.string('userFloat("pfChargedPU") < userFloat("pfChargedMediumCut") && '
                                'userFloat("pfNeutralPU") < userFloat("pfNeutralMediumCut") && '
                                'userFloat("pfGammaPU")   < userFloat("pfGammaMediumCut")'
                               )

photonISOCutLoose  = cms.string('userFloat("pfChargedPU") < userFloat("pfChargedLooseCut") && '
                                'userFloat("pfNeutralPU") < userFloat("pfNeutralLooseCut") && '
                                'userFloat("pfGammaPU")   < userFloat("pfGammaLooseCut")'
                               )

photonISOCutVeryLoose  = cms.string('userFloat("pfChargedPU") < userFloat("pfChargedVeryLooseCut") && '
                                    'userFloat("pfNeutralPU") < userFloat("pfNeutralVeryLooseCut") && '
                                    'userFloat("pfGammaPU")   < userFloat("pfGammaVeryLooseCut")'
                                    )

fitTemplatePhotonCut  = cms.string('et > 100.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && '
                                   'hadTowOverEm < userFloat("hadTowOverEmTightCut") && '
                                   'userFloat("pfChargedPU") < userFloat("pfChargedTightCut") &&'
                                   'userFloat("pfNeutralPU") < userFloat("pfNeutralTightCut") &&'
                                   'userFloat("pfGammaPU")   < userFloat("pfGammaTightCut")'
                                   )
bkgdTemplatePhotonCut  = cms.string('et > 100.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && '
                                    'hadTowOverEm < userFloat("hadTowOverEmTightCut") && '
                                    'userFloat("pfChargedPU") < min(userFloat("pfChargedVeryLooseCut"), 0.2*pt) && '
                                    'userFloat("pfNeutralPU") < min(userFloat("pfNeutralVeryLooseCut"), 0.2*pt) && '
                                    'userFloat("pfGammaPU")   < min(userFloat("pfGammaVeryLooseCut")  , 0.2*pt) && '
                                    
                                    '(userFloat("pfChargedPU") > userFloat("pfChargedTightCut") || '
                                    'userFloat("pfNeutralPU") > userFloat("pfNeutralTightCut") || '
                                    'userFloat("pfGammaPU")   > userFloat("pfGammaTightCut") )'
                               )
jetFakePhotonCut  = cms.string('et > 100.0 && (abs(eta) < 1.4442 || (abs(eta) > 1.566 && abs(eta) < 2.5)) && '
                               'hadTowOverEm < userFloat("hadTowOverEmVeryLooseCut") && '
                               'sigmaIetaIeta < userFloat("showerShapeVeryLooseCut") && '
                               'userFloat("pfChargedPU") < min(userFloat("pfChargedVeryLooseCut"), 0.2*pt) && '
                               'userFloat("pfNeutralPU") < min(userFloat("pfNeutralVeryLooseCut"), 0.2*pt) && '
                               'userFloat("pfGammaPU")   < min(userFloat("pfGammaVeryLooseCut")  , 0.2*pt) && '
                               
                               '( sigmaIetaIeta > userFloat("showerShapeTightCut") || '
                               'userFloat("pfChargedPU") > userFloat("pfChargedTightCut") || '
                               'userFloat("pfNeutralPU") > userFloat("pfNeutralTightCut") || '
                               'userFloat("pfGammaPU")   > userFloat("pfGammaTightCut") )'
                               )


patPhotonsID = cms.EDFilter(
  "PATPhotonSelector",
   src = cms.InputTag('patPhotonsUserData'),
   cut = photonIDCutTight,
   filter = cms.bool(False),
)
patPhotonsIDPFIsoLoose = patPhotonsID.clone(
    src = cms.InputTag('patPhotonsID'),
    cut = photonISOCutLoose,
)
patPhotonsIDPFIsoTight = patPhotonsID.clone(
    src = cms.InputTag('patPhotonsID'),
    cut = photonISOCutTight,
)

patPhotonsIDPFIso = patPhotonsID.clone(
    src = cms.InputTag('patPhotonsID'),
    cut = photonISOCutTight,
)

patJetFakePhotons = patPhotonsID.clone(
    src = cms.InputTag('patPhotonsUserData'),
    cut = jetFakePhotonCut,
)

patFitTemplatePhotons = patPhotonsID.clone(
    src = cms.InputTag('patPhotonsUserData'),
    cut = fitTemplatePhotonCut,
)

patBkgdTemplatePhotons = patPhotonsID.clone(
    src = cms.InputTag('patPhotonsUserData'),
    cut = bkgdTemplatePhotonCut,
)

patPhotonsIDCombIsoR03 = patPhotonsID.clone(
    src = cms.InputTag('patPhotonsID'),
    cut = photoncombiso04cut,
)

patPhotonsIDCombIsoR04 = patPhotonsID.clone(
   src = cms.InputTag('patPhotonsID'),
   cut = photoncombiso04cut,
)

patPhotonRefsID = cms.EDFilter(
  "PATPhotonRefSelector",
   src = cms.InputTag('patPhotonsUserData'),
   cut = photonIDCutLoose,
   filter = cms.bool(False)
)

zinvPhotons = cms.Sequence(
    patPhotonsID
  * patPhotonsIDPFIso
#  * patPhotonsIDPFIsoLoose
#  * patPhotonsIDPFIsoTight
  * patFitTemplatePhotons
  #* patBkgdTemplatePhotons
  * patJetFakePhotons
#  * patPhotonsIDCombIsoR03
#  * patPhotonsIDCombIsoR04
)

countPhotonsID = countPatPhotons.clone()
countPhotonsID.src = cms.InputTag("patPhotonsID")
countPhotonsID.minNumber = cms.uint32(1)

countMaxPhotonsID = countPhotonsID.clone()
countMaxPhotonsID.maxNumber = cms.uint32(1)

countPhotonsIDPFIso = countPhotonsID.clone()
countPhotonsIDPFIso.src = cms.InputTag("patPhotonsIDPFIso")

countMaxPhotonsIDPFIso = countPhotonsIDPFIso.clone()
countMaxPhotonsIDPFIso.maxNumber = cms.uint32(1)

countFitTemplatePhotons = countPhotonsID.clone()
countFitTemplatePhotons.src = cms.InputTag("patFitTemplatePhotons")

countMaxFitTemplatePhotons = countFitTemplatePhotons.clone()
countMaxFitTemplatePhotons.maxNumber = cms.uint32(1)

countBkgdTemplatePhotons = countPhotonsID.clone()
countBkgdTemplatePhotons.src = cms.InputTag("patBkgdTemplatePhotons")

countMaxBkgdTemplatePhotons = countBkgdTemplatePhotons.clone()
countMaxBkgdTemplatePhotons.maxNumber = cms.uint32(1)

countJetFakePhotons = countPhotonsID.clone()
countJetFakePhotons.src = cms.InputTag("patJetFakePhotons")

countMaxJetFakePhotons = countJetFakePhotons.clone()
countMaxJetFakePhotons.maxNumber = cms.uint32(1)

countPhotonRefsID = countPhotonsID.clone()
countPhotonRefsID.src = cms.InputTag("patPhotonsRefsID")
