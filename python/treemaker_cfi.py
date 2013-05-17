import FWCore.ParameterSet.Config as cms

photonTree = cms.EDAnalyzer('RA2ZInvPhotonTreeMaker',
                          Debug           = cms.bool(False),
                          DebugString     = cms.string("photonID"),
                          Data            = cms.bool(True),
                          ScaleFactor     = cms.double(1.),

                          PhotonSrc       = cms.InputTag("patPhotonsID"),
                          TightPhotonSrc  = cms.InputTag("patPhotonsIDPFIso"),

                          VertexSrc       = cms.InputTag("goodVertices"),

                          JetSrc          = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30"),
                          htJetSrc        = cms.InputTag("patJetsPFNoPhotonIDSpecialPt50Eta25"),
                          bJetSrc         = cms.InputTag("patCSVMJetsPFNoPhotonIDSpecialPt30Eta24"),

                          htSource        = cms.InputTag("htPFchsNoPhotID"),
                          mhtSource       = cms.InputTag("mhtPFchsNoPhotID"),
                          ra2HTSource     = cms.InputTag("htPFchs"),
                          ra2MHTSource    = cms.InputTag("mhtPFchs"),

                          computeMET = cms.bool(True),
                          metSource       = cms.InputTag("pfType1MetNoPhotID","pfcand"),
                          ra2METSource    = cms.InputTag("newMETwPhiCorr"),

                          TriggerResults  = cms.InputTag("hltTriggerSummaryAOD"),
                          DoPUReweight    = cms.bool(False),
                          PUWeightSource    = cms.InputTag("puWeight"   ,"weight"),
                          EventWeightSource = cms.InputTag("eventWeight","weight"),

                          ##runTopTagger    = cms.bool(False),
                          ##topTaggerSource = cms.string("topTaggerType3"),

                          ##objects for the gen study
                          runGenStudy = cms.bool(False),
                          genSrc      = cms.InputTag("genParticles"),
                          maxDR       = cms.double(0.2),
                          genJetSrc   = cms.InputTag("ak5GenJets"),
                          ra2JetSrc   = cms.InputTag("patJetsPFchsPt30"),
                          genMETSrc   = cms.InputTag("genMetTrue"),

                          storeExtraVetos    = cms.bool(True),
                          ra2ElectronForVeto = cms.InputTag("ra2ElectronVeto"),
                          ra2MuonForVeto     = cms.InputTag("ra2PFMuonVeto"),
                          electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
                          muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
                          isoTrkVetoSource   = cms.InputTag("sTopTrkIsolationMaker","trkIsoVeto"),
)

dimuonTree = cms.EDAnalyzer('RA2ZInvDiMuonTreeMaker',
                          Debug           = cms.bool(False),
                          Data            = cms.bool(True),
                          ScaleFactor     = cms.double(1.),
                          MuonSrc         = cms.InputTag("specialMuonCollection"),
                          VertexSrc       = cms.InputTag("goodVertices"),

                          JetSrc          = cms.InputTag("patJetsPFNoMuonPt30"),
                          htJetSrc        = cms.InputTag("patJetsPFNoMuonPt50Eta25"),
                          bJetSrc         = cms.InputTag("patCSVMJetsPFNoMuonPt30Eta24"),

                          htSource        = cms.InputTag("htPFchsNoMuon"),
                          mhtSource       = cms.InputTag("mhtPFchsNoMuon"),
                          ra2HTSource     = cms.InputTag("htPFchs"),
                          ra2MHTSource    = cms.InputTag("mhtPFchs"),

                          computeMET = cms.bool(True),
                          metSource       = cms.InputTag("pfType1MetNoMuon","pfcand"),
                          ra2METSource    = cms.InputTag("newMETwPhiCorr"),

                          TriggerResults  = cms.InputTag("hltTriggerSummaryAOD"),
                          DoPUReweight    = cms.bool(False),
                          PUWeightSource    = cms.InputTag("puWeight"   ,"weight"),
                          EventWeightSource = cms.InputTag("eventWeight","weight"),

                          runTopTagger    = cms.bool(False),
                          topTaggerSource = cms.string("topTaggerType3"),

                          ##objects for the gen study
                          runGenStudy = cms.bool(False),
                          genSrc      = cms.InputTag("zinvBkgdst3ZMuMuBosons"),
                          genMuonSrc  = cms.InputTag("zinvBkgdGenMuons"),
                          maxDR       = cms.double(0.2),
                          genJetSrc   = cms.InputTag("ak5GenJets"),
                          ra2JetSrc   = cms.InputTag("patJetsPFchsPt30"),
                          genMETSrc   = cms.InputTag("genMetTrue"),

                          storeExtraVetos    = cms.bool(True),
                          ra2ElectronForVeto = cms.InputTag("ra2ElectronVeto"),
                          ra2MuonForVeto     = cms.InputTag("ra2PFMuonVetoDiMuon"),
                          electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
                          muonVetoSource     = cms.InputTag("sTopPFMuonVetoDiMuon"),
                          isoTrkVetoSource   = cms.InputTag("sTopTrkIsolationMakerDiLeptons","trkIsoVeto"),
)

dielectronTree = cms.EDAnalyzer('RA2ZInvDiElectronTreeMaker',
                          Debug           = cms.bool(False),
                          Data            = cms.bool(True),
                          ScaleFactor     = cms.double(1.),
                          ElectronSrc     = cms.InputTag("specialElectronCollection"),
                          VertexSrc       = cms.InputTag("goodVertices"),
                          JetSrc          = cms.InputTag("patJetsPFNoElectronPt30"),
                          htJetSrc        = cms.InputTag("patJetsPFNoElectronPt50Eta25"),
                          bJetSrc         = cms.InputTag("patCSVMJetsPFNoElectronPt30Eta24"),
                          htSource        = cms.InputTag("htPFchsNoElectron"),
                          mhtSource       = cms.InputTag("mhtPFchsNoElectron"),
                          metSource       = cms.InputTag("pfType1MetNoElectron","pfcand"),
                          ra2HTSource     = cms.InputTag("htPFchs"),
                          ra2MHTSource    = cms.InputTag("mhtPFchs"),
                          ra2METSource    = cms.InputTag("newMETwPhiCorr"),
                          TriggerResults  = cms.InputTag("hltTriggerSummaryAOD"),

                          runTopTagger    = cms.bool(False),
                          topTaggerSource = cms.string("topTaggerType3"),
                          DoPUReweight    = cms.bool(False),
                          PUWeightSource    = cms.InputTag("puWeight"   ,"weight"),
                          EventWeightSource = cms.InputTag("eventWeight","weight"),

                          storeExtraVetos    = cms.bool(True),
                          ra2ElectronForVeto = cms.InputTag("ra2ElectronVetoDiElectron"),
                          ra2MuonForVeto     = cms.InputTag("ra2PFMuonVeto"),
                          electronVetoSource = cms.InputTag("sTopPFElectronVetoDiElectron"),
                          muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
                          isoTrkVetoSource   = cms.InputTag("sTopTrkIsolationMakerDiLeptons","trkIsoVeto"),
)

zvvTree = cms.EDAnalyzer('RA2ZInvTreeMaker',
                         Debug           = cms.bool(False),
                         ScaleFactor     = cms.double(1.),
                         
                         genLabel        = cms.InputTag("zinvBkgdst3ZBosons"),
                         #genSrc          = cms.InputTag("zinvBkgdst3ZBosons"),
                         genJetSrc       = cms.InputTag("ak5GenJets"),
                         genMETSrc       = cms.InputTag("genMetTrue"),
                         
                         VertexSrc       = cms.InputTag("goodVertices"),
                         
                         JetSrc          = cms.InputTag("patJetsPFchsPt30"),
                         htJetSrc        = cms.InputTag("patJetsPFchsPt50Eta25"),
                         bJetSrc         = cms.InputTag("patCSVMJetsPFPt30Eta24"),
                         
                         htSource        = cms.InputTag("htPFchs"),
                         mhtSource       = cms.InputTag("mhtPFchs"),
                         metSource       = cms.InputTag("newMETwPhiCorr"),
                         
                         DoPUReweight    = cms.bool(True),
                         PUWeightSource    = cms.InputTag("puWeight"   ,"weight"),
                         EventWeightSource = cms.InputTag("eventWeight","weight"),
                         
                         runTopTagger    = cms.bool(False),
                         topTaggerSource = cms.string("topTaggerType3"),
                         
                         storeExtraVetos    = cms.bool(True),
                         ra2ElectronForVeto = cms.InputTag("ra2ElectronVeto"),
                         ra2MuonForVeto     = cms.InputTag("ra2PFMuonVeto"),
                         electronVetoSource = cms.InputTag("sTopPFElectronVeto"),
                         muonVetoSource     = cms.InputTag("sTopPFMuonVeto"),
                         isoTrkVetoSource   = cms.InputTag("sTopTrkIsolationMaker","trkIsoVeto"),
)
