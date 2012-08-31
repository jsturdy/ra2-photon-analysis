import FWCore.ParameterSet.Config as cms

photonTree = cms.EDAnalyzer('RA2ZInvPhotonTreeMaker',
                          Debug           = cms.bool(False),
                          Data            = cms.bool(True),
                          ScaleFactor     = cms.double(1.),
                          PhotonSrc       = cms.InputTag("patPhotonsUserData"),
                          VertexSrc       = cms.InputTag("goodVertices"),
                          JetSrc          = cms.InputTag("patJetsPFNoPhotonPt30"),
                          htJetSrc        = cms.InputTag("patJetsPFNoPhotonPt50Eta25"),
                          bJetSrc         = cms.InputTag("patCSVJetsPFNoPhotonPt30Eta24"),
                          htSource        = cms.InputTag("htPFchsNoPhot"),
                          mhtSource       = cms.InputTag("mhtPFchsNoPhot"),
                          DoPUReweight    = cms.bool(True),
                          PUWeightSource  = cms.InputTag("puWeight")
)

dimuonTree = cms.EDAnalyzer('RA2ZInvDiMuonTreeMaker',
                          Debug           = cms.bool(False),
                          Data            = cms.bool(True),
                          ScaleFactor     = cms.double(1.),
                          MuonSrc         = cms.InputTag("specialMuonCollection"),
                          VertexSrc       = cms.InputTag("goodVertices"),
                          JetSrc          = cms.InputTag("patJetsPFNoMuonPt30"),
                          htJetSrc        = cms.InputTag("patJetsPFNoMuonPt50Eta25"),
                          bJetSrc         = cms.InputTag("patCSVJetsPFNoMuonPt30Eta24"),
                          htSource        = cms.InputTag("htPFchsNoMuon"),
                          mhtSource       = cms.InputTag("mhtPFchsNoMuon"),
                          DoPUReweight    = cms.bool(True),
                          PUWeightSource  = cms.InputTag("puWeight")
)

zvvTree = cms.EDAnalyzer('RA2ZInvTreeMaker',
                          Debug           = cms.bool(False),
                          ScaleFactor     = cms.double(1.),
                          genLabel        = cms.InputTag("zinvBkgdst3ZBosons"),
                          VertexSrc       = cms.InputTag("goodVertices"),
                          JetSrc          = cms.InputTag("patJetsPFchsPt30"),
                          htJetSrc        = cms.InputTag("patJetsPFchsPt50Eta25"),
                          bJetSrc         = cms.InputTag("patCSVJetsPFPt30Eta24"),
                          htSource        = cms.InputTag("htPFchs"),
                          mhtSource       = cms.InputTag("mhtPFchs"),
                          DoPUReweight    = cms.bool(True),
                          PUWeightSource  = cms.InputTag("puWeight")
)
