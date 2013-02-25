import FWCore.ParameterSet.Config as cms

photonTemplate = cms.EDAnalyzer('RA2ZInvPhotonTemplateMaker',
                          Debug           = cms.bool(False),
                          DebugString     = cms.string("photonTemplate"),
                          Data            = cms.bool(True),
                          ScaleFactor     = cms.double(1.),
                          PhotonSrc       = cms.InputTag("patPhotonsID"),
                          TightPhotonSrc  = cms.InputTag("patPhotonsIDPFIso"),
                          VertexSrc       = cms.InputTag("goodVertices"),
                          JetSrc          = cms.InputTag("patJetsPFNoPhotonIDSpecialPt30"),
                          htJetSrc        = cms.InputTag("patJetsPFNoPhotonIDSpecialPt50Eta25"),
                          htSource        = cms.InputTag("htPFchsNoPhotID"),
                          mhtSource       = cms.InputTag("mhtPFchsNoPhotID"),
                          metSource       = cms.InputTag("patMETsPF"),
                          DoPUReweight    = cms.bool(False),
                          PUWeightSource    = cms.InputTag("puWeight"   ,"weight"),
                          EventWeightSource = cms.InputTag("eventWeight","weight"),
)
