import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer('RA2ZInvPhotonTemplateMaker',
                          Debug           = cms.bool(False),
                          Data            = cms.bool(True),
                          ScaleFactor     = cms.double(1.),
                          PhotonSrc       = cms.InputTag("patPhotonsUserData"),
                          JetSrc          = cms.InputTag("patJetsAK5PFPt30"),
                          bJetSrc         = cms.InputTag("patCSVJetsAK5PFPt30Eta24"),
                          JetHTSource     = cms.InputTag("patJetsAK5PFPt50Eta25"),
                          DoPUReweight    = cms.bool(False),
                          PUWeightSource  = cms.InputTag("puWeight")
)
