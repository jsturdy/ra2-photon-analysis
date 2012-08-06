import FWCore.ParameterSet.Config as cms

analysis = cms.EDAnalyzer('RA2ZInvPhotonTreeMaker',
                          Debug           = cms.bool(False),
                          Data            = cms.bool(True),
                          ScaleFactor     = cms.double(1.),
                          PhotonSrc       = cms.InputTag("patPhotonsUserData"),
                          JetSrc          = cms.InputTag("patJetsAK5PFPt30"),
                          bJetSrc         = cms.InputTag("patCSVJetsAK5PFPt30Eta24"),
                          JetHTSource     = cms.InputTag("patJetsAK5PFPt50Eta25"),
                          #                                  RA2NJets        = cms.uint32(3),
                          #                                  RA2HT           = cms.double(350.0),
                          #                                  RA2MHT          = cms.double(200.0),
                          #                                  RA2ApplyDphiCuts= cms.bool(True),
                          DoPUReweight    = cms.bool(False),
                          PUWeightSource  = cms.InputTag("puWeight")
)
