import FWCore.ParameterSet.Config as cms

selectedRA2PatJets = cms.EDFilter("RA2BasicJetSelector",
                                  src = cms.InputTag("newJetsMET"),
                                  cut = cms.string("")
                                  )
