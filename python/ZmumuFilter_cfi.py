import FWCore.ParameterSet.Config as cms

ZmumuFilter = cms.EDFilter('ZmumuFilter',
                           muonInputTag = cms.InputTag("muons"),
                           csvFileName = cms.string("Zmumu.csv"),
                           minMuonPt = cms.double(25.0),
                           maxMuonEta = cms.double(2.1),
                           maxRelIso = cms.double(0.15),
                           invariantMassMin = cms.double(60.0),
                           invariantMassMax = cms.double(120.0)
                           )
