import FWCore.ParameterSet.Config as cms

hltobject = cms.EDAnalyzer("TriggerObjectAnalyzer",
                           processName = cms.string("HLT"),
                           treeName = cms.string("JetTriggers"),
                           triggerL1Name1 = cms.untracked.string("hltL1sL1DoubleEG2NotMinimumBiasHF2AND"),
                           triggerL1Name2 = cms.untracked.string("hltL1sL1EG5NotMinimumBiasHF2AND"),
                           triggerNames = cms.vstring(
                            "HLT_HIUPCL1SingleEG5NotHF2_v1",
                            "HLT_HIUPCL1DoubleEG2NotHF2_v1"
                            ),
                           triggerResults = cms.InputTag("TriggerResults","","HLT"),
                           triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")
)
