#
#  SUSY-PAT configuration file adapted for RA2 workflow
#
#  PAT configuration for the SUSY group - 52X series
#  More information here:
#  https://twiki.cern.ch/twiki/bin/view/CMS/SusyPatLayer1DefV12
#

# Starting with a skeleton process which gets imported with the following line
from PhysicsTools.PatAlgos.patTemplate_cfg import *

#runningOnMC = True 
runningOnMC = False

#-- Message Logger ------------------------------------------------------------
process.MessageLogger.categories.append('PATSummaryTables')
process.MessageLogger.cerr.PATSummaryTables = cms.untracked.PSet(
    limit = cms.untracked.int32(10),
    reportEvery = cms.untracked.int32(1)
    )
process.MessageLogger.cerr.FwkReport.reportEvery = 1


#-- Input Source --------------------------------------------------------------
process.maxEvents.input = 100

# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

############################# START SUSYPAT specifics ####################################
from PhysicsTools.Configuration.SUSY_pattuple_cff import addDefaultSUSYPAT, getSUSY_pattuple_outputCommands

hltMenu = 'HLT'

theJetColls = ['AK5PF']

jetMetCorr = ['L1FastJet', 'L2Relative', 'L3Absolute']
if runningOnMC == False: jetMetCorr.append('L2L3Residual')  

process.GlobalTag.globaltag = "START52_V5::All"
if runningOnMC == False:
    process.GlobalTag.globaltag = "GR_R_52_V7::All"

process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring('')
)
if runningOnMC == False:
    process.source = cms.Source("PoolSource",
     fileNames = cms.untracked.vstring('/store/data/Run2012B/HTMHT/AOD/PromptReco-v1/000/194/052/9A12DE52-699E-E111-BC9F-5404A638869E.root')
    )

# Due to problem in production of LM samples: same event number appears multiple times
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

addDefaultSUSYPAT(process,mcInfo=runningOnMC,HLTMenu=hltMenu,jetMetCorrections=jetMetCorr,mcVersion='',theJetNames=theJetColls, doSusyTopProjection=True)

#from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso, setupPFPhotonIso
from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFPhotonIso
#process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
#process.muIsoSequence = setupPFMuonIso(process, 'muons')
process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

process.patPhotons.isoDeposits = cms.PSet(
    user = cms.VInputTag(
        cms.InputTag("phPFIsoDepositChargedPFIso"),
        cms.InputTag("phPFIsoDepositChargedAllPFIso"),
        cms.InputTag("phPFIsoDepositNeutralPFIso"),
        cms.InputTag("phPFIsoDepositGammaPFIso"),
        cms.InputTag("phPFIsoDepositPUPFIso")
    ),
)
process.patPhotons.userIsolation = cms.PSet(
    user = cms.VPSet(
        cms.PSet( src = cms.InputTag("phPFIsoValueCharged03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValueChargedAll03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValueNeutral03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValueGamma03PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValuePU03PFIdPFIso"))
    )
)
process.patPhotons.userData.userFloats = cms.PSet(
    user = cms.VPSet(
        cms.PSet( src = cms.InputTag("kt6PFJetsForIsolation","rho")),
    )
)
process.patPhotonsR04 = getattr(process,'patPhotons').clone()
process.patPhotonsR04.userIsolation = cms.PSet(
    user = cms.VPSet(
        cms.PSet( src = cms.InputTag("phPFIsoValueCharged04PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValueChargedAll04PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValueNeutral04PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValueGamma04PFIdPFIso")),
        cms.PSet( src = cms.InputTag("phPFIsoValuePU04PFIdPFIso"))
    )
)


#SUSY_pattuple_outputCommands = getSUSY_pattuple_outputCommands( process )

# overwrite default output content
from SandBox.Skims.RA2Content_cff import getRA2PATOutput
process.out.outputCommands = getRA2PATOutput(process)
process.out.dropMetaData = cms.untracked.string('DROPPED')
############################## END SUSYPAT specifics ####################################

#-- HLT selection ------------------------------------------------------------
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
hltSelection=''
#if options.hltSelection:
if hltSelection:
    process.hltFilter = hlt.hltHighLevel.clone(
        HLTPaths = cms.vstring(hltSelection),
        TriggerResultsTag = cms.InputTag("TriggerResults","",hltMenu),
        throw = False
    )
    process.susyPatDefaultSequence.replace(process.eventCountProducer, process.eventCountProducer * process.hltFilter)


#-- check RA2 recipe here ------------------------------------------------------------
process.prefilterCounter        = cms.EDProducer("EventCountProducer")
process.postStdCleaningCounter  = cms.EDProducer("EventCountProducer")

#-- Output module configuration -----------------------------------------------
process.load('SandBox.Skims.RA2Objects_cff')
process.load('SandBox.Skims.RA2Selection_cff')
process.load('SandBox.Skims.RA2Cleaning_cff')
process.load('SandBox.Skims.allMETs_cff')

process.cleanpatseq = cms.Sequence(
                      process.susyPatDefaultSequence  *
                      process.kt6PFJetsForIsolation   *
                      process.phoIsoSequence          *
                      process.patPhotonsR03           *
                      process.patPhotonsR04           *
                      process.prefilterCounter        *
                      process.ra2StdCleaning          *
                      process.postStdCleaningCounter  *
                      process.ra2Objects              *
                      process.produceAllPFMETWithCHSCorrections *
                      process.produceAllPFMETWithC2VCorrections *
                      process.produceAllPatPFMETCorrections     
                      #process.ra2PostCleaning         *
                      #process.ra2FullPFSelectionNoMHT *
                      #process.mhtPFFilter
                      )

process.ppf = cms.Path(
              process.cleanpatseq
              )

#-- Output module configuration ---------------------------------------
process.out.fileName = cms.untracked.string('susypat_default.root')

process.out.SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('ppf') )
process.outpath = cms.EndPath( process.out )

###-- Dump config ------------------------------------------------------------
file = open('SusyPAT_RA2_data525p1_new_cfg.py','w')
file.write(str(process.dumpPython()))
file.close()
