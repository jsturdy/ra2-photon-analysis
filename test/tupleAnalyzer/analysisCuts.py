import sys,os
import ROOT as r
from array import array
import math

def hadStopBaseline(event,nBjets,nJets):
    ra2_Vertices        = event.ra2_Vertices      
    ra2_nJetsPt30Eta24  = event.ra2_nJetsPt30Eta24 
    ra2_bJetsPt30Eta24  = event.ra2_bJetsPt30Eta24 
    ra2_MET             = event.ra2_MET            
    ra2_dPhi1           = event.ra2_dPhiMET1
    ra2_dPhi2           = event.ra2_dPhiMET2
    ra2_dPhi3           = event.ra2_dPhiMET3
    ra2_Jet1Pt          = event.ra2_Jet1Pt         
    ra2_Jet1Eta         = event.ra2_Jet1Eta        
    ra2_Jet2Pt          = event.ra2_Jet2Pt         
    ra2_Jet2Eta         = event.ra2_Jet2Eta        
    
    jet1Cuts = (ra2_Jet1Pt>70 and (abs(ra2_Jet1Eta)<2.4))
    jet2Cuts = (ra2_Jet2Pt>70 and (abs(ra2_Jet2Eta)<2.4))
    dPhiCuts = (ra2_dPhi1>0.5 and ra2_dPhi2>0.5 and ra2_dPhi3>0.3)
    bJetCuts = (ra2_bJetsPt30Eta24>(nBjets-1))
    basecuts = ((ra2_MET>175 and ra2_nJetsPt30Eta24>(nJets-1)) and jet1Cuts and jet2Cuts and dPhiCuts and bJetCuts)

    return basecuts

####
def topTaggerCuts(event):
    ##top tagger variables 
    ra2_bestTopJetMass  = event.ra2_bestTopJetMass 
    ra2_TbestTopJet     = event.ra2_TbestTopJet    
    ra2_TbJet           = event.ra2_TbJet          
    ra2_MT2             = event.ra2_MT2            
    
    #Cuts
    topTaggerCuts = (ra2_bestTopJetMass>80 and ra2_bestTopJetMass<270 and (ra2_TbJet+0.5*ra2_TbestTopJet)>500 and ra2_MT2>300)

    return topTaggerCuts

####

def ra2Baseline(event):
    ra2_Vertices        = event.ra2_Vertices      
    ra2_nJetsPt50Eta25  = event.ra2_nJetsPt50Eta25 
    ra2_HT              = event.ra2_HT            
    ra2_MHT             = event.ra2_MHT            
    ra2_dPhi1           = event.ra2_dPhiMHT1
    ra2_dPhi2           = event.ra2_dPhiMHT2
    ra2_dPhi3           = event.ra2_dPhiMHT3
    
    dPhiCuts = (ra2_dPhi1>0.5 and ra2_dPhi2>0.5 and ra2_dPhi3>0.3)
    basecuts = ((ra2_MHT>200 and ra2_HT>500 and ra2_nJetsPt50Eta25>1) and dPhiCuts)
    return basecuts
####
def ra2JetCuts(event,min,max):
    ra2_nJetsPt50Eta25  = event.ra2_nJetsPt50Eta25 
    
    return (ra2_nJetsPt50Eta25>=min and ra2_nJetsPt50Eta25<=max)
####
def ra2HTCuts(event,min,max):
    ra2_HT  = event.ra2_HT 
    
    return (ra2_HT>=min and ra2_HT<=max)
####
def ra2MHTCuts(event,min,max):
    ra2_MHT  = event.ra2_MHT 
    
    return (ra2_MHT>=min and ra2_MHT<=max)
####
def genAcceptance(event,minpt):
    ra2_genBosonPt  = event.ra2_genBoson1Pt
    ra2_genBosonEta = event.ra2_genBoson1Eta

    return (ra2_genBosonPt>minpt and abs(ra2_genBosonEta) < 2.5 and ((abs(ra2_genBosonEta) < 1.4442) or (abs(ra2_genBosonEta) > 1.566)))
####
def genRecoIsoEff(event,minpt,recoPt,passReco):
    ra2_nBosons  = event.ra2_nBosons
    ra2_bosonPt  = event.ra2_boson1Pt
    ra2_bosonEta = event.ra2_boson1Eta
    
    return (recoPt>minpt and passReco and abs(ra2_bosonEta) < 2.5 and ((abs(ra2_bosonEta) < 1.4442) or (abs(ra2_bosonEta) > 1.566)))

####
def genCuts(event,selection,minpt,recoPt,passReco,nJets,jetCuts,htVal,htCuts,mhtCuts):
    if nJets>=jetCuts[0] and nJets<=jetCuts[1]:
        if selection==0:
            return event.ra2_genBoson1Pt > minpt
        elif selection==1:
            return genAcceptance(event,minpt)
        elif selection==2:
            return genAcceptance(event,minpt) and genRecoIsoEff(event,minpt,recoPt,passReco)
        elif selection==3:
            return genAcceptance(event,minpt) and genRecoIsoEff(event,minpt,recoPt,passReco)
        elif selection==4:
            return genAcceptance(event,minpt) and genRecoIsoEff(event,minpt,recoPt,passReco) and htVal>350 and event.ra2_genMHT>150 and passGenDPhiCuts(event)
        else:
            return False
    else:
        return False

####
def genCutsFull(event,selection,minpt,recoPt,passReco,nJets,jetCuts,htVal,htCuts,mhtCuts):
    passDPhi = passGenDPhiCuts(event)

    passNJets  = False
    if nJets>=jetCuts[0] and nJets<=jetCuts[1]:
        passNJets  = True
    ###
    passHT  = False
    if htVal > htCuts[0] and htVal < htCuts[1]:
        passHT  = True
    ###
    passMHT = False
    if event.ra2_genMHT > mhtCuts[0] and event.ra2_genMHT < mhtCuts[1]:
        passMHT  = True

    ######    
    if passNJets and passHT and passMHT and passDPhi:
        if selection==0:
            return event.ra2_genBoson1Pt > 0
        elif selection==1:
            return event.ra2_genBoson1Pt > minpt
        elif selection==2:
            return genAcceptance(event,minpt)
        elif selection==3:
            return genAcceptance(event,minpt) and genRecoIsoEff(event,minpt,recoPt,passReco)
        elif selection==4:
            return genAcceptance(event,minpt) and genRecoIsoEff(event,minpt,recoPt,passReco)
        else:
            return False
    else:
        return False

####
def phenoCuts(event,cutLevel,minpt,nJets,jetCuts,htVal,htCuts,mhtCuts):
    passDPhi = True
    passBase = True
    if cutLevel > 0:
        passDPhi = passGenDPhiCuts(event)
    if cutLevel > 1:
        passBase = htVal > 500 and event.ra2_genMHT > 200

    passNJets  = False
    if nJets>=jetCuts[0] and nJets<=jetCuts[1]:
        passNJets  = True
    ###
    passHT  = False
    if htVal > htCuts[0] and htVal < htCuts[1]:
        passHT  = True
    ###
    passMHT = False
    if event.ra2_genMHT > mhtCuts[0] and event.ra2_genMHT < mhtCuts[1]:
        passMHT  = True
    #print "passDphi %d(%d,%d), nJets %d, bosonPt %2.2f"%(passGenDPhiCuts(event),dphi,passDPhi, nJets, event.ra2_genBoson1Pt)
    if passNJets and passHT and passMHT and passDPhi and passBase:
        return event.ra2_genBoson1Pt > minpt
    else:
        return False

def passGenDPhiCuts(event):
    return event.ra2_genDPhiMHT1 > 0.5 and event.ra2_genDPhiMHT2 > 0.5 and event.ra2_genDPhiMHT3 > 0.3
