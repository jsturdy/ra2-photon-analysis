import sys,os
import ROOT as r
from array import array
import math

def hadStopBaseline(event,nBjets,tightCSV,nJets,invert):
    ra2_Vertices        = event.ra2_Vertices      
    ra2_nJetsPt30Eta24  = event.ra2_nJetsPt30Eta24 
    ra2_nJetsCSVM       = event.ra2_nJetsCSVM 
    ra2_nJetsCSVT       = event.ra2_nJetsCSVT 

    ra2_MET             = event.ra2_MET            
    ra2_dPhi1           = event.ra2_dPhiMET1
    ra2_dPhi2           = event.ra2_dPhiMET2
    ra2_dPhi3           = event.ra2_dPhiMET3

    ra2_Jet1Pt          = event.ra2_Jet1Pt         
    ra2_Jet1Eta         = event.ra2_Jet1Eta        
    ra2_Jet2Pt          = event.ra2_Jet2Pt         
    ra2_Jet2Eta         = event.ra2_Jet2Eta        
    ra2_Jet3Pt          = event.ra2_Jet3Pt         
    ra2_Jet3Eta         = event.ra2_Jet3Eta        
    ra2_Jet4Pt          = event.ra2_Jet4Pt         
    ra2_Jet4Eta         = event.ra2_Jet4Eta        
    
    jet1Cuts = (ra2_Jet1Pt>70 and (abs(ra2_Jet1Eta)<2.4))
    jet2Cuts = (ra2_Jet2Pt>70 and (abs(ra2_Jet2Eta)<2.4))
    jet3Cuts = (ra2_Jet3Pt>50 and (abs(ra2_Jet3Eta)<2.4))
    jet4Cuts = (ra2_Jet4Pt>50 and (abs(ra2_Jet4Eta)<2.4))
    jetCuts = jet1Cuts and jet2Cuts and jet3Cuts and jet4Cuts and ra2_nJetsPt30Eta24 > (nJets-1)

    dPhiCuts = (ra2_dPhi1>0.5 and ra2_dPhi2>0.5 and ra2_dPhi3>0.3)

    bJetCuts = (ra2_nJetsCSVM>(nBjets-1))
    if tightCSV:
        bJetCuts = (ra2_nJetsCSVT>(nBjets-1))
    invertcuts = (ra2_MET<175 and jetCuts and dPhiCuts and bJetCuts)
    basecuts   = (ra2_MET>175 and jetCuts and dPhiCuts and bJetCuts)

    if invert:
        return invertcuts
    else:
        return basecuts
####
def hadStopControlSample(event,cutMET):
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
    jetCuts = jet1Cuts and jet2Cuts

    dPhiCuts = (ra2_dPhi1>0.5 and ra2_dPhi2>0.5 and ra2_dPhi3>0.3)
    cscuts   = (jetCuts and dPhiCuts)
    if cutMET:
        return cscuts and (ra2_MET > 175)
    else:
        return cscuts
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
def ra2PhotonSelection(event,minpt,cutDR):
    ra2_nPhotons  = event.ra2_nPhotonsIso
    ra2_photonPt  = event.ra2_Photon1Pt
    ra2_photonEta = event.ra2_Photon1Eta
    ra2_photonDR  = event.ra2_Photon1MinDR
    
    return (ra2_photonPt>minpt and ra2_photonDR > cutDR and abs(ra2_photonEta) < 2.5 and ((abs(ra2_photonEta) < 1.4442) or (abs(ra2_photonEta) > 1.566)))

####
def ra2MuonSelection(event,minpt,cutDR):
    ra2_nMuons  = event.ra2_nMuonsIso
    ra2_muon1Pt  = event.ra2_Muon1Pt
    ra2_muon1Eta = event.ra2_Muon1Eta
    ra2_muon2Pt  = event.ra2_Muon2Pt
    ra2_muon2Eta = event.ra2_Muon2Eta
    ra2_dimuonPt  = event.ra2_DiMuonPt
    ra2_dimuonEta = event.ra2_DiMuonEta
    ra2_dimuonDR  = event.ra2_DiMuonMinDR

    return ra2_nMuons > 1
    #return ( ra2_muon1DR > cutDR and ra2_muon2DR > cutDR and ra2_dimuonDR > cutDR and abs(ra2_muon1Eta) < 2.1 and abs(ra2_muon2Eta) < 2.1)
####
def ra2DPhiCuts(event):
    return (event.ra2_dPhiMHT1>0.5 and event.ra2_dPhiMHT2>0.5 and event.ra2_dPhiMHT3>0.3)
####
def ra2JetCuts(nJets,min,max):
    return (nJets>=min and nJets<=max)
####
def ra2HTCuts(htVal,min,max):
    return (htVal>min and htVal<max)
####
def ra2MHTCuts(mhtVal,min,max):
    return (mhtVal>min and mhtVal<max)
####
def genAcceptance(event,minpt):
    ra2_genBosonPt  = event.ra2_genBoson1Pt
    ra2_genBosonEta = event.ra2_genBoson1Eta

    return (ra2_genBosonPt>minpt and abs(ra2_genBosonEta) < 2.5 and ((abs(ra2_genBosonEta) < 1.4442) or (abs(ra2_genBosonEta) > 1.566)))
####
def genRecoIsoEff(event,minpt,recoPt,passReco,cutDR):
    ra2_nBosons  = event.ra2_nBosons
    ra2_bosonPt  = event.ra2_boson1Pt
    ra2_bosonEta = event.ra2_boson1Eta
    ra2_bosonDR  = event.ra2_boson1MinDR
    
    return (recoPt>minpt and passReco and ra2_bosonDR > cutDR and abs(ra2_bosonEta) < 2.5 and ((abs(ra2_bosonEta) < 1.4442) or (abs(ra2_bosonEta) > 1.566)))

####
def genCuts(event,selection,minpt,recoPt,passReco,nJets,jetCuts,htVal,htCuts,mhtCuts,cutDR):
    if nJets>=jetCuts[0] and nJets<=jetCuts[1]:
        if selection==0:
            return event.ra2_genBoson1Pt > minpt
        elif selection==1:
            return genAcceptance(event,minpt)
        elif selection==2:
            return genAcceptance(event,minpt) and genRecoIsoEff(event,minpt,recoPt,passReco,cutDR)
        elif selection==3:
            return genAcceptance(event,minpt) and genRecoIsoEff(event,minpt,recoPt,passReco,cutDR)
        elif selection==4:
            return genAcceptance(event,minpt) and genRecoIsoEff(event,minpt,recoPt,passReco,cutDR) and htVal>350 and event.ra2_genMHT>150 and passGenDPhiCuts(event)
        else:
            return False
    else:
        return False

####
def genCutsFull(event,selection,minpt,recoPt,passReco,nJets,jetCuts,htVal,htCuts,mhtCuts,cutDR):
    passDPhi = passGenDPhiCuts(event)

    passNJets  = False
    if ra2JetCuts(nJets,jetCuts[0],jetCuts[1]):
        passNJets  = True
    #
    passHT  = False
    if ra2HTCuts(htVal,htCuts[0],htCuts[1]):
        passHT  = True
    #
    passMHT = False
    if ra2MHTCuts(event.ra2_genMHT,mhtCuts[0],mhtCuts[1]):
        passMHT  = True

    #
    if passNJets and passHT and passMHT and passDPhi:
        if selection==0:
            return event.ra2_genBoson1Pt > 0
        elif selection==1:
            return event.ra2_genBoson1Pt > minpt
        elif selection==2:
            return genAcceptance(event,minpt)
        elif selection==3:
            return genAcceptance(event,minpt) and genRecoIsoEff(event,minpt,recoPt,passReco,cutDR)
        elif selection==4:
            return genAcceptance(event,minpt) and genRecoIsoEff(event,minpt,recoPt,passReco,cutDR)
        else:
            return False
    else:
        return False

#
def phenoCuts(event,cutLevel,minpt,nJets,jetCuts,htVal,htCuts,mhtCuts,cutDR):
    passDPhi = True
    passBase = True
    if cutLevel > 0:
        passDPhi = passGenDPhiCuts(event)

    #
    passNJets  = False
    if ra2JetCuts(nJets,jetCuts[0],jetCuts[1]):
        passNJets  = True
    #
    passHT  = False
    if ra2HTCuts(htVal,htCuts[0],htCuts[1]):
        passHT  = True
    #
    passMHT = False
    if ra2MHTCuts(event.ra2_genMHT,mhtCuts[0],mhtCuts[1]):
        passMHT  = True

    if passNJets and passHT and passMHT and passDPhi and passBase:
        return event.ra2_genBoson1Pt > minpt and event.ra2_genBoson1MinDR > cutDR
    else:
        return False
#
def passGenDPhiCuts(event):
    return event.ra2_genDPhiMHT1 > 0.5 and event.ra2_genDPhiMHT2 > 0.5 and event.ra2_genDPhiMHT3 > 0.3
