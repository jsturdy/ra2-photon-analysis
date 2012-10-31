import sys,os
import ROOT as r
from array import array
import math
import string
import re
#import argparse
import optparse

import analysisCuts as cutF

def main() :
    parser = optparse.OptionParser(description="Switch for data/MC running")
    parser.add_option('-m', action="store_true", default=False, dest="isMC")
    parser.add_option('-d', action="store_true", default=False, dest="debug")
    parser.add_option('-j', action="store",      default=5,     dest="numJets",type="int")
    parser.add_option('-f', action="store",      default=0,     dest="subsec", type="int")
    parser.add_option('-s', action="store",      default="zinv",dest="sample", type="string")
    parser.add_option('-t', action="store",      default="analysis5M",dest="treeName", type="string")
    options, args = parser.parse_args()

    r.gROOT.SetBatch(True)
    #debug = False
    myWorkingDir = os.getcwd()
    
    outFileName = "topTaggingPlots_%s_%s_min%d_job%d.root"%(options.sample,options.treeName,options.numJets,options.subsec)
    if options.debug:
        outFileName = "topTaggingPlots_%s_%s_min%d_test.root"%(options.sample,options.treeName,options.numJets)
        
    print outFileName
    outputFile = r.TFile(outFileName,"RECREATE")
    
    myChain = r.TChain('%s/RA2Values'%(options.treeName))
    
    if options.debug:
        if options.sample=="ttbar":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/ttjets_tunez2star_reco_tree_topTagged_csvt/res/*_?_?_???.root")
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht200_reco_tree_topTagged_csvt/res/*_?_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht400_reco_tree_topTagged_csvt/res/*_?_?_???.root")
        elif options.sample=="gjets":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht200_reco_tree_topTagged_csvt/res/*_?_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht400_reco_tree_topTagged_csvt/res/*_?_?_???.root")
        elif options.sample=="data":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012a_reco_tree_topTagged_csvt/res/*_?_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012b_reco_tree_topTagged_csvt/res/*_?_?_???.root")
    else:
        if options.sample=="ttbar":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/ttjets_tunez2star_reco_tree_topTagged_csvt/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht200_reco_tree_topTagged_csvt/res/%s.root"%(subfiles[options.subsec]))
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht400_reco_tree_topTagged_csvt/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="gjets":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht200_reco_tree_topTagged_csvt/res/%s.root"%(subfiles[options.subsec]))
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht400_reco_tree_topTagged_csvt/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="data":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012a_reco_tree_topTagged_csvt/res/%s.root"%(subfiles[options.subsec]))
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012b_reco_tree_topTagged_csvt/res/%s.root"%(subfiles[options.subsec]))
            

    baselineDir = outputFile.mkdir("baseline")
    invertedDir = outputFile.mkdir("inverted")
    baselineTopTagged = baselineDir.mkdir("topTagged")
    invertedTopTagged = invertedDir.mkdir("topTagged")

    analysisDirs = [
        [baselineDir,       "baseline" ],
        [invertedDir,       "inverted" ],
        [baselineTopTagged, "topTagged"],
        [invertedTopTagged, "topTagged"],
        ]
    
    nJets  = []
    bJets  = []
    jet1Pt = []
    jet2Pt = []
    jet3Pt = []
    jet4Pt = []
    met    = []
    
    bestTopJetMass = []
    mTbestTopJet   = []
    mTbJet         = []
    mT2            = []
    
    dPhi1    = []
    dPhi2    = []
    dPhi3    = []
    dPhi4    = []
    dPhiMin  = []
    dPhiMinB = []

    ###variables
    for d,dir in enumerate(analysisDirs):
        outputFile.cd()
        dir[0].cd()

        nJets .append(r.TH1D("h_nJets", "h_nJets", 15, -0.5, 14.5))
        bJets .append(r.TH1D("h_bJets", "h_bJets", 10,  -0.5, 9.5))
        jet1Pt.append(r.TH1D("h_jet1Pt","h_jet1Pt",2*50,  0, 1000))
        jet2Pt.append(r.TH1D("h_jet2Pt","h_jet2Pt",2*50,  0, 1000))
        jet3Pt.append(r.TH1D("h_jet3Pt","h_jet3Pt",2*50,  0, 1000))
        jet4Pt.append(r.TH1D("h_jet4Pt","h_jet4Pt",2*50,  0, 1000))
        met   .append(r.TH1D("h_met",   "h_met",   2*50,  0, 1000))
        
        bestTopJetMass.append(r.TH1D("h_bestTopJetMass", "h_bestTopJetMass", 2*50, 0, 1000))
        mTbestTopJet  .append(r.TH1D("h_mTbestTopJet"  , "h_mTbestTopJet"  , 2*75, 0, 1500))
        mTbJet        .append(r.TH1D("h_mTbJet"        , "h_mTbJet"        , 2*75, 0, 1500))
        mT2           .append(r.TH1D("h_mT2"           , "h_mT2"           , 2*50, 0, 1000))
    
        dPhi1   .append(r.TH1D("h_dPhi1","h_dPhi1",50,  0, 3.2))
        dPhi2   .append(r.TH1D("h_dPhi2","h_dPhi2",50,  0, 3.2))
        dPhi3   .append(r.TH1D("h_dPhi3","h_dPhi3",50,  0, 3.2))
        dPhi4   .append(r.TH1D("h_dPhi4","h_dPhi4",50,  0, 3.2))
        dPhiMin .append(r.TH1D("h_dPhiMin", "h_dPhiMin", 50, 0, 3.2))
        dPhiMinB.append(r.TH1D("h_dPhiMinB","h_dPhiMinB",50, 0, 3.2))

        
    ####
        nJets[d] .Sumw2()
        bJets[d] .Sumw2()
        jet1Pt[d].Sumw2()
        jet2Pt[d].Sumw2()
        jet3Pt[d].Sumw2()
        jet4Pt[d].Sumw2()
        met[d]   .Sumw2()

        bestTopJetMass[d].Sumw2()
        mTbestTopJet[d]  .Sumw2()
        mTbJet[d]        .Sumw2()
        mT2[d]           .Sumw2()

        dPhi1[d]   .Sumw2()
        dPhi2[d]   .Sumw2()
        dPhi3[d]   .Sumw2()
        dPhi4[d]   .Sumw2()
        dPhiMin[d] .Sumw2()
        dPhiMinB[d].Sumw2()

    ##################
    fChain = myChain

    ###Timing information
    decade  = 0
    century = 0
    tsw = r.TStopwatch()
    tenpcount = 1
    onepcount = 1


    nentries = fChain.GetEntries()
    print "nentries %d"%(nentries)
    sys.stdout.flush()
    i = 0
    #for i,event in enumerate(fChain):
    for event in fChain:

        # ==============print number of events done == == == == == == == =
        #double progress1 = 100.0 * jentry / (1.0 * nentries)
        #double progress10 = 10.0 * jentry / (1.0 * nentries)
        #int d = int (progress1)
        #int k = int (progress10)
        #if (d > century)
        #  cout << "."
        #century = d
        #if (k > decade)
        #  cout << 10 * k << " %" << endl
        #decade = k
        
        if ( i==0):
            tsw.Start()
            #print('.', end='')
            sys.stdout.write('.')
            sys.stdout.flush()
        if ((i*10)/nentries == tenpcount ) :
            tsw.Stop() 
            time = tsw.RealTime() 
            tsw.Start(r.kFALSE) 
            finTime = 0.
            frac = (i*1.0)/(nentries*1.0) 
            if (frac>0):
                finTime = time / frac - time 
                finMin = finTime / 60. 
                sys.stdout.write("%d%% done.  "%(tenpcount*10))
                #sys.stdout.write("t=7.2f"%(time))
                sys.stdout.write("t="+str(time))
                sys.stdout.write(" projected finish=%7d s("%(finTime))
                sys.stdout.write("%2.2f min).   "%(finMin))
                sys.stdout.write("\n")
                sys.stdout.flush()
                tenpcount = tenpcount + 1
        
        elif ( (i*100)/nentries == onepcount ) :
            #print('.', end='')
            sys.stdout.write('.')
            sys.stdout.flush()
            onepcount = onepcount + 1

        #sys.stdout.flush()
        #print "PU_et %2.2f"%(event.ra2_PUWt)
        ra2_nJetsPt30Eta24  = event.ra2_nJetsPt30Eta24 
        ra2_bJetsPt30Eta24  = event.ra2_bJetsPt30Eta24 
        ra2_PUWt            = event.ra2_PUWt           
        ra2_EventWt         = event.ra2_EventWt        
        ra2_MET             = event.ra2_MET            
        ra2_dPhi1           = event.ra2_dPhiMET1         
        ra2_dPhi2           = event.ra2_dPhiMET2         
        ra2_dPhi3           = event.ra2_dPhiMET3         
        ra2_dPhi4           = event.ra2_dPhiMET4         
        ra2_dPhiMin         = event.ra2_dPhiMETMin         
        ra2_dPhiMinB        = event.ra2_dPhiMETMinB         
        ra2_Jet1Pt          = event.ra2_Jet1Pt         
        ra2_Jet1Eta         = event.ra2_Jet1Eta        
        ra2_Jet2Pt          = event.ra2_Jet2Pt         
        ra2_Jet2Eta         = event.ra2_Jet2Eta        
        ra2_Jet3Pt          = event.ra2_Jet3Pt         
        ra2_Jet3Eta         = event.ra2_Jet3Eta        
        ra2_Jet4Pt          = event.ra2_Jet4Pt         
        ra2_Jet4Eta         = event.ra2_Jet4Eta        

        ##top tagger variables 
        ra2_bestTopJetMass  = event.ra2_bestTopJetMass 
        ra2_TbestTopJet     = event.ra2_TbestTopJet    
        ra2_TbJet           = event.ra2_TbJet          
        ra2_MT2             = event.ra2_MT2            

        ####Fill plots
        for d,dir in enumerate(analysisDirs):
            outputFile.cd()
            dir[0].cd()
            fillPlots = False
            if d == 0 and cutF.hadStopBaseline(event,1,options.numJets,False):
                fillPlots = True
            if d == 1 and cutF.hadStopBaseline(event,1,options.numJets,True):
                fillPlots = True
            if d == 2 and cutF.hadStopBaseline(event,1,options.numJets,False):
                if cutF.topTaggerCuts(event):
                    fillPlots = True
            if d == 3 and cutF.hadStopBaseline(event,1,options.numJets,True):
                if cutF.topTaggerCuts(event):
                    fillPlots = True
            if fillPlots:
                nJets[d] .Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
                bJets[d] .Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                jet1Pt[d].Fill(ra2_Jet1Pt        ,ra2_EventWt)
                jet2Pt[d].Fill(ra2_Jet2Pt        ,ra2_EventWt)
                jet3Pt[d].Fill(ra2_Jet3Pt        ,ra2_EventWt)
                jet4Pt[d].Fill(ra2_Jet4Pt        ,ra2_EventWt)
                met[d]   .Fill(ra2_MET           ,ra2_EventWt)
                
                bestTopJetMass[d].Fill(ra2_bestTopJetMass,ra2_EventWt)
                mTbestTopJet[d]  .Fill(ra2_TbestTopJet   ,ra2_EventWt)
                mTbJet[d]        .Fill(ra2_TbJet         ,ra2_EventWt)
                mT2[d]           .Fill(ra2_MT2           ,ra2_EventWt)
                
                dPhi1[d]   .Fill(ra2_dPhi1         ,ra2_EventWt)
                dPhi2[d]   .Fill(ra2_dPhi2         ,ra2_EventWt)
                dPhi3[d]   .Fill(ra2_dPhi3         ,ra2_EventWt)
                dPhi4[d]   .Fill(ra2_dPhi4         ,ra2_EventWt)
                dPhiMin[d] .Fill(ra2_dPhiMin       ,ra2_EventWt)
                dPhiMinB[d].Fill(ra2_dPhiMinB      ,ra2_EventWt)
                
        ### at least one b-jet
                
        #########
        i = i + 1

    #####        
    outputFile.cd()

    ####
    for d,dir in enumerate(analysisDirs):
        outputFile.cd()
        dir[0].cd()

        nJets[d] .Write()
        bJets[d] .Write()
        jet1Pt[d].Write()
        jet2Pt[d].Write()
        jet3Pt[d].Write()
        jet4Pt[d].Write()
        met[d]   .Write()

        bestTopJetMass[d].Write()
        mTbestTopJet[d]  .Write()
        mTbJet[d]        .Write()
        mT2[d]           .Write()

        dPhi1[d].Write()
        dPhi2[d].Write()
        dPhi3[d].Write()
        dPhi4[d].Write()
        dPhiMin[d] .Write()
        dPhiMinB[d].Write()


    outputFile.Write()
    outputFile.Close()

##########################
if __name__ == '__main__':
    main()
    
