import sys,os
import ROOT as r
from array import array
import math
import string
import re
#import argparse
import optparse

import itertools

import analysisCuts as cutF

def main() :
    parser = optparse.OptionParser(description="Switch for data/MC running")
    parser.add_option('-m', action="store_true", default=False, dest="isMC")
    parser.add_option('-d', action="store_true", default=False, dest="debug")
    parser.add_option('-s', action="store",      default="data",dest="sample", type="string")
    parser.add_option('-j', action="store",      default=5,     dest="numJets",type="int")
    parser.add_option('-t', action="store",      default="analysisIDPFIso5T",dest="treeName", type="string")
    parser.add_option('-c', action="store",      default=100.0,dest="minpt", type="float")
    options, args = parser.parse_args()

    r.gROOT.SetBatch(True)
    #debug = False
    myWorkingDir = os.getcwd()

#    htBins = [
#        [500,750,1000,1250,1500,10000],
#        [500,750,1000,1250,1500,10000],
#        [500,750,1000,1250,10000],
#        [500,750,1000,1250,10000],
#        ]
#    mhtBins = [
#        [[200,300,450,600,10000],
#        [100,200,300,450,600,10000],
#        [100,200,300,450,600,10000],
#        [100,200,300,450,10000],
#        [100,200,350,10000]],
#        
#        [[200,300,450,600,10000],
#        [100,200,300,450,600,10000],
#        [100,200,300,450,600,10000],
#        [100,200,300,450,10000],
#        [100,200,350,10000]],
#
#        [[200,300,450,10000],
#        [100,200,300,450,10000],
#        [100,200,300,10000],
#        [100,200,300,10000]],
#
#        [[200,300,10000],
#        [100,200,10000],
#        [100,200,10000],
#        [100,200,10000]],
#
#        ]

    
    outFileName = "%sYieldPlots.root"%(options.sample)
    if options.debug:
        outFileName = "%sYieldPlots_test.root"%(options.sample)

    print outFileName
    sys.stdout.flush()
    outputFile = r.TFile(outFileName,"RECREATE")

    inclusiveDir = outputFile.mkdir("inclusive")
    jet2Dir      = outputFile.mkdir("2jets")
    jet3Dir      = outputFile.mkdir("3jets")
    jet4Dir      = outputFile.mkdir("4jets")
    jet5Dir      = outputFile.mkdir("5jets")
    jet6Dir      = outputFile.mkdir("6jets")
    jet7Dir      = outputFile.mkdir("7jets")
    jet8Dir      = outputFile.mkdir("8jets")
    jet3to5Dir   = outputFile.mkdir("3to5jets")
    jet6to7Dir   = outputFile.mkdir("6to7jets")

    njbins = [
        ["inclusive",[2,100]],
        ["2jets"    ,[2,2]],
        ["3jets"    ,[3,3]],
        ["4jets"    ,[4,4]],
        ["5jets"    ,[5,5]],
        ["6jets"    ,[6,6]],
        ["7jets"    ,[7,7]],
        ["8jets"    ,[8,100]],
        ["3to5jets" ,[3,5]],
        ["6to7jets" ,[6,7]],
        
        ]
    fChain = r.TChain('%s/RA2Values'%(options.treeName))

    if options.debug:
        if options.sample=="gjets":
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht200_tree_topTagged_csvt/res/*_?_?_???.root")
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht400_tree_topTagged_csvt/res/*_?_?_???.root")
        elif options.sample=="zinv":
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht200_tree_topTagged_csvt/res/*_?_?_???.root")
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht400_tree_topTagged_csvt/res/*_?_?_???.root")
        elif options.sample=="data":
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012a_tree_topTagged_csvt/res/*_?_?_???.root")
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012b_tree_topTagged_csvt/res/*_?_?_???.root")
    else:
        if options.sample=="gjets":
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht200_tree_topTagged_csvt/res/*.root")
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht400_tree_topTagged_csvt/res/*.root")
        elif options.sample=="zinv":
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht200_tree_topTagged_csvt/res/*.root")
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht400_tree_topTagged_csvt/res/*.root")
        elif options.sample=="data":
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012a_tree_topTagged_csvt/res/*.root")
            fChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012b_tree_topTagged_csvt/res/*.root")
    
    ###variables
    htvsmht = []
    htminvvsmht = []

    for j,jetdir in enumerate(njbins):
        outputFile.cd()
        outputFile.cd(jetdir[0])
        
        htvsmht    .append(r.TH2D("h_%s_htvsmht"%(jetdir[0])    ,"h_%s_htvsmht"%(jetdir[0])    ,60,  0, 3000, 20, 0, 1000))
        htminvvsmht.append(r.TH2D("h_%s_htminvvsmht"%(jetdir[0]),"h_%s_htminvvsmht"%(jetdir[0]),60,  0, 3000, 20, 0, 1000))
    ####
        htvsmht[j]    .Sumw2()
        htminvvsmht[j].Sumw2()



    ##################
    #fChain = myChain

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
        ra2_Vertices            = event.ra2_Vertices
        ra2_nJetsPt50Eta25      = event.ra2_nJetsPt50Eta25
        ra2_nJetsPt50Eta25MInv  = event.ra2_nJetsPt50Eta25MInv
        ra2_PUWt                = event.ra2_PUWt
        ra2_EventWt             = event.ra2_EventWt
        ra2_HT                  = event.ra2_HT
        ra2_HTMInv              = event.ra2_HTMInv
        ra2_MHT                 = event.ra2_MHT
        ra2_dPhi1               = event.ra2_dPhiMHT1
        ra2_dPhi2               = event.ra2_dPhiMHT2
        ra2_dPhi3               = event.ra2_dPhiMHT3
        ra2_nPhotonsIso         = event.ra2_nPhotonsIso
        ra2_Photon1Pt           = event.ra2_Photon1Pt
        ra2_Photon70PFHT400     = event.ra2_Photon70PFHT400
        ra2_Photon70PFNoPUHT400 = event.ra2_Photon70PFNoPUHT400

        triggers = (ra2_Photon70PFHT400 or ra2_Photon70PFNoPUHT400)
        if options.isMC:
            triggers = True

        if (triggers) and (ra2_nPhotonsIso>0 and ra2_Photon1Pt>100) and cutF.ra2DPhiCuts(event):
        ####Fill plots
            for j,jetdir in enumerate(njbins):
                outputFile.cd()
                outputFile.cd(jetdir[0])

                if cutF.ra2JetCuts(ra2_nJetsPt50Eta25,jetdir[1][0],jetdir[1][1]):
                    htvsmht[j].Fill(ra2_HT,ra2_MHT    ,ra2_EventWt)

                ##Passing MInv cut
                if cutF.ra2JetCuts(ra2_nJetsPt50Eta25MInv,jetdir[1][0],jetdir[1][1]):
                    htminvvsmht[j].Fill(ra2_HTMInv,ra2_MHT,ra2_EventWt)
                
        #########
        i = i + 1

    #####        
    outputFile.cd()

    ####
    for j,jetdir in enumerate(njbins):
        outputFile.cd()
        outputFile.cd(jetdir[0])
        
        htvsmht[j].Write()
        htminvvsmht[j].Write()
        
    outputFile.Write()
    outputFile.Close()

##########################
if __name__ == '__main__':
    main()
    
