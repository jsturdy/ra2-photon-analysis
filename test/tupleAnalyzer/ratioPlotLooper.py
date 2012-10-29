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
    ##parser = argparse.ArgumentParser(description="Switch for data/MC running")
    ##parser.add_argument('-m', action="store_true", dest="isMC", default=False)
    ##parser.add_argument('-d', action="store_true", dest="debug", default=False)
    ##results = parser.parse_args()

    parser = optparse.OptionParser(description="Switch for data/MC running")
    parser.add_option('-m', action="store_true", default=False, dest="isMC")
    parser.add_option('-d', action="store_true", default=False, dest="debug")
    options, args = parser.parse_args()

    r.gROOT.SetBatch(True)
    #debug = False
    #isMC  = True
    myWorkingDir = os.getcwd()
    suffix = "data"
    if options.isMC:
        suffix = "mc"
    outFileName = "photonPlots_%s.root"%(suffix)
    if options.debug:
        outFileName = "photonPlots_%s_test.root"%(suffix)
        print "output file %s"%(outFileName)

    outputFile = r.TFile(outFileName,"RECREATE")

    myChain = r.TChain('analysisIDPFIso/RA2Values')
    if options.isMC:
        if options.debug:
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht200_reco_tree_topTagged_v3/res/*_?_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht400_reco_tree_topTagged_v3/res/*_?_?_???.root")
        else:
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht200_reco_tree_topTagged_v3/res/*.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/gjetsht400_reco_tree_topTagged_v3/res/*.root")
    else:
        if options.debug:
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012a_tree_topTagged_v2/res/*_?_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012b_tree_topTagged_v2/res/*_?_?_???.root")
        else:
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012a_tree_topTagged_v2/res/*.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/run2012b_tree_topTagged_v2/res/*.root")
    ###variables
    nVertices = r.TH1D("h_nVertices","h_nVertices",50, -0.5, 49.5)
    nJetsAll = r.TH1D("h_nJetsAll","h_nJetsAll",15, -0.5, 14.5)
    phot1Pt = r.TH1D("h_phot1Pt","h_phot1Pt",50,  0, 500)
    jet1Pt = r.TH1D("h_jet1Pt","h_jet1Pt",50,  0, 1000)
    jet2Pt = r.TH1D("h_jet2Pt","h_jet2Pt",50,  0, 500)
    jet3Pt = r.TH1D("h_jet3Pt","h_jet3Pt",50,  0, 500)
    jet4Pt = r.TH1D("h_jet4Pt","h_jet4Pt",50,  0, 500)
    ht     = r.TH1D("h_ht","h_ht",30,  0, 3000)
    ht2    = r.TH1D("h_ht2","h_ht2",30,  0, 3000)
    ht3    = r.TH1D("h_ht3","h_ht3",30,  0, 3000)
    mht    = r.TH1D("h_mht","h_mht",20,  0, 1000)
    mht2   = r.TH1D("h_mht2","h_mht2",20,  0, 1000)
    mht3   = r.TH1D("h_mht3","h_mht3",20,  0, 1000)
    dPhi1  = r.TH1D("h_dPhi1","h_dPhi1",50,  0, 3.2)
    dPhi2  = r.TH1D("h_dPhi2","h_dPhi2",50,  0, 3.2)
    dPhi3  = r.TH1D("h_dPhi3","h_dPhi3",50,  0, 3.2)

    nVerticesReWt = r.TH1D("h_nVerticesReWt","h_nVerticesReWt",50, -0.5, 49.5)
    nJetsAllPU = r.TH1D("h_nJetsAllPU","h_nJetsAllPU",15, -0.5, 14.5)
    phot1PtPU = r.TH1D("h_phot1PtPU","h_phot1PtPU",50,  0, 500)
    jet1PtPU = r.TH1D("h_jet1PtPU","h_jet1PtPU",100,  0, 100)
    jet2PtPU = r.TH1D("h_jet2PtPU","h_jet2PtPU",50,  0, 500)
    jet3PtPU = r.TH1D("h_jet3PtPU","h_jet3PtPU",50,  0, 500)
    jet4PtPU = r.TH1D("h_jet4PtPU","h_jet4PtPU",50,  0, 500)
    htPU     = r.TH1D("h_htPU","h_htPU",30,  0, 3000)
    ht2PU    = r.TH1D("h_ht2PU","h_ht2PU",30,  0, 3000)
    ht3PU    = r.TH1D("h_ht3PU","h_ht3PU",30,  0, 3000)
    mhtPU    = r.TH1D("h_mhtPU","h_mhtPU",20,  0, 1000)
    mht2PU   = r.TH1D("h_mht2PU","h_mht2PU",20,  0, 1000)
    mht3PU   = r.TH1D("h_mht3PU","h_mht3PU",20,  0, 1000)
    dPhi1PU  = r.TH1D("h_dPhi1PU","h_dPhi1PU",50,  0, 3.2)
    dPhi2PU  = r.TH1D("h_dPhi2PU","h_dPhi2PU",50,  0, 3.2)
    dPhi3PU  = r.TH1D("h_dPhi3PU","h_dPhi3PU",50,  0, 3.2)


    nVertices    .Sumw2()
    nJetsAll.Sumw2()
    phot1Pt.Sumw2()
    jet1Pt.Sumw2()
    jet2Pt.Sumw2()
    jet3Pt.Sumw2()
    jet4Pt.Sumw2()
    ht   .Sumw2()
    ht2  .Sumw2()
    ht3  .Sumw2()
    mht  .Sumw2()
    mht2 .Sumw2()
    mht3 .Sumw2()
    dPhi1.Sumw2()
    dPhi2.Sumw2()
    dPhi3.Sumw2()

    nVerticesReWt.Sumw2()
    nJetsAllPU.Sumw2()
    phot1PtPU.Sumw2()
    jet1PtPU.Sumw2()
    jet2PtPU.Sumw2()
    jet3PtPU.Sumw2()
    jet4PtPU.Sumw2()
    htPU   .Sumw2()
    ht2PU  .Sumw2()
    ht3PU  .Sumw2()
    mhtPU  .Sumw2()
    mht2PU .Sumw2()
    mht3PU .Sumw2()
    dPhi1PU.Sumw2()
    dPhi2PU.Sumw2()
    dPhi3PU.Sumw2()

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
        ra2_Vertices        = event.ra2_Vertices      
        ra2_nJetsPt50Eta25  = event.ra2_nJetsPt50Eta25 
        ra2_PUWt            = event.ra2_PUWt           
        ra2_EventWt         = event.ra2_EventWt        
        ra2_HT              = event.ra2_HT            
        ra2_MHT             = event.ra2_MHT            
        ra2_dPhi1           = event.ra2_dPhiMHT1
        ra2_dPhi2           = event.ra2_dPhiMHT2
        ra2_dPhi3           = event.ra2_dPhiMHT3
        ra2_Jet1Pt          = event.ra2_Jet1Pt         
        ra2_Jet1Eta         = event.ra2_Jet1Eta        
        ra2_Jet2Pt          = event.ra2_Jet2Pt         
        ra2_Jet2Eta         = event.ra2_Jet2Eta        
        ra2_Jet3Pt          = event.ra2_Jet3Pt         
        ra2_Jet3Eta         = event.ra2_Jet3Eta        
        ra2_Jet4Pt          = event.ra2_Jet4Pt         
        ra2_Jet4Eta         = event.ra2_Jet4Eta        
        ra2_nPhotonsIso     = event.ra2_nPhotonsIso
        ra2_Photon1Pt       = event.ra2_Photon1Pt
        ra2_Photon70PFHT400 = event.ra2_Photon70PFHT400
        ra2_Photon70PFNoPUHT400 = event.ra2_Photon70PFNoPUHT400

        triggers = (ra2_Photon70PFHT400 or ra2_Photon70PFNoPUHT400)
        if options.isMC:
            triggers = True
        if (triggers) and (ra2_nPhotonsIso>0 and ra2_Photon1Pt>100) and cutF.ra2Baseline(event):
        ####Fill plots
            nVertices.Fill(ra2_Vertices,ra2_EventWt)
            nJetsAll.Fill(ra2_nJetsPt50Eta25,ra2_EventWt)
            phot1Pt.Fill(ra2_Photon1Pt,ra2_EventWt)
            jet1Pt.Fill(ra2_Jet1Pt,ra2_EventWt)
            jet2Pt.Fill(ra2_Jet2Pt,ra2_EventWt)
            jet3Pt.Fill(ra2_Jet3Pt,ra2_EventWt)
            jet4Pt.Fill(ra2_Jet4Pt,ra2_EventWt)
            ht   .Fill(ra2_HT,ra2_EventWt)
            mht  .Fill(ra2_MHT,ra2_EventWt)
            dPhi1.Fill(ra2_dPhi1,ra2_EventWt)
            dPhi2.Fill(ra2_dPhi2,ra2_EventWt)
            dPhi3.Fill(ra2_dPhi3,ra2_EventWt)
            
            nVerticesReWt.Fill(ra2_Vertices,ra2_EventWt*ra2_PUWt)
            nJetsAllPU.Fill(ra2_nJetsPt50Eta25,ra2_EventWt*ra2_PUWt)
            phot1PtPU.Fill(ra2_Photon1Pt,ra2_EventWt*ra2_PUWt)
            jet1PtPU.Fill(ra2_Jet1Pt,ra2_EventWt*ra2_PUWt)
            jet2PtPU.Fill(ra2_Jet2Pt,ra2_EventWt*ra2_PUWt)
            jet3PtPU.Fill(ra2_Jet3Pt,ra2_EventWt*ra2_PUWt)
            jet4PtPU.Fill(ra2_Jet4Pt,ra2_EventWt*ra2_PUWt)
            htPU   .Fill(ra2_HT,ra2_EventWt*ra2_PUWt)
            mhtPU  .Fill(ra2_MHT,ra2_EventWt*ra2_PUWt)
            dPhi1PU.Fill(ra2_dPhi1,ra2_EventWt*ra2_PUWt)
            dPhi2PU.Fill(ra2_dPhi2,ra2_EventWt*ra2_PUWt)
            dPhi3PU.Fill(ra2_dPhi3,ra2_EventWt*ra2_PUWt)

            if ra2_nJetsPt50Eta25 == 2:
                ht2   .Fill(ra2_HT,ra2_EventWt)
                mht2  .Fill(ra2_MHT,ra2_EventWt)
                ht2PU .Fill(ra2_HT,ra2_EventWt*ra2_PUWt)
                mht2PU.Fill(ra2_MHT,ra2_EventWt*ra2_PUWt)

            if ra2_nJetsPt50Eta25 > 2:
                ht3   .Fill(ra2_HT,ra2_EventWt)
                mht3  .Fill(ra2_MHT,ra2_EventWt)
                ht3PU .Fill(ra2_HT,ra2_EventWt*ra2_PUWt)
                mht3PU.Fill(ra2_MHT,ra2_EventWt*ra2_PUWt)

        #########
        i = i + 1

    #####        
    outputFile.cd()
    nVertices    .Write()
    nJetsAll.Write()
    phot1Pt.Write()
    jet1Pt.Write()
    jet2Pt.Write()
    jet3Pt.Write()
    jet4Pt.Write()
    ht   .Write()
    ht2  .Write()
    ht3  .Write()
    mht  .Write()
    mht2 .Write()
    mht3 .Write()
    dPhi1.Write()
    dPhi2.Write()
    dPhi3.Write()

    nVerticesReWt.Write()
    nJetsAllPU.Write()
    phot1PtPU.Write()
    jet1PtPU.Write()
    jet2PtPU.Write()
    jet3PtPU.Write()
    jet4PtPU.Write()
    htPU   .Write()
    ht2PU  .Write()
    ht3PU  .Write()
    mhtPU  .Write()
    mht2PU .Write()
    mht3PU .Write()
    dPhi1PU.Write()
    dPhi2PU.Write()
    dPhi3PU.Write()

    outputFile.Write()
    outputFile.Close()

##########################
if __name__ == '__main__':
    main()
    
