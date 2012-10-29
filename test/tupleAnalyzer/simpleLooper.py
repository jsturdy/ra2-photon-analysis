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
    options, args = parser.parse_args()

    r.gROOT.SetBatch(True)
    #debug = False
    myWorkingDir = os.getcwd()
    
    outFileName = "topTaggingPlots_zinv.root"
    #outFileName = "topTaggingPlots_DY_ht_no_muon_removal.root"
    #outFileName = "topTaggingPlots_DY_ht_dimuons_removed.root"
    if options.debug:
        outFileName = "topTaggingPlots_test.root"

    outputFile = r.TFile(outFileName,"RECREATE")

    myChain = r.TChain('analysis/RA2Values')
    #myChain = r.TChain('analysisWithMuon/RA2Values')
    #myChain = r.TChain('analysisNoMuon/RA2Values')
    if options.debug:
        myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht200_reco_tree_topTagged_v2/res/*_?_?_???.root")
        myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht400_reco_tree_topTagged_v2/res/*_?_?_???.root")
    else:
        myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht200_reco_tree_topTagged_v2/res/*.root")
        myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht400_reco_tree_topTagged_v2/res/*.root")
    
##    if options.debug:
##        myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/dyjetstoll_ht200_reco_tree_topTagged_v2/res/*_?_?_???.root")
##        myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/dyjetstoll_ht400_reco_tree_topTagged_v2/res/*_?_?_???.root")
##    else:
##        myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/dyjetstoll_ht200_reco_tree_topTagged_v2/res/*.root")
##        myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/dyjetstoll_ht400_reco_tree_topTagged_v2/res/*.root")
    
    ###variables
    nVertices = r.TH1D("h_nVertices","h_nVertices",50, -0.5, 49.5)
    nVerticesReWt = r.TH1D("h_nVerticesReWt","h_nVerticesReWt",50, -0.5, 49.5)
    nJetsAll = r.TH1D("h_nJetsAll","h_nJetsAll",15, -0.5, 14.5)
    bJetsAll = r.TH1D("h_bJetsAll","h_bJetsAll",5,  -0.5, 4.5)
    jet1Pt = r.TH1D("h_jet1Pt","h_jet1Pt",50,  0, 1000)
    jet2Pt = r.TH1D("h_jet2Pt","h_jet2Pt",50,  0, 500)
    met    = r.TH1D("h_met","h_met",50,  0, 1000)
    mht    = r.TH1D("h_mht","h_mht",50,  0, 1000)
    dPhi1  = r.TH1D("h_dPhi1","h_dPhi1",50,  0, 3.2)
    dPhi2  = r.TH1D("h_dPhi2","h_dPhi2",50,  0, 3.2)
    dPhi3  = r.TH1D("h_dPhi3","h_dPhi3",50,  0, 3.2)

    nJetsAllPU = r.TH1D("h_nJetsAllPU","h_nJetsAllPU",15, -0.5, 14.5)
    bJetsAllPU = r.TH1D("h_bJetsAllPU","h_bJetsAllPU",5,  -0.5, 4.5)
    jet1PtPU = r.TH1D("h_jet1PtPU","h_jet1PtPU",50,  0, 1000)
    jet2PtPU = r.TH1D("h_jet2PtPU","h_jet2PtPU",50,  0, 500)
    metPU    = r.TH1D("h_metPU","h_metPU",50,  0, 1000)
    mhtPU    = r.TH1D("h_mhtPU","h_mhtPU",50,  0, 1000)
    dPhi1PU  = r.TH1D("h_dPhi1PU","h_dPhi1PU",50,  0, 3.2)
    dPhi2PU  = r.TH1D("h_dPhi2PU","h_dPhi2PU",50,  0, 3.2)
    dPhi3PU  = r.TH1D("h_dPhi3PU","h_dPhi3PU",50,  0, 3.2)

    ####geq 4 jets
    nJets4 = r.TH1D("h_nJets4","h_nJets4",15, -0.5, 14.5)
    bJets4 = r.TH1D("h_bJets4","h_bJets4",5,  -0.5, 4.5)
    h4_jet1Pt = r.TH1D("h4_jet1Pt","h4_jet1Pt",50,  0, 1000)
    h4_jet2Pt = r.TH1D("h4_jet2Pt","h4_jet2Pt",50,  0, 500)
    h4_met    = r.TH1D("h4_met","h4_met",50,  0, 1000)
    h4_mht    = r.TH1D("h4_mht","h4_mht",50,  0, 1000)
    h4_dPhi1  = r.TH1D("h4_dPhi1","h4_dPhi1",50,  0, 3.2)
    h4_dPhi2  = r.TH1D("h4_dPhi2","h4_dPhi2",50,  0, 3.2)
    h4_dPhi3  = r.TH1D("h4_dPhi3","h4_dPhi3",50,  0, 3.2)
    nJets4TopTagged = r.TH1D("h_nJets4TopTagged","h_nJets4TopTagged",15, -0.5, 14.5)
    bJets4TopTagged = r.TH1D("h_bJets4TopTagged","h_bJets4TopTagged",5,  -0.5, 4.5)
    h4tt_jet1Pt = r.TH1D("h4tt_jet1Pt","h4tt_jet1Pt",50,  0, 1000)
    h4tt_jet2Pt = r.TH1D("h4tt_jet2Pt","h4tt_jet2Pt",50,  0, 500)
    h4tt_met    = r.TH1D("h4tt_met","h4tt_met",50,  0, 1000)
    h4tt_mht    = r.TH1D("h4tt_mht","h4tt_mht",50,  0, 1000)
    h4tt_dPhi1  = r.TH1D("h4tt_dPhi1","h4tt_dPhi1",50,  0, 3.2)
    h4tt_dPhi2  = r.TH1D("h4tt_dPhi2","h4tt_dPhi2",50,  0, 3.2)
    h4tt_dPhi3  = r.TH1D("h4tt_dPhi3","h4tt_dPhi3",50,  0, 3.2)

    nJets4PU = r.TH1D("h_nJets4PU","h_nJets4PU",15, -0.5, 14.5)
    bJets4PU = r.TH1D("h_bJets4PU","h_bJets4PU",5,  -0.5, 4.5)
    h4_jet1PtPU = r.TH1D("h4_jet1PtPU","h4_jet1PtPU",50,  0, 1000)
    h4_jet2PtPU = r.TH1D("h4_jet2PtPU","h4_jet2PtPU",50,  0, 500)
    h4_metPU    = r.TH1D("h4_metPU","h4_metPU",50,  0, 1000)
    h4_mhtPU    = r.TH1D("h4_mhtPU","h4_mhtPU",50,  0, 1000)
    h4_dPhi1PU  = r.TH1D("h4_dPhi1PU","h4_dPhi1PU",50,  0, 3.2)
    h4_dPhi2PU  = r.TH1D("h4_dPhi2PU","h4_dPhi2PU",50,  0, 3.2)
    h4_dPhi3PU  = r.TH1D("h4_dPhi3PU","h4_dPhi3PU",50,  0, 3.2)
    nJets4PUTopTagged = r.TH1D("h_nJets4PUTopTagged","h_nJets4PUTopTagged",15, -0.5, 14.5)
    bJets4PUTopTagged = r.TH1D("h_bJets4PUTopTagged","h_bJets4PUTopTagged",5,  -0.5, 4.5)
    h4tt_jet1PtPU = r.TH1D("h4tt_jet1PtPU","h4tt_jet1PtPU",50,  0, 1000)
    h4tt_jet2PtPU = r.TH1D("h4tt_jet2PtPU","h4tt_jet2PtPU",50,  0, 500)
    h4tt_metPU    = r.TH1D("h4tt_metPU","h4tt_metPU",50,  0, 1000)
    h4tt_mhtPU    = r.TH1D("h4tt_mhtPU","h4tt_mhtPU",50,  0, 1000)
    h4tt_dPhi1PU  = r.TH1D("h4tt_dPhi1PU","h4tt_dPhi1PU",50,  0, 3.2)
    h4tt_dPhi2PU  = r.TH1D("h4tt_dPhi2PU","h4tt_dPhi2PU",50,  0, 3.2)
    h4tt_dPhi3PU  = r.TH1D("h4tt_dPhi3PU","h4tt_dPhi3PU",50,  0, 3.2)

    nJetsNobReq4 = r.TH1D("h_nJetsNobReq4","h_nJetsNobReq4",15, -0.5, 14.5)
    bJetsNobReq4 = r.TH1D("h_bJetsNobReq4","h_bJetsNobReq4",5,  -0.5, 4.5)
    h4nob_jet1Pt = r.TH1D("h4nob_jet1Pt","h4nob_jet1Pt",50,  0, 1000)
    h4nob_jet2Pt = r.TH1D("h4nob_jet2Pt","h4nob_jet2Pt",50,  0, 500)
    h4nob_met    = r.TH1D("h4nob_met","h4nob_met",50,  0, 1000)
    h4nob_mht    = r.TH1D("h4nob_mht","h4nob_mht",50,  0, 1000)
    h4nob_dPhi1  = r.TH1D("h4nob_dPhi1","h4nob_dPhi1",50,  0, 3.2)
    h4nob_dPhi2  = r.TH1D("h4nob_dPhi2","h4nob_dPhi2",50,  0, 3.2)
    h4nob_dPhi3  = r.TH1D("h4nob_dPhi3","h4nob_dPhi3",50,  0, 3.2)
    nJetsNobReq4TopTagged = r.TH1D("h_nJetsNobReq4TopTagged","h_nJetsNobReq4TopTagged",15, -0.5, 14.5)
    bJetsNobReq4TopTagged = r.TH1D("h_bJetsNobReq4TopTagged","h_bJetsNobReq4TopTagged",5,  -0.5, 4.5)
    h4nobtt_jet1Pt = r.TH1D("h4nobtt_jet1Pt","h4nobtt_jet1Pt",50,  0, 1000)
    h4nobtt_jet2Pt = r.TH1D("h4nobtt_jet2Pt","h4nobtt_jet2Pt",50,  0, 500)
    h4nobtt_met    = r.TH1D("h4nobtt_met","h4nobtt_met",50,  0, 1000)
    h4nobtt_mht    = r.TH1D("h4nobtt_mht","h4nobtt_mht",50,  0, 1000)
    h4nobtt_dPhi1  = r.TH1D("h4nobtt_dPhi1","h4nobtt_dPhi1",50,  0, 3.2)
    h4nobtt_dPhi2  = r.TH1D("h4nobtt_dPhi2","h4nobtt_dPhi2",50,  0, 3.2)
    h4nobtt_dPhi3  = r.TH1D("h4nobtt_dPhi3","h4nobtt_dPhi3",50,  0, 3.2)

    nJetsNobReq4PU = r.TH1D("h_nJetsNobReq4PU","h_nJetsNobReq4PU",15, -0.5, 14.5)
    bJetsNobReq4PU = r.TH1D("h_bJetsNobReq4PU","h_bJetsNobReq4PU",5,  -0.5, 4.5)
    h4nob_jet1PtPU = r.TH1D("h4nob_jet1PtPU","h4nob_jet1PtPU",50,  0, 1000)
    h4nob_jet2PtPU = r.TH1D("h4nob_jet2PtPU","h4nob_jet2PtPU",50,  0, 500)
    h4nob_metPU    = r.TH1D("h4nob_metPU","h4nob_metPU",50,  0, 1000)
    h4nob_mhtPU    = r.TH1D("h4nob_mhtPU","h4nob_mhtPU",50,  0, 1000)
    h4nob_dPhi1PU  = r.TH1D("h4nob_dPhi1PU","h4nob_dPhi1PU",50,  0, 3.2)
    h4nob_dPhi2PU  = r.TH1D("h4nob_dPhi2PU","h4nob_dPhi2PU",50,  0, 3.2)
    h4nob_dPhi3PU  = r.TH1D("h4nob_dPhi3PU","h4nob_dPhi3PU",50,  0, 3.2)
    nJetsNobReq4PUTopTagged = r.TH1D("h_nJetsNobReq4PUTopTagged","h_nJetsNobReq4PUTopTagged",15, -0.5, 14.5)
    bJetsNobReq4PUTopTagged = r.TH1D("h_bJetsNobReq4PUTopTagged","h_bJetsNobReq4PUTopTagged",5,  -0.5, 4.5)
    h4nobtt_jet1PtPU = r.TH1D("h4nobtt_jet1PtPU","h4nobtt_jet1PtPU",50,  0, 1000)
    h4nobtt_jet2PtPU = r.TH1D("h4nobtt_jet2PtPU","h4nobtt_jet2PtPU",50,  0, 500)
    h4nobtt_metPU    = r.TH1D("h4nobtt_metPU","h4nobtt_metPU",50,  0, 1000)
    h4nobtt_mhtPU    = r.TH1D("h4nobtt_mhtPU","h4nobtt_mhtPU",50,  0, 1000)
    h4nobtt_dPhi1PU  = r.TH1D("h4nobtt_dPhi1PU","h4nobtt_dPhi1PU",50,  0, 3.2)
    h4nobtt_dPhi2PU  = r.TH1D("h4nobtt_dPhi2PU","h4nobtt_dPhi2PU",50,  0, 3.2)
    h4nobtt_dPhi3PU  = r.TH1D("h4nobtt_dPhi3PU","h4nobtt_dPhi3PU",50,  0, 3.2)

    ####geq 5 jets
    nJets5 = r.TH1D("h_nJets5","h_nJets5",15, -0.5, 14.5)
    bJets5 = r.TH1D("h_bJets5","h_bJets5",5,  -0.5, 4.5)
    nJets5TopTagged = r.TH1D("h_nJets5TopTagged","h_nJets5TopTagged",15, -0.5, 14.5)
    bJets5TopTagged = r.TH1D("h_bJets5TopTagged","h_bJets5TopTagged",5,  -0.5, 4.5)

    nJets5PU = r.TH1D("h_nJets5PU","h_nJets5PU",15, -0.5, 14.5)
    bJets5PU = r.TH1D("h_bJets5PU","h_bJets5PU",5,  -0.5, 4.5)
    nJets5PUTopTagged = r.TH1D("h_nJets5PUTopTagged","h_nJets5PUTopTagged",15, -0.5, 14.5)
    bJets5PUTopTagged = r.TH1D("h_bJets5PUTopTagged","h_bJets5PUTopTagged",5,  -0.5, 4.5)

    nJetsNobReq5 = r.TH1D("h_nJetsNobReq5","h_nJetsNobReq5",15, -0.5, 14.5)
    bJetsNobReq5 = r.TH1D("h_bJetsNobReq5","h_bJetsNobReq5",5,  -0.5, 4.5)
    nJetsNobReq5TopTagged = r.TH1D("h_nJetsNobReq5TopTagged","h_nJetsNobReq5TopTagged",15, -0.5, 14.5)
    bJetsNobReq5TopTagged = r.TH1D("h_bJetsNobReq5TopTagged","h_bJetsNobReq5TopTagged",5,  -0.5, 4.5)

    nJetsNobReq5PU = r.TH1D("h_nJetsNobReq5PU","h_nJetsNobReq5PU",15, -0.5, 14.5)
    bJetsNobReq5PU = r.TH1D("h_bJetsNobReq5PU","h_bJetsNobReq5PU",5,  -0.5, 4.5)
    nJetsNobReq5PUTopTagged = r.TH1D("h_nJetsNobReq5PUTopTagged","h_nJetsNobReq5PUTopTagged",15, -0.5, 14.5)
    bJetsNobReq5PUTopTagged = r.TH1D("h_bJetsNobReq5PUTopTagged","h_bJetsNobReq5PUTopTagged",5,  -0.5, 4.5)

    ####geq 6 jets
    nJets6 = r.TH1D("h_nJets6","h_nJets6",15, -0.5, 14.5)
    bJets6 = r.TH1D("h_bJets6","h_bJets6",5,  -0.5, 4.5)
    nJets6TopTagged = r.TH1D("h_nJets6TopTagged","h_nJets6TopTagged",15, -0.5, 14.5)
    bJets6TopTagged = r.TH1D("h_bJets6TopTagged","h_bJets6TopTagged",5,  -0.5, 4.5)

    nJets6PU = r.TH1D("h_nJets6PU","h_nJets6PU",15, -0.5, 14.5)
    bJets6PU = r.TH1D("h_bJets6PU","h_bJets6PU",5,  -0.5, 4.5)
    nJets6PUTopTagged = r.TH1D("h_nJets6PUTopTagged","h_nJets6PUTopTagged",15, -0.5, 14.5)
    bJets6PUTopTagged = r.TH1D("h_bJets6PUTopTagged","h_bJets6PUTopTagged",5,  -0.5, 4.5)

    nJetsNobReq6 = r.TH1D("h_nJetsNobReq6","h_nJetsNobReq6",15, -0.5, 14.5)
    bJetsNobReq6 = r.TH1D("h_bJetsNobReq6","h_bJetsNobReq6",5,  -0.5, 4.5)
    nJetsNobReq6TopTagged = r.TH1D("h_nJetsNobReq6TopTagged","h_nJetsNobReq6TopTagged",15, -0.5, 14.5)
    bJetsNobReq6TopTagged = r.TH1D("h_bJetsNobReq6TopTagged","h_bJetsNobReq6TopTagged",5,  -0.5, 4.5)

    nJetsNobReq6PU = r.TH1D("h_nJetsNobReq6PU","h_nJetsNobReq6PU",15, -0.5, 14.5)
    bJetsNobReq6PU = r.TH1D("h_bJetsNobReq6PU","h_bJetsNobReq6PU",5,  -0.5, 4.5)
    nJetsNobReq6PUTopTagged = r.TH1D("h_nJetsNobReq6PUTopTagged","h_nJetsNobReq6PUTopTagged",15, -0.5, 14.5)
    bJetsNobReq6PUTopTagged = r.TH1D("h_bJetsNobReq6PUTopTagged","h_bJetsNobReq6PUTopTagged",5,  -0.5, 4.5)

    nVertices    .Sumw2()
    nVerticesReWt.Sumw2()
    nJetsAll.Sumw2()
    bJetsAll.Sumw2()
    jet1Pt.Sumw2()
    jet2Pt.Sumw2()
    met  .Sumw2()
    mht  .Sumw2()
    dPhi1.Sumw2()
    dPhi2.Sumw2()
    dPhi3.Sumw2()

    nJetsAllPU.Sumw2()
    bJetsAllPU.Sumw2()
    jet1PtPU.Sumw2()
    jet2PtPU.Sumw2()
    metPU  .Sumw2()
    mhtPU  .Sumw2()
    dPhi1PU.Sumw2()
    dPhi2PU.Sumw2()
    dPhi3PU.Sumw2()

    nJets4.Sumw2()
    bJets4.Sumw2()
    nJets4PU.Sumw2()
    bJets4PU.Sumw2()
    nJets4TopTagged.Sumw2()
    bJets4TopTagged.Sumw2()
    nJets4PUTopTagged.Sumw2()
    bJets4PUTopTagged.Sumw2()

    nJetsNobReq4.Sumw2()
    bJetsNobReq4.Sumw2()
    nJetsNobReq4PU.Sumw2()
    bJetsNobReq4PU.Sumw2()
    nJetsNobReq4TopTagged.Sumw2()
    bJetsNobReq4TopTagged.Sumw2()
    nJetsNobReq4PUTopTagged.Sumw2()
    bJetsNobReq4PUTopTagged.Sumw2()

    nJets5.Sumw2()
    bJets5.Sumw2()
    nJets5PU.Sumw2()
    bJets5PU.Sumw2()
    nJets5TopTagged.Sumw2()
    bJets5TopTagged.Sumw2()
    nJets5PUTopTagged.Sumw2()
    bJets5PUTopTagged.Sumw2()

    nJetsNobReq5.Sumw2()
    bJetsNobReq5.Sumw2()
    nJetsNobReq5PU.Sumw2()
    bJetsNobReq5PU.Sumw2()
    nJetsNobReq5TopTagged.Sumw2()
    bJetsNobReq5TopTagged.Sumw2()
    nJetsNobReq5PUTopTagged.Sumw2()
    bJetsNobReq5PUTopTagged.Sumw2()

    nJets6.Sumw2()
    bJets6.Sumw2()
    nJets6PU.Sumw2()
    bJets6PU.Sumw2()
    nJets6TopTagged.Sumw2()
    bJets6TopTagged.Sumw2()
    nJets6PUTopTagged.Sumw2()
    bJets6PUTopTagged.Sumw2()

    nJetsNobReq6.Sumw2()
    bJetsNobReq6.Sumw2()
    nJetsNobReq6PU.Sumw2()
    bJetsNobReq6PU.Sumw2()
    nJetsNobReq6TopTagged.Sumw2()
    bJetsNobReq6TopTagged.Sumw2()
    nJetsNobReq6PUTopTagged.Sumw2()
    bJetsNobReq6PUTopTagged.Sumw2()

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
        ra2_nJetsPt30Eta24  = event.ra2_nJetsPt30Eta24 
        ra2_bJetsPt30Eta24  = event.ra2_bJetsPt30Eta24 
        ra2_PUWt            = event.ra2_PUWt           
        ra2_EventWt         = event.ra2_EventWt        
        ra2_MET             = event.ra2_MET            
        ra2_MHT             = event.ra2_MHT            
        ra2_dPhi1           = event.ra2_dPhiMET1         
        ra2_dPhi2           = event.ra2_dPhiMET2         
        ra2_dPhi3           = event.ra2_dPhiMET3         
        ra2_Jet1Pt          = event.ra2_Jet1Pt         
        ra2_Jet1Eta         = event.ra2_Jet1Eta        
        ra2_Jet2Pt          = event.ra2_Jet2Pt         
        ra2_Jet2Eta         = event.ra2_Jet2Eta        

    ##top tagger variables 
        ra2_bestTopJetMass  = event.ra2_bestTopJetMass 
        ra2_TbestTopJet     = event.ra2_TbestTopJet    
        ra2_TbJet           = event.ra2_TbJet          
        ra2_MT2             = event.ra2_MT2            

        #Cuts
        #jet1Cuts = (ra2_Jet1Pt>70 and (abs(ra2_Jet1Eta)<2.4))
        #jet2Cuts = (ra2_Jet2Pt>70 and (abs(ra2_Jet2Eta)<2.4))
        #dPhiCuts = (ra2_dPhi1>0.5 and ra2_dPhi2>0.5 and ra2_dPhi3>0.3)
        #bJetCuts = (ra2_bJetsPt30Eta24>0)
        #base4cuts = ((ra2_MHT>175 and ra2_nJetsPt30Eta24>3) and jet1Cuts and jet2Cuts and dPhiCuts and bJetCuts)
        #base5cuts = ((ra2_MHT>175 and ra2_nJetsPt30Eta24>4) and jet1Cuts and jet2Cuts and dPhiCuts and bJetCuts)
        #topTaggerCuts = (ra2_bestTopJetMass>80 and ra2_bestTopJetMass<270 and (ra2_TbJet+0.5*ra2_TbestTopJet)>500 and ra2_MT2>300)

        ####Fill plots
        nVertices.Fill(ra2_Vertices,ra2_EventWt)
        nVerticesReWt.Fill(ra2_Vertices,ra2_EventWt*ra2_PUWt)
        nJetsAll.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
        bJetsAll.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
        jet1Pt.Fill(ra2_Jet1Pt,ra2_EventWt)
        jet2Pt.Fill(ra2_Jet2Pt,ra2_EventWt)
        met  .Fill(ra2_MET,ra2_EventWt)
        mht  .Fill(ra2_MHT,ra2_EventWt)
        dPhi1.Fill(ra2_dPhi1,ra2_EventWt)
        dPhi2.Fill(ra2_dPhi2,ra2_EventWt)
        dPhi3.Fill(ra2_dPhi3,ra2_EventWt)
    
        nJetsAllPU.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
        bJetsAllPU.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
        jet1PtPU.Fill(ra2_Jet1Pt,ra2_EventWt*ra2_PUWt)
        jet2PtPU.Fill(ra2_Jet2Pt,ra2_EventWt*ra2_PUWt)
        metPU  .Fill(ra2_MET,ra2_EventWt*ra2_PUWt)
        mhtPU  .Fill(ra2_MHT,ra2_EventWt*ra2_PUWt)
        dPhi1PU.Fill(ra2_dPhi1,ra2_EventWt*ra2_PUWt)
        dPhi2PU.Fill(ra2_dPhi2,ra2_EventWt*ra2_PUWt)
        dPhi3PU.Fill(ra2_dPhi3,ra2_EventWt*ra2_PUWt)
    

        if cutF.hadStopBaseline(event,1,4):
            nJets4.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
            bJets4.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
            nJets4PU.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            bJets4PU.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            if cutF.topTaggerCuts(event):
                nJets4TopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
                bJets4TopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                nJets4PUTopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                bJets4PUTopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                
        if cutF.hadStopBaseline(event,0,4):
            nJetsNobReq4.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
            bJetsNobReq4.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
            nJetsNobReq4PU.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            bJetsNobReq4PU.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            if cutF.topTaggerCuts(event):
                nJetsNobReq4TopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
                bJetsNobReq4TopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                nJetsNobReq4PUTopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                bJetsNobReq4PUTopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                
        if cutF.hadStopBaseline(event,1,5):
            nJets5.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
            bJets5.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
            nJets5PU.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            bJets5PU.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            if cutF.topTaggerCuts(event):
                nJets5TopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
                bJets5TopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                nJets5PUTopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                bJets5PUTopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                
        if cutF.hadStopBaseline(event,0,5):
            nJetsNobReq5.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
            bJetsNobReq5.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
            nJetsNobReq5PU.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            bJetsNobReq5PU.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            if cutF.topTaggerCuts(event):
                nJetsNobReq5TopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
                bJetsNobReq5TopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                nJetsNobReq5PUTopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                bJetsNobReq5PUTopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)

        if cutF.hadStopBaseline(event,1,6):
            nJets6.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
            bJets6.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
            nJets6PU.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            bJets6PU.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            if cutF.topTaggerCuts(event):
                nJets6TopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
                bJets6TopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                nJets6PUTopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                bJets6PUTopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                
        if cutF.hadStopBaseline(event,0,6):
            nJetsNobReq6.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
            bJetsNobReq6.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
            nJetsNobReq6PU.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            bJetsNobReq6PU.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            if cutF.topTaggerCuts(event):
                nJetsNobReq6TopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
                bJetsNobReq6TopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                nJetsNobReq6PUTopTagged.Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                bJetsNobReq6PUTopTagged.Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)

        #########
        i = i + 1

    #####        
    outputFile.cd()
    nVertices    .Write()
    nVerticesReWt.Write()
    nJetsAll.Write()
    bJetsAll.Write()
    jet1Pt.Write()
    jet2Pt.Write()
    met  .Write()
    mht  .Write()
    dPhi1.Write()
    dPhi2.Write()
    dPhi3.Write()
    nJetsAllPU.Write()
    bJetsAllPU.Write()
    jet1PtPU.Write()
    jet2PtPU.Write()
    metPU  .Write()
    mhtPU  .Write()
    dPhi1PU.Write()
    dPhi2PU.Write()
    dPhi3PU.Write()


    nJets4.Write()
    bJets4.Write()
    nJets4PU.Write()
    bJets4PU.Write()
    nJets4TopTagged.Write()
    bJets4TopTagged.Write()
    nJets4PUTopTagged.Write()
    bJets4PUTopTagged.Write()

    nJets5.Write()
    bJets5.Write()
    nJets5PU.Write()
    bJets5PU.Write()
    nJets5TopTagged.Write()
    bJets5TopTagged.Write()
    nJets5PUTopTagged.Write()
    bJets5PUTopTagged.Write()

    nJetsNobReq4.Write()
    bJetsNobReq4.Write()
    nJetsNobReq4PU.Write()
    bJetsNobReq4PU.Write()
    nJetsNobReq4TopTagged.Write()
    bJetsNobReq4TopTagged.Write()
    nJetsNobReq4PUTopTagged.Write()
    bJetsNobReq4PUTopTagged.Write()

    nJetsNobReq5.Write()
    bJetsNobReq5.Write()
    nJetsNobReq5PU.Write()
    bJetsNobReq5PU.Write()
    nJetsNobReq5TopTagged.Write()
    bJetsNobReq5TopTagged.Write()
    nJetsNobReq5PUTopTagged.Write()
    bJetsNobReq5PUTopTagged.Write()

    outputFile.Write()
    outputFile.Close()

##########################
if __name__ == '__main__':
    main()
    
