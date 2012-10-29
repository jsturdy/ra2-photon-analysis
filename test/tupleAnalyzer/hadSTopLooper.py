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
    parser.add_option('-s', action="store",      default="zinv",dest="sample", type="string")
    options, args = parser.parse_args()

    r.gROOT.SetBatch(True)
    #debug = False
    myWorkingDir = os.getcwd()
    
    outFileName = "topTaggingPlots_%s_min%d.root"%(options.sample,options.numJets)
    #outFileName = "topTaggingPlots_DY_ht_no_muon_removal.root"
    #outFileName = "topTaggingPlots_DY_ht_dimuons_removed.root"
    if options.debug:
        outFileName = "topTaggingPlots_%s_test_min%d.root"%(options.sample,options.numJets)

    print outFileName
    outputFile = r.TFile(outFileName,"RECREATE")

    myChain = r.TChain('analysis/RA2Values')
    #myChain = r.TChain('analysisWithMuon/RA2Values')
    #myChain = r.TChain('analysisNoMuon/RA2Values')

    if options.debug:
        if options.sample=="ttbar":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/ttjets_tunez2star_reco_tree_topTagged_v2/res/*_?_?_???.root")
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht200_reco_tree_topTagged_v2/res/*_?_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht400_reco_tree_topTagged_v2/res/*_?_?_???.root")
    else:
        if options.sample=="ttbar":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/ttjets_tunez2star_reco_tree_topTagged_v2/res/*.root")
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht200_reco_tree_topTagged_v2/res/*.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/treeMaker/zinvht400_reco_tree_topTagged_v2/res/*.root")
    
    ###variables
    nVertices     = r.TH1D("h_nVertices",    "h_nVertices",    50, -0.5, 49.5)
    nVerticesReWt = r.TH1D("h_nVerticesReWt","h_nVerticesReWt",50, -0.5, 49.5)

    nJets  = r.TH1D("h_nJets", "h_nJets", 15, -0.5, 14.5)
    bJets  = r.TH1D("h_bJets", "h_bJets", 10,  -0.5, 9.5)
    jet1Pt = r.TH1D("h_jet1Pt","h_jet1Pt",2*50,  0, 1000)
    jet2Pt = r.TH1D("h_jet2Pt","h_jet2Pt",2*50,  0, 500)
    met    = r.TH1D("h_met",   "h_met",   2*50,  0, 1000)
    mht    = r.TH1D("h_mht",   "h_mht",   2*50,  0, 1000)

    bestTopJetMass  = r.TH1D("h_bestTopJetMass", "h_bestTopJetMass", 2*50, 0, 1000)
    TbestTopJet     = r.TH1D("h_TbestTopJet"   , "h_TbestTopJet"   , 2*75, 0, 1500)
    TbJet           = r.TH1D("h_TbJet"         , "h_TbJet"         , 2*75, 0, 1500)
    mt2             = r.TH1D("h_mt2"           , "h_mt2"           , 2*50, 0, 1000)
    
    dPhi1  = r.TH1D("h_dPhi1","h_dPhi1",50,  0, 3.2)
    dPhi2  = r.TH1D("h_dPhi2","h_dPhi2",50,  0, 3.2)
    dPhi3  = r.TH1D("h_dPhi3","h_dPhi3",50,  0, 3.2)

    nJetsPU  = r.TH1D("h_nJetsPU", "h_nJetsPU", 15, -0.5, 14.5)
    bJetsPU  = r.TH1D("h_bJetsPU", "h_bJetsPU", 10,  -0.5, 9.5)
    jet1PtPU = r.TH1D("h_jet1PtPU","h_jet1PtPU",2*50,  0, 1000)
    jet2PtPU = r.TH1D("h_jet2PtPU","h_jet2PtPU",2*50,  0, 500)
    metPU    = r.TH1D("h_metPU",   "h_metPU",   2*50,  0, 1000)
    mhtPU    = r.TH1D("h_mhtPU",   "h_mhtPU",   2*50,  0, 1000)

    bestTopJetMassPU  = r.TH1D("h_bestTopJetMassPU", "h_bestTopJetMassPU", 2*50, 0, 1000)
    TbestTopJetPU     = r.TH1D("h_TbestTopJetPU"   , "h_TbestTopJetPU"   , 2*75, 0, 1500)
    TbJetPU           = r.TH1D("h_TbJetPU"         , "h_TbJetPU"         , 2*75, 0, 1500)
    mt2PU             = r.TH1D("h_mt2PU"           , "h_mt2PU"           , 2*50, 0, 1000)
    
    dPhi1PU  = r.TH1D("h_dPhi1PU","h_dPhi1PU",50,  0, 3.2)
    dPhi2PU  = r.TH1D("h_dPhi2PU","h_dPhi2PU",50,  0, 3.2)
    dPhi3PU  = r.TH1D("h_dPhi3PU","h_dPhi3PU",50,  0, 3.2)

    ####geq 4 jets
    hj_nJets  = r.TH1D("hj_nJets", "hj_nJets", 15, -0.5, 14.5)
    hj_bJets  = r.TH1D("hj_bJets", "hj_bJets", 10,  -0.5, 9.5)
    hj_jet1Pt = r.TH1D("hj_jet1Pt","hj_jet1Pt",2*50,  0, 1000)
    hj_jet2Pt = r.TH1D("hj_jet2Pt","hj_jet2Pt",2*50,  0, 500)
    hj_met    = r.TH1D("hj_met",   "hj_met",   2*50,  0, 1000)
    hj_mht    = r.TH1D("hj_mht",   "hj_mht",   2*50,  0, 1000)
    hj_bestTopJetMass  = r.TH1D("hj_bestTopJetMass", "hj_bestTopJetMass", 2*50, 0, 1000)
    hj_TbestTopJet     = r.TH1D("hj_TbestTopJet"   , "hj_TbestTopJet"   , 2*75, 0, 1500)
    hj_TbJet           = r.TH1D("hj_TbJet"         , "hj_TbJet"         , 2*75, 0, 1500)
    hj_mt2             = r.TH1D("hj_mt2"           , "hj_mt2"           , 2*50, 0, 1000)
    hj_dPhi1  = r.TH1D("hj_dPhi1","hj_dPhi1",50,  0, 3.2)
    hj_dPhi2  = r.TH1D("hj_dPhi2","hj_dPhi2",50,  0, 3.2)
    hj_dPhi3  = r.TH1D("hj_dPhi3","hj_dPhi3",50,  0, 3.2)

    hjtt_nJets  = r.TH1D("hjtt_nJets", "hjtt_nJets", 15, -0.5, 14.5)
    hjtt_bJets  = r.TH1D("hjtt_bJets", "hjtt_bJets", 10,  -0.5, 9.5)
    hjtt_jet1Pt = r.TH1D("hjtt_jet1Pt","hjtt_jet1Pt",2*50,  0, 1000)
    hjtt_jet2Pt = r.TH1D("hjtt_jet2Pt","hjtt_jet2Pt",2*50,  0, 500)
    hjtt_met    = r.TH1D("hjtt_met",   "hjtt_met",   2*50,  0, 1000)
    hjtt_mht    = r.TH1D("hjtt_mht",   "hjtt_mht",   2*50,  0, 1000)
    hjtt_bestTopJetMass  = r.TH1D("hjtt_bestTopJetMass", "hjtt_bestTopJetMass", 2*50, 0, 1000)
    hjtt_TbestTopJet     = r.TH1D("hjtt_TbestTopJet"   , "hjtt_TbestTopJet"   , 2*75, 0, 1500)
    hjtt_TbJet           = r.TH1D("hjtt_TbJet"         , "hjtt_TbJet"         , 2*75, 0, 1500)
    hjtt_mt2             = r.TH1D("hjtt_mt2"           , "hjtt_mt2"           , 2*50, 0, 1000)
    hjtt_dPhi1  = r.TH1D("hjtt_dPhi1","hjtt_dPhi1",50,  0, 3.2)
    hjtt_dPhi2  = r.TH1D("hjtt_dPhi2","hjtt_dPhi2",50,  0, 3.2)
    hjtt_dPhi3  = r.TH1D("hjtt_dPhi3","hjtt_dPhi3",50,  0, 3.2)

    hj_nJetsPU  = r.TH1D("hj_nJetsPU", "hj_nJetsPU", 15, -0.5, 14.5)
    hj_bJetsPU  = r.TH1D("hj_bJetsPU", "hj_bJetsPU", 10,  -0.5, 9.5)
    hj_jet1PtPU = r.TH1D("hj_jet1PtPU","hj_jet1PtPU",2*50,  0, 1000)
    hj_jet2PtPU = r.TH1D("hj_jet2PtPU","hj_jet2PtPU",2*50,  0, 500)
    hj_metPU    = r.TH1D("hj_metPU",   "hj_metPU",   2*50,  0, 1000)
    hj_mhtPU    = r.TH1D("hj_mhtPU",   "hj_mhtPU",   2*50,  0, 1000)
    hj_bestTopJetMassPU  = r.TH1D("hj_bestTopJetMassPU", "hj_bestTopJetMassPU", 2*50, 0, 1000)
    hj_TbestTopJetPU     = r.TH1D("hj_TbestTopJetPU"   , "hj_TbestTopJetPU"   , 2*75, 0, 1500)
    hj_TbJetPU           = r.TH1D("hj_TbJetPU"         , "hj_TbJetPU"         , 2*75, 0, 1500)
    hj_mt2PU             = r.TH1D("hj_mt2PU"           , "hj_mt2PU"           , 2*50, 0, 1000)
    hj_dPhi1PU  = r.TH1D("hj_dPhi1PU","hj_dPhi1PU",50,  0, 3.2)
    hj_dPhi2PU  = r.TH1D("hj_dPhi2PU","hj_dPhi2PU",50,  0, 3.2)
    hj_dPhi3PU  = r.TH1D("hj_dPhi3PU","hj_dPhi3PU",50,  0, 3.2)

    hjtt_nJetsPU  = r.TH1D("hjtt_nJetsPU", "hjtt_nJetsPU", 15, -0.5, 14.5)
    hjtt_bJetsPU  = r.TH1D("hjtt_bJetsPU", "hjtt_bJetsPU", 10,  -0.5, 9.5)
    hjtt_jet1PtPU = r.TH1D("hjtt_jet1PtPU","hjtt_jet1PtPU",2*50,  0, 1000)
    hjtt_jet2PtPU = r.TH1D("hjtt_jet2PtPU","hjtt_jet2PtPU",2*50,  0, 500)
    hjtt_metPU    = r.TH1D("hjtt_metPU",   "hjtt_metPU",   2*50,  0, 1000)
    hjtt_mhtPU    = r.TH1D("hjtt_mhtPU",   "hjtt_mhtPU",   2*50,  0, 1000)
    hjtt_bestTopJetMassPU  = r.TH1D("hjtt_bestTopJetMassPU", "hjtt_bestTopJetMassPU", 2*50, 0, 1000)
    hjtt_TbestTopJetPU     = r.TH1D("hjtt_TbestTopJetPU"   , "hjtt_TbestTopJetPU"   , 2*75, 0, 1500)
    hjtt_TbJetPU           = r.TH1D("hjtt_TbJetPU"         , "hjtt_TbJetPU"         , 2*75, 0, 1500)
    hjtt_mt2PU             = r.TH1D("hjtt_mt2PU"           , "hjtt_mt2PU"           , 2*50, 0, 1000)
    hjtt_dPhi1PU  = r.TH1D("hjtt_dPhi1PU","hjtt_dPhi1PU",50,  0, 3.2)
    hjtt_dPhi2PU  = r.TH1D("hjtt_dPhi2PU","hjtt_dPhi2PU",50,  0, 3.2)
    hjtt_dPhi3PU  = r.TH1D("hjtt_dPhi3PU","hjtt_dPhi3PU",50,  0, 3.2)

    ##no requirement of a b-jet
    hjnob_nJets  = r.TH1D("hjnob_nJets", "hjnob_nJets", 15, -0.5, 14.5)
    hjnob_bJets  = r.TH1D("hjnob_bJets", "hjnob_bJets", 10,  -0.5, 9.5)
    hjnob_jet1Pt = r.TH1D("hjnob_jet1Pt","hjnob_jet1Pt",2*50,  0, 1000)
    hjnob_jet2Pt = r.TH1D("hjnob_jet2Pt","hjnob_jet2Pt",2*50,  0, 500)
    hjnob_met    = r.TH1D("hjnob_met",   "hjnob_met",   2*50,  0, 1000)
    hjnob_mht    = r.TH1D("hjnob_mht",   "hjnob_mht",   2*50,  0, 1000)
    hjnob_bestTopJetMass  = r.TH1D("hjnob_bestTopJetMass", "hjnob_bestTopJetMass", 2*50, 0, 1000)
    hjnob_TbestTopJet     = r.TH1D("hjnob_TbestTopJet"   , "hjnob_TbestTopJet"   , 2*75, 0, 1500)
    hjnob_TbJet           = r.TH1D("hjnob_TbJet"         , "hjnob_TbJet"         , 2*75, 0, 1500)
    hjnob_mt2             = r.TH1D("hjnob_mt2"           , "hjnob_mt2"           , 2*50, 0, 1000)
    hjnob_dPhi1  = r.TH1D("hjnob_dPhi1","hjnob_dPhi1",50,  0, 3.2)
    hjnob_dPhi2  = r.TH1D("hjnob_dPhi2","hjnob_dPhi2",50,  0, 3.2)
    hjnob_dPhi3  = r.TH1D("hjnob_dPhi3","hjnob_dPhi3",50,  0, 3.2)

    hjnobtt_nJets  = r.TH1D("hjnobtt_nJets", "hjnobtt_nJets", 15, -0.5, 14.5)
    hjnobtt_bJets  = r.TH1D("hjnobtt_bJets", "hjnobtt_bJets", 10,  -0.5, 9.5)
    hjnobtt_jet1Pt = r.TH1D("hjnobtt_jet1Pt","hjnobtt_jet1Pt",2*50,  0, 1000)
    hjnobtt_jet2Pt = r.TH1D("hjnobtt_jet2Pt","hjnobtt_jet2Pt",2*50,  0, 500)
    hjnobtt_met    = r.TH1D("hjnobtt_met",   "hjnobtt_met",   2*50,  0, 1000)
    hjnobtt_mht    = r.TH1D("hjnobtt_mht",   "hjnobtt_mht",   2*50,  0, 1000)
    hjnobtt_bestTopJetMass  = r.TH1D("hjnobtt_bestTopJetMass", "hjnobtt_bestTopJetMass", 2*50, 0, 1000)
    hjnobtt_TbestTopJet     = r.TH1D("hjnobtt_TbestTopJet"   , "hjnobtt_TbestTopJet"   , 2*75, 0, 1500)
    hjnobtt_TbJet           = r.TH1D("hjnobtt_TbJet"         , "hjnobtt_TbJet"         , 2*75, 0, 1500)
    hjnobtt_mt2             = r.TH1D("hjnobtt_mt2"           , "hjnobtt_mt2"           , 2*50, 0, 1000)
    hjnobtt_dPhi1  = r.TH1D("hjnobtt_dPhi1","hjnobtt_dPhi1",50,  0, 3.2)
    hjnobtt_dPhi2  = r.TH1D("hjnobtt_dPhi2","hjnobtt_dPhi2",50,  0, 3.2)
    hjnobtt_dPhi3  = r.TH1D("hjnobtt_dPhi3","hjnobtt_dPhi3",50,  0, 3.2)

    ##no requirement of a b-jet PU re-weighted
    hjnob_nJetsPU  = r.TH1D("hjnob_nJetsPU", "hjnob_nJetsPU", 15, -0.5, 14.5)
    hjnob_bJetsPU  = r.TH1D("hjnob_bJetsPU", "hjnob_bJetsPU", 10,  -0.5, 9.5)
    hjnob_jet1PtPU = r.TH1D("hjnob_jet1PtPU","hjnob_jet1PtPU",2*50,  0, 1000)
    hjnob_jet2PtPU = r.TH1D("hjnob_jet2PtPU","hjnob_jet2PtPU",2*50,  0, 500)
    hjnob_metPU    = r.TH1D("hjnob_metPU",   "hjnob_metPU",   2*50,  0, 1000)
    hjnob_mhtPU    = r.TH1D("hjnob_mhtPU",   "hjnob_mhtPU",   2*50,  0, 1000)
    hjnob_bestTopJetMassPU  = r.TH1D("hjnob_bestTopJetMassPU", "hjnob_bestTopJetMassPU", 2*50, 0, 1000)
    hjnob_TbestTopJetPU     = r.TH1D("hjnob_TbestTopJetPU"   , "hjnob_TbestTopJetPU"   , 2*75, 0, 1500)
    hjnob_TbJetPU           = r.TH1D("hjnob_TbJetPU"         , "hjnob_TbJetPU"         , 2*75, 0, 1500)
    hjnob_mt2PU             = r.TH1D("hjnob_mt2PU"           , "hjnob_mt2PU"           , 2*50, 0, 1000)
    hjnob_dPhi1PU  = r.TH1D("hjnob_dPhi1PU","hjnob_dPhi1PU",50,  0, 3.2)
    hjnob_dPhi2PU  = r.TH1D("hjnob_dPhi2PU","hjnob_dPhi2PU",50,  0, 3.2)
    hjnob_dPhi3PU  = r.TH1D("hjnob_dPhi3PU","hjnob_dPhi3PU",50,  0, 3.2)

    hjnobtt_nJetsPU  = r.TH1D("hjnobtt_nJetsPU", "hjnobtt_nJetsPU", 15, -0.5, 14.5)
    hjnobtt_bJetsPU  = r.TH1D("hjnobtt_bJetsPU", "hjnobtt_bJetsPU", 10,  -0.5, 9.5)
    hjnobtt_jet1PtPU = r.TH1D("hjnobtt_jet1PtPU","hjnobtt_jet1PtPU",2*50,  0, 1000)
    hjnobtt_jet2PtPU = r.TH1D("hjnobtt_jet2PtPU","hjnobtt_jet2PtPU",2*50,  0, 500)
    hjnobtt_metPU    = r.TH1D("hjnobtt_metPU",   "hjnobtt_metPU",   2*50,  0, 1000)
    hjnobtt_mhtPU    = r.TH1D("hjnobtt_mhtPU",   "hjnobtt_mhtPU",   2*50,  0, 1000)
    hjnobtt_bestTopJetMassPU  = r.TH1D("hjnobtt_bestTopJetMassPU", "hjnobtt_bestTopJetMassPU", 2*50, 0, 1000)
    hjnobtt_TbestTopJetPU     = r.TH1D("hjnobtt_TbestTopJetPU"   , "hjnobtt_TbestTopJetPU"   , 2*75, 0, 1500)
    hjnobtt_TbJetPU           = r.TH1D("hjnobtt_TbJetPU"         , "hjnobtt_TbJetPU"         , 2*75, 0, 1500)
    hjnobtt_mt2PU             = r.TH1D("hjnobtt_mt2PU"           , "hjnobtt_mt2PU"           , 2*50, 0, 1000)
    hjnobtt_dPhi1PU  = r.TH1D("hjnobtt_dPhi1PU","hjnobtt_dPhi1PU",50,  0, 3.2)
    hjnobtt_dPhi2PU  = r.TH1D("hjnobtt_dPhi2PU","hjnobtt_dPhi2PU",50,  0, 3.2)
    hjnobtt_dPhi3PU  = r.TH1D("hjnobtt_dPhi3PU","hjnobtt_dPhi3PU",50,  0, 3.2)

    ####
    nVertices    .Sumw2()
    nVerticesReWt.Sumw2()
    nJets.Sumw2()
    bJets.Sumw2()
    jet1Pt.Sumw2()
    jet2Pt.Sumw2()
    met  .Sumw2()
    mht  .Sumw2()
    bestTopJetMass.Sumw2()
    TbestTopJet   .Sumw2()
    TbJet         .Sumw2()
    mt2           .Sumw2()
    dPhi1.Sumw2()
    dPhi2.Sumw2()
    dPhi3.Sumw2()

    nJetsPU.Sumw2()
    bJetsPU.Sumw2()
    jet1PtPU.Sumw2()
    jet2PtPU.Sumw2()
    metPU  .Sumw2()
    mhtPU  .Sumw2()
    bestTopJetMassPU.Sumw2()
    TbestTopJetPU   .Sumw2()
    TbJetPU         .Sumw2()
    mt2PU           .Sumw2()
    dPhi1PU.Sumw2()
    dPhi2PU.Sumw2()
    dPhi3PU.Sumw2()

    ###
    hj_nJets.Sumw2()
    hj_bJets.Sumw2()
    hj_jet1Pt.Sumw2()
    hj_jet2Pt.Sumw2()
    hj_met  .Sumw2()
    hj_mht  .Sumw2()
    hj_bestTopJetMass.Sumw2()
    hj_TbestTopJet   .Sumw2()
    hj_TbJet         .Sumw2()
    hj_mt2           .Sumw2()
    hj_dPhi1.Sumw2()
    hj_dPhi2.Sumw2()
    hj_dPhi3.Sumw2()
    
    hj_nJetsPU.Sumw2()
    hj_bJetsPU.Sumw2()
    hj_jet1PtPU.Sumw2()
    hj_jet2PtPU.Sumw2()
    hj_metPU  .Sumw2()
    hj_mhtPU  .Sumw2()
    hj_bestTopJetMassPU.Sumw2()
    hj_TbestTopJetPU   .Sumw2()
    hj_TbJetPU         .Sumw2()
    hj_mt2PU           .Sumw2()
    hj_dPhi1PU.Sumw2()
    hj_dPhi2PU.Sumw2()
    hj_dPhi3PU.Sumw2()

    ###
    hjtt_nJets.Sumw2()
    hjtt_bJets.Sumw2()
    hjtt_jet1Pt.Sumw2()
    hjtt_jet2Pt.Sumw2()
    hjtt_met  .Sumw2()
    hjtt_mht  .Sumw2()
    hjtt_bestTopJetMass.Sumw2()
    hjtt_TbestTopJet   .Sumw2()
    hjtt_TbJet         .Sumw2()
    hjtt_mt2           .Sumw2()
    hjtt_dPhi1.Sumw2()
    hjtt_dPhi2.Sumw2()
    hjtt_dPhi3.Sumw2()
    
    hjtt_nJetsPU.Sumw2()
    hjtt_bJetsPU.Sumw2()
    hjtt_jet1PtPU.Sumw2()
    hjtt_jet2PtPU.Sumw2()
    hjtt_metPU  .Sumw2()
    hjtt_mhtPU  .Sumw2()
    hjtt_bestTopJetMassPU.Sumw2()
    hjtt_TbestTopJetPU   .Sumw2()
    hjtt_TbJetPU         .Sumw2()
    hjtt_mt2PU           .Sumw2()
    hjtt_dPhi1PU.Sumw2()
    hjtt_dPhi2PU.Sumw2()
    hjtt_dPhi3PU.Sumw2()

    ###
    hjnob_nJets.Sumw2()
    hjnob_bJets.Sumw2()
    hjnob_jet1Pt.Sumw2()
    hjnob_jet2Pt.Sumw2()
    hjnob_met  .Sumw2()
    hjnob_mht  .Sumw2()
    hjnob_bestTopJetMass.Sumw2()
    hjnob_TbestTopJet   .Sumw2()
    hjnob_TbJet         .Sumw2()
    hjnob_mt2           .Sumw2()
    hjnob_dPhi1.Sumw2()
    hjnob_dPhi2.Sumw2()
    hjnob_dPhi3.Sumw2()
    
    hjnob_nJetsPU.Sumw2()
    hjnob_bJetsPU.Sumw2()
    hjnob_jet1PtPU.Sumw2()
    hjnob_jet2PtPU.Sumw2()
    hjnob_metPU  .Sumw2()
    hjnob_mhtPU  .Sumw2()
    hjnob_bestTopJetMassPU.Sumw2()
    hjnob_TbestTopJetPU   .Sumw2()
    hjnob_TbJetPU         .Sumw2()
    hjnob_mt2PU           .Sumw2()
    hjnob_dPhi1PU.Sumw2()
    hjnob_dPhi2PU.Sumw2()
    hjnob_dPhi3PU.Sumw2()

    ###
    hjnobtt_nJets.Sumw2()
    hjnobtt_bJets.Sumw2()
    hjnobtt_jet1Pt.Sumw2()
    hjnobtt_jet2Pt.Sumw2()
    hjnobtt_met  .Sumw2()
    hjnobtt_mht  .Sumw2()
    hjnobtt_bestTopJetMass.Sumw2()
    hjnobtt_TbestTopJet   .Sumw2()
    hjnobtt_TbJet         .Sumw2()
    hjnobtt_mt2           .Sumw2()
    hjnobtt_dPhi1.Sumw2()
    hjnobtt_dPhi2.Sumw2()
    hjnobtt_dPhi3.Sumw2()
    
    hjnobtt_nJetsPU.Sumw2()
    hjnobtt_bJetsPU.Sumw2()
    hjnobtt_jet1PtPU.Sumw2()
    hjnobtt_jet2PtPU.Sumw2()
    hjnobtt_metPU  .Sumw2()
    hjnobtt_mhtPU  .Sumw2()
    hjnobtt_bestTopJetMassPU.Sumw2()
    hjnobtt_TbestTopJetPU   .Sumw2()
    hjnobtt_TbJetPU         .Sumw2()
    hjnobtt_mt2PU           .Sumw2()
    hjnobtt_dPhi1PU.Sumw2()
    hjnobtt_dPhi2PU.Sumw2()
    hjnobtt_dPhi3PU.Sumw2()


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

        ####Fill plots
        nVertices.Fill(ra2_Vertices,ra2_EventWt)
        nVerticesReWt.Fill(ra2_Vertices,ra2_EventWt*ra2_PUWt)
    
        nJets      .Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
        bJets      .Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
        jet1Pt        .Fill(ra2_Jet1Pt        ,ra2_EventWt)
        jet2Pt        .Fill(ra2_Jet2Pt        ,ra2_EventWt)
        met           .Fill(ra2_MET           ,ra2_EventWt)
        mht           .Fill(ra2_MHT           ,ra2_EventWt)
        bestTopJetMass.Fill(ra2_bestTopJetMass,ra2_EventWt)
        TbestTopJet   .Fill(ra2_TbestTopJet   ,ra2_EventWt)
        TbJet         .Fill(ra2_TbJet         ,ra2_EventWt)
        mt2           .Fill(ra2_MT2           ,ra2_EventWt)
        dPhi1         .Fill(ra2_dPhi1         ,ra2_EventWt)
        dPhi2         .Fill(ra2_dPhi2         ,ra2_EventWt)
        dPhi3         .Fill(ra2_dPhi3         ,ra2_EventWt)
        
        nJetsPU      .Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
        bJetsPU      .Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
        jet1PtPU        .Fill(ra2_Jet1Pt        ,ra2_EventWt*ra2_PUWt)
        jet2PtPU        .Fill(ra2_Jet2Pt        ,ra2_EventWt*ra2_PUWt)
        metPU           .Fill(ra2_MET           ,ra2_EventWt*ra2_PUWt)
        mhtPU           .Fill(ra2_MHT           ,ra2_EventWt*ra2_PUWt)
        bestTopJetMassPU.Fill(ra2_bestTopJetMass,ra2_EventWt*ra2_PUWt)
        TbestTopJetPU   .Fill(ra2_TbestTopJet   ,ra2_EventWt*ra2_PUWt)
        TbJetPU         .Fill(ra2_TbJet         ,ra2_EventWt*ra2_PUWt)
        mt2PU           .Fill(ra2_MT2           ,ra2_EventWt*ra2_PUWt)
        dPhi1PU         .Fill(ra2_dPhi1         ,ra2_EventWt*ra2_PUWt)
        dPhi2PU         .Fill(ra2_dPhi2         ,ra2_EventWt*ra2_PUWt)
        dPhi3PU         .Fill(ra2_dPhi3         ,ra2_EventWt*ra2_PUWt)
        
        ### at least one b-jet
        if cutF.hadStopBaseline(event,1,options.numJets):
            hj_nJets      .Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
            hj_bJets      .Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
            hj_jet1Pt        .Fill(ra2_Jet1Pt        ,ra2_EventWt)
            hj_jet2Pt        .Fill(ra2_Jet2Pt        ,ra2_EventWt)
            hj_met           .Fill(ra2_MET           ,ra2_EventWt)
            hj_mht           .Fill(ra2_MHT           ,ra2_EventWt)
            hj_bestTopJetMass.Fill(ra2_bestTopJetMass,ra2_EventWt)
            hj_TbestTopJet   .Fill(ra2_TbestTopJet   ,ra2_EventWt)
            hj_TbJet         .Fill(ra2_TbJet         ,ra2_EventWt)
            hj_mt2           .Fill(ra2_MT2           ,ra2_EventWt)
            hj_dPhi1         .Fill(ra2_dPhi1         ,ra2_EventWt)
            hj_dPhi2         .Fill(ra2_dPhi2         ,ra2_EventWt)
            hj_dPhi3         .Fill(ra2_dPhi3         ,ra2_EventWt)
            
            hj_nJetsPU      .Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            hj_bJetsPU      .Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            hj_jet1PtPU        .Fill(ra2_Jet1Pt        ,ra2_EventWt*ra2_PUWt)
            hj_jet2PtPU        .Fill(ra2_Jet2Pt        ,ra2_EventWt*ra2_PUWt)
            hj_metPU           .Fill(ra2_MET           ,ra2_EventWt*ra2_PUWt)
            hj_mhtPU           .Fill(ra2_MHT           ,ra2_EventWt*ra2_PUWt)
            hj_bestTopJetMassPU.Fill(ra2_bestTopJetMass,ra2_EventWt*ra2_PUWt)
            hj_TbestTopJetPU   .Fill(ra2_TbestTopJet   ,ra2_EventWt*ra2_PUWt)
            hj_TbJetPU         .Fill(ra2_TbJet         ,ra2_EventWt*ra2_PUWt)
            hj_mt2PU           .Fill(ra2_MT2           ,ra2_EventWt*ra2_PUWt)
            hj_dPhi1PU         .Fill(ra2_dPhi1         ,ra2_EventWt*ra2_PUWt)
            hj_dPhi2PU         .Fill(ra2_dPhi2         ,ra2_EventWt*ra2_PUWt)
            hj_dPhi3PU         .Fill(ra2_dPhi3         ,ra2_EventWt*ra2_PUWt)

            if cutF.topTaggerCuts(event):
                hjtt_nJets      .Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
                hjtt_bJets      .Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                hjtt_jet1Pt        .Fill(ra2_Jet1Pt        ,ra2_EventWt)
                hjtt_jet2Pt        .Fill(ra2_Jet2Pt        ,ra2_EventWt)
                hjtt_met           .Fill(ra2_MET           ,ra2_EventWt)
                hjtt_mht           .Fill(ra2_MHT           ,ra2_EventWt)
                hjtt_bestTopJetMass.Fill(ra2_bestTopJetMass,ra2_EventWt)
                hjtt_TbestTopJet   .Fill(ra2_TbestTopJet   ,ra2_EventWt)
                hjtt_TbJet         .Fill(ra2_TbJet         ,ra2_EventWt)
                hjtt_mt2           .Fill(ra2_MT2           ,ra2_EventWt)
                hjtt_dPhi1         .Fill(ra2_dPhi1         ,ra2_EventWt)
                hjtt_dPhi2         .Fill(ra2_dPhi2         ,ra2_EventWt)
                hjtt_dPhi3         .Fill(ra2_dPhi3         ,ra2_EventWt)
                
                hjtt_nJetsPU      .Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                hjtt_bJetsPU      .Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                hjtt_jet1PtPU        .Fill(ra2_Jet1Pt        ,ra2_EventWt*ra2_PUWt)
                hjtt_jet2PtPU        .Fill(ra2_Jet2Pt        ,ra2_EventWt*ra2_PUWt)
                hjtt_metPU           .Fill(ra2_MET           ,ra2_EventWt*ra2_PUWt)
                hjtt_mhtPU           .Fill(ra2_MHT           ,ra2_EventWt*ra2_PUWt)
                hjtt_bestTopJetMassPU.Fill(ra2_bestTopJetMass,ra2_EventWt*ra2_PUWt)
                hjtt_TbestTopJetPU   .Fill(ra2_TbestTopJet   ,ra2_EventWt*ra2_PUWt)
                hjtt_TbJetPU         .Fill(ra2_TbJet         ,ra2_EventWt*ra2_PUWt)
                hjtt_mt2PU           .Fill(ra2_MT2           ,ra2_EventWt*ra2_PUWt)
                hjtt_dPhi1PU         .Fill(ra2_dPhi1         ,ra2_EventWt*ra2_PUWt)
                hjtt_dPhi2PU         .Fill(ra2_dPhi2         ,ra2_EventWt*ra2_PUWt)
                hjtt_dPhi3PU         .Fill(ra2_dPhi3         ,ra2_EventWt*ra2_PUWt)
                
        ### no b-jet requirement
        if cutF.hadStopBaseline(event,0,options.numJets):

            hjnob_nJets      .Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
            hjnob_bJets      .Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
            hjnob_jet1Pt        .Fill(ra2_Jet1Pt        ,ra2_EventWt)
            hjnob_jet2Pt        .Fill(ra2_Jet2Pt        ,ra2_EventWt)
            hjnob_met           .Fill(ra2_MET           ,ra2_EventWt)
            hjnob_mht           .Fill(ra2_MHT           ,ra2_EventWt)
            hjnob_bestTopJetMass.Fill(ra2_bestTopJetMass,ra2_EventWt)
            hjnob_TbestTopJet   .Fill(ra2_TbestTopJet   ,ra2_EventWt)
            hjnob_TbJet         .Fill(ra2_TbJet         ,ra2_EventWt)
            hjnob_mt2           .Fill(ra2_MT2           ,ra2_EventWt)
            hjnob_dPhi1         .Fill(ra2_dPhi1         ,ra2_EventWt)
            hjnob_dPhi2         .Fill(ra2_dPhi2         ,ra2_EventWt)
            hjnob_dPhi3         .Fill(ra2_dPhi3         ,ra2_EventWt)
            
            hjnob_nJetsPU      .Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            hjnob_bJetsPU      .Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
            hjnob_jet1PtPU        .Fill(ra2_Jet1Pt        ,ra2_EventWt*ra2_PUWt)
            hjnob_jet2PtPU        .Fill(ra2_Jet2Pt        ,ra2_EventWt*ra2_PUWt)
            hjnob_metPU           .Fill(ra2_MET           ,ra2_EventWt*ra2_PUWt)
            hjnob_mhtPU           .Fill(ra2_MHT           ,ra2_EventWt*ra2_PUWt)
            hjnob_bestTopJetMassPU.Fill(ra2_bestTopJetMass,ra2_EventWt*ra2_PUWt)
            hjnob_TbestTopJetPU   .Fill(ra2_TbestTopJet   ,ra2_EventWt*ra2_PUWt)
            hjnob_TbJetPU         .Fill(ra2_TbJet         ,ra2_EventWt*ra2_PUWt)
            hjnob_mt2PU           .Fill(ra2_MT2           ,ra2_EventWt*ra2_PUWt)
            hjnob_dPhi1PU         .Fill(ra2_dPhi1         ,ra2_EventWt*ra2_PUWt)
            hjnob_dPhi2PU         .Fill(ra2_dPhi2         ,ra2_EventWt*ra2_PUWt)
            hjnob_dPhi3PU         .Fill(ra2_dPhi3         ,ra2_EventWt*ra2_PUWt)

            if cutF.topTaggerCuts(event):
                hjnobtt_nJets      .Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
                hjnobtt_bJets      .Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                hjnobtt_jet1Pt        .Fill(ra2_Jet1Pt        ,ra2_EventWt)
                hjnobtt_jet2Pt        .Fill(ra2_Jet2Pt        ,ra2_EventWt)
                hjnobtt_met           .Fill(ra2_MET           ,ra2_EventWt)
                hjnobtt_mht           .Fill(ra2_MHT           ,ra2_EventWt)
                hjnobtt_bestTopJetMass.Fill(ra2_bestTopJetMass,ra2_EventWt)
                hjnobtt_TbestTopJet   .Fill(ra2_TbestTopJet   ,ra2_EventWt)
                hjnobtt_TbJet         .Fill(ra2_TbJet         ,ra2_EventWt)
                hjnobtt_mt2           .Fill(ra2_MT2           ,ra2_EventWt)
                hjnobtt_dPhi1         .Fill(ra2_dPhi1         ,ra2_EventWt)
                hjnobtt_dPhi2         .Fill(ra2_dPhi2         ,ra2_EventWt)
                hjnobtt_dPhi3         .Fill(ra2_dPhi3         ,ra2_EventWt)
                
                hjnobtt_nJetsPU      .Fill(ra2_nJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                hjnobtt_bJetsPU      .Fill(ra2_bJetsPt30Eta24,ra2_EventWt*ra2_PUWt)
                hjnobtt_jet1PtPU        .Fill(ra2_Jet1Pt        ,ra2_EventWt*ra2_PUWt)
                hjnobtt_jet2PtPU        .Fill(ra2_Jet2Pt        ,ra2_EventWt*ra2_PUWt)
                hjnobtt_metPU           .Fill(ra2_MET           ,ra2_EventWt*ra2_PUWt)
                hjnobtt_mhtPU           .Fill(ra2_MHT           ,ra2_EventWt*ra2_PUWt)
                hjnobtt_bestTopJetMassPU.Fill(ra2_bestTopJetMass,ra2_EventWt*ra2_PUWt)
                hjnobtt_TbestTopJetPU   .Fill(ra2_TbestTopJet   ,ra2_EventWt*ra2_PUWt)
                hjnobtt_TbJetPU         .Fill(ra2_TbJet         ,ra2_EventWt*ra2_PUWt)
                hjnobtt_mt2PU           .Fill(ra2_MT2           ,ra2_EventWt*ra2_PUWt)
                hjnobtt_dPhi1PU         .Fill(ra2_dPhi1         ,ra2_EventWt*ra2_PUWt)
                hjnobtt_dPhi2PU         .Fill(ra2_dPhi2         ,ra2_EventWt*ra2_PUWt)
                hjnobtt_dPhi3PU         .Fill(ra2_dPhi3         ,ra2_EventWt*ra2_PUWt)
                
        #########
        i = i + 1

    #####        
    outputFile.cd()

    ####
    nVertices    .Write()
    nVerticesReWt.Write()
    nJets.Write()
    bJets.Write()
    jet1Pt.Write()
    jet2Pt.Write()
    met  .Write()
    mht  .Write()
    bestTopJetMass.Write()
    TbestTopJet   .Write()
    TbJet         .Write()
    mt2           .Write()
    dPhi1.Write()
    dPhi2.Write()
    dPhi3.Write()

    nJetsPU.Write()
    bJetsPU.Write()
    jet1PtPU.Write()
    jet2PtPU.Write()
    metPU  .Write()
    mhtPU  .Write()
    bestTopJetMassPU.Write()
    TbestTopJetPU   .Write()
    TbJetPU         .Write()
    mt2PU           .Write()
    dPhi1PU.Write()
    dPhi2PU.Write()
    dPhi3PU.Write()

    ###
    hj_nJets.Write()
    hj_bJets.Write()
    hj_jet1Pt.Write()
    hj_jet2Pt.Write()
    hj_met  .Write()
    hj_mht  .Write()
    hj_bestTopJetMass.Write()
    hj_TbestTopJet   .Write()
    hj_TbJet         .Write()
    hj_mt2           .Write()
    hj_dPhi1.Write()
    hj_dPhi2.Write()
    hj_dPhi3.Write()
    
    hj_nJetsPU.Write()
    hj_bJetsPU.Write()
    hj_jet1PtPU.Write()
    hj_jet2PtPU.Write()
    hj_metPU  .Write()
    hj_mhtPU  .Write()
    hj_bestTopJetMassPU.Write()
    hj_TbestTopJetPU   .Write()
    hj_TbJetPU         .Write()
    hj_mt2PU           .Write()
    hj_dPhi1PU.Write()
    hj_dPhi2PU.Write()
    hj_dPhi3PU.Write()

    ###
    hjtt_nJets.Write()
    hjtt_bJets.Write()
    hjtt_jet1Pt.Write()
    hjtt_jet2Pt.Write()
    hjtt_met  .Write()
    hjtt_mht  .Write()
    hjtt_bestTopJetMass.Write()
    hjtt_TbestTopJet   .Write()
    hjtt_TbJet         .Write()
    hjtt_mt2           .Write()
    hjtt_dPhi1.Write()
    hjtt_dPhi2.Write()
    hjtt_dPhi3.Write()
    
    hjtt_nJetsPU.Write()
    hjtt_bJetsPU.Write()
    hjtt_jet1PtPU.Write()
    hjtt_jet2PtPU.Write()
    hjtt_metPU  .Write()
    hjtt_mhtPU  .Write()
    hjtt_bestTopJetMassPU.Write()
    hjtt_TbestTopJetPU   .Write()
    hjtt_TbJetPU         .Write()
    hjtt_mt2PU           .Write()
    hjtt_dPhi1PU.Write()
    hjtt_dPhi2PU.Write()
    hjtt_dPhi3PU.Write()

    ###
    hjnob_nJets.Write()
    hjnob_bJets.Write()
    hjnob_jet1Pt.Write()
    hjnob_jet2Pt.Write()
    hjnob_met  .Write()
    hjnob_mht  .Write()
    hjnob_bestTopJetMass.Write()
    hjnob_TbestTopJet   .Write()
    hjnob_TbJet         .Write()
    hjnob_mt2           .Write()
    hjnob_dPhi1.Write()
    hjnob_dPhi2.Write()
    hjnob_dPhi3.Write()
    
    hjnob_nJetsPU.Write()
    hjnob_bJetsPU.Write()
    hjnob_jet1PtPU.Write()
    hjnob_jet2PtPU.Write()
    hjnob_metPU  .Write()
    hjnob_mhtPU  .Write()
    hjnob_bestTopJetMassPU.Write()
    hjnob_TbestTopJetPU   .Write()
    hjnob_TbJetPU         .Write()
    hjnob_mt2PU           .Write()
    hjnob_dPhi1PU.Write()
    hjnob_dPhi2PU.Write()
    hjnob_dPhi3PU.Write()

    ###
    hjnobtt_nJets.Write()
    hjnobtt_bJets.Write()
    hjnobtt_jet1Pt.Write()
    hjnobtt_jet2Pt.Write()
    hjnobtt_met  .Write()
    hjnobtt_mht  .Write()
    hjnobtt_bestTopJetMass.Write()
    hjnobtt_TbestTopJet   .Write()
    hjnobtt_TbJet         .Write()
    hjnobtt_mt2           .Write()
    hjnobtt_dPhi1.Write()
    hjnobtt_dPhi2.Write()
    hjnobtt_dPhi3.Write()
    
    hjnobtt_nJetsPU.Write()
    hjnobtt_bJetsPU.Write()
    hjnobtt_jet1PtPU.Write()
    hjnobtt_jet2PtPU.Write()
    hjnobtt_metPU  .Write()
    hjnobtt_mhtPU  .Write()
    hjnobtt_bestTopJetMassPU.Write()
    hjnobtt_TbestTopJetPU   .Write()
    hjnobtt_TbJetPU         .Write()
    hjnobtt_mt2PU           .Write()
    hjnobtt_dPhi1PU.Write()
    hjnobtt_dPhi2PU.Write()
    hjnobtt_dPhi3PU.Write()

    outputFile.Write()
    outputFile.Close()

##########################
if __name__ == '__main__':
    main()
    
