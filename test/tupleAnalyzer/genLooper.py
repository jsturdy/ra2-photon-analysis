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
    parser.add_option('-j', action="store",      default=5,     dest="numJets",type="int")
    parser.add_option('-s', action="store",      default="gjets",dest="sample", type="string")
    parser.add_option('-f', action="store",      default="0",    dest="subsec", type="int")
    parser.add_option('-t', action="store",      default="directPhotonsID",dest="treeName", type="string")
    parser.add_option('-i', action="store",      default="ID",dest="idiso", type="string")
    parser.add_option('-c', action="store",      default=100.0,dest="minpt", type="float")
    options, args = parser.parse_args()

    r.gROOT.SetBatch(True)
    #debug = False
    myWorkingDir = os.getcwd()

    subfiles = [
        "*",
        "*_?_?_???",
        "*_??_?_???",
        "*_1??_?_???",
        "*_2??_?_???",
        "*_3??_?_???",
        "*_4??_?_???",
        "*_5??_?_???",
        "*_6??_?_???",
        "*_7??_?_???",
        "*_8??_?_???",
        "*_9??_?_???",
        "*_10??_?_???",
        ]
    outFileName = "genLevelPlots_%s_%s_%d.root"%(options.sample,options.idiso,options.subsec)
    if options.debug:
        outFileName = "genLevelPlots_%s_test_%s.root"%(options.sample,options.idiso)

    print outFileName
    sys.stdout.flush()
    outputFile = r.TFile(outFileName,"RECREATE")
    precutsDir    = outputFile.mkdir("precuts")
    acceptanceDir = outputFile.mkdir("acceptance")
    #recoidDir     = outputFile.mkdir("reco%s"%(options.idiso))
    recoidDir     = outputFile.mkdir("recoid")
    recopfisoDir  = outputFile.mkdir("recopfiso")
    postcutsDir   = outputFile.mkdir("postcuts")

    analysisDirs = [
        [precutsDir   ,"precuts"],
        [acceptanceDir,"acceptance"],
        #[recoisoDir   ,"reco%s"%(options.idiso)],
        [recoidDir    ,"recoid"],
        [recopfisoDir ,"recopfiso"],
        [postcutsDir  ,"postcuts"],

        ]
    for anDir in analysisDirs:
        inclusiveDir = anDir[0].mkdir("inclusive")
        jet2Dir      = anDir[0].mkdir("2jets")
        jet3Dir      = anDir[0].mkdir("3jets")
        jet4Dir      = anDir[0].mkdir("4jets")
        jet5Dir      = anDir[0].mkdir("5jets")
        jet6Dir      = anDir[0].mkdir("6jets")
        jet7Dir      = anDir[0].mkdir("7jets")
        jet8Dir      = anDir[0].mkdir("8jets")
        jet3to5Dir   = anDir[0].mkdir("3to5jets")
        jet6to7Dir   = anDir[0].mkdir("6to7jets")

    njbins = [
#        ["precuts"  ,2,100],
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
    htbins = [
        ["inclusive",[500,10000]],
        ["bin1",[500,900]],
        ["bin2",[900,1300]],
        ["bin3",[1300,10000]],
        ]
    mhtbins = [
        ["inclusive",[200,10000]],
        ["bin1",[200,350]],
        ["bin2",[350,500]],
        ["bin3",[500,10000]],
        ]
    #idChain = r.TChain('%s/RA2Values'%(options.treeName))
    idChain = r.TChain('directPhotonsID/RA2Values')
    pfisoChain = r.TChain('directPhotonsIDPFIso/RA2Values')

    if options.debug:
        if options.sample=="gjets":
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v6/res/*_?_?_???.root")
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v6/res/*_?_?_???.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v6/res/*_?_?_???.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v6/res/*_?_?_???.root")
        elif options.sample=="zinv":
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v6/res/*_?_?_???.root")
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v6/res/*_?_?_???.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v6/res/*_?_?_???.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v6/res/*_?_?_???.root")
    else:
        if options.sample=="gjets":
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv":
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v6/res/*.root")
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v6/res/*.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v6/res/*.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v6/res/*.root")
    
    ###variables
    nJets  = []
    bJets  = []
    jet1Pt = []
    jet2Pt = []
    ht     = []
    mht    = []
    htvsmht     = []
    dPhi1  = []
    dPhi2  = []
    dPhi3  = []
    genBosonPt = []
    genBosonEta = []
    bosonPt = []
    bosonEta = []

    nJetsminv  = []
    bJetsminv  = []
    jet1Ptminv = []
    jet2Ptminv = []
    htminv  = []
    mhtminv = []
    htminvvsmht = []
    dPhi1minv  = []
    dPhi2minv  = []
    dPhi3minv  = []
    genBosonPtminv  = []
    genBosonEtaminv = []
    bosonPtminv  = []
    bosonEtaminv = []

    for ad,anDir in enumerate(analysisDirs):
        for j,jetdir in enumerate(njbins):
            for h,htdir in enumerate(htbins):
                for m,mhtdir in enumerate(mhtbins):
                    outputFile.cd()
                    outputFile.cd(anDir[1]+"/"+subdir[0])
                    
                    nJets .append(r.TH1D("h_%s_ht%s_mht%s_nJets_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_nJets_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 15, -0.5, 14.5))
                    bJets .append(r.TH1D("h_%s_ht%s_mht%s_bJets_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_bJets_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 10,  -0.5, 9.5))
                    jet1Pt.append(r.TH1D("h_%s_ht%s_mht%s_jet1Pt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet1Pt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),2*50,  0, 1000))
                    jet2Pt.append(r.TH1D("h_%s_ht%s_mht%s_jet2Pt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet2Pt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),2*50,  0, 500))
                    ht    .append(r.TH1D("h_%s_ht%s_mht%s_ht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    "h_%s_ht%s_mht%s_ht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    30,  0, 3000))
                    mht   .append(r.TH1D("h_%s_ht%s_mht%s_mht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   "h_%s_ht%s_mht%s_mht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   20,  0, 1000))
                    htvsmht.append(r.TH2D("h_%s_ht%s_mht%s_htvsmht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_htvsmht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),30,  0, 3000, 20, 0, 1000))
                    dPhi1 .append(r.TH1D("h_%s_ht%s_mht%s_dPhi1_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_dPhi1_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 3.2))
                    dPhi2 .append(r.TH1D("h_%s_ht%s_mht%s_dPhi2_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_dPhi2_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 3.2))
                    dPhi3 .append(r.TH1D("h_%s_ht%s_mht%s_dPhi3_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_dPhi3_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 3.2))
                    genBosonPt .append(r.TH1D("h_%s_ht%s_mht%s_genBosonPt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_genBosonPt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  100,   0, 1000))
                    genBosonEta.append(r.TH1D("h_%s_ht%s_mht%s_genBosonEta_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_genBosonEta_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),2*50,  -5, 5))
                    bosonPt .append(r.TH1D("h_%s_ht%s_mht%s_bosonPt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_bosonPt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  100,   0, 1000))
                    bosonEta.append(r.TH1D("h_%s_ht%s_mht%s_bosonEta_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_bosonEta_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),2*50,  -5, 5))
                    
                    nJetsminv .append(r.TH1D("h_%s_ht%s_mht%s_nJetsminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_nJetsminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 15, -0.5, 14.5))
                    bJetsminv .append(r.TH1D("h_%s_ht%s_mht%s_bJetsminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_bJetsminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 10,  -0.5, 9.5))
                    jet1Ptminv.append(r.TH1D("h_%s_ht%s_mht%s_jet1Ptminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet1Ptminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),2*50,  0, 1000))
                    jet2Ptminv.append(r.TH1D("h_%s_ht%s_mht%s_jet2Ptminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet2Ptminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),2*50,  0, 500))
                    htminv.append(r.TH1D("h_%s_ht%s_mht%s_htminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_htminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),30,  0, 3000))
                    mhtminv.append(r.TH1D("h_%s_ht%s_mht%s_mhtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_mhtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),30,  0, 3000))
                    htminvvsmht.append(r.TH2D("h_%s_ht%s_mht%s_htminvvsmht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_htminvvsmht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),30,  0, 3000, 20, 0, 1000))
                    dPhi1minv .append(r.TH1D("h_%s_ht%s_mht%s_dPhi1minv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_dPhi1minv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 3.2))
                    dPhi2minv .append(r.TH1D("h_%s_ht%s_mht%s_dPhi2minv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_dPhi2minv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 3.2))
                    dPhi3minv .append(r.TH1D("h_%s_ht%s_mht%s_dPhi3minv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_dPhi3minv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 3.2))
                    genBosonPtminv .append(r.TH1D("h_%s_ht%s_mht%s_genBosonPtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_genBosonPtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  100,   0, 1000))
                    genBosonEtaminv.append(r.TH1D("h_%s_ht%s_mht%s_genBosonEtaminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_genBosonEtaminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),2*50,  -5, 5))
                    bosonPtminv .append(r.TH1D("h_%s_ht%s_mht%s_bosonPtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_bosonPtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  100,   0, 1000))
                    bosonEtaminv.append(r.TH1D("h_%s_ht%s_mht%s_bosonEtaminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_bosonEtaminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),2*50,  -5, 5))
                    
    ####
                    #(ad*sdsize)+sd
                    nJets[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bJets[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet1Pt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet2Pt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    ht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    mht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    htvsmht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    dPhi1[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    dPhi2[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    dPhi3[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonPt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonEta[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonPt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonEta[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    
                    nJetsminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bJetsminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet1Ptminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet2Ptminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    htminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    mhtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    htminvvsmht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    dPhi1minv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    dPhi2minv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    dPhi3minv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonPtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonEtaminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonPtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonEtaminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()



    ##################
    #fChain = myChain
    #idChain    = idChain
    #pfisoChain = myChain

    ###Timing information
    decade  = 0
    century = 0
    tsw = r.TStopwatch()
    tenpcount = 1
    onepcount = 1


    nentries = idChain.GetEntries()
    print "nentries %d"%(nentries)
    sys.stdout.flush()
    i = 0
    #for i,event in enumerate(fChain):
    for idevent,pfisoevent in itertools.izip(idChain,pfisoChain):

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
        ra2_Vertices           = idevent.ra2_Vertices      
        ra2_nJetsPt50Eta25     = idevent.ra2_nJetsGenPt50Eta25 
        ra2_nJetsPt50Eta25MInv = idevent.ra2_nJetsGenPt50Eta25MInv
        ra2_bJetsPt30Eta24     = idevent.ra2_bJetsGenPt30Eta24 
        ra2_PUWt               = idevent.ra2_PUWt           
        ra2_EventWt            = idevent.ra2_EventWt        
        ra2_HT                 = idevent.ra2_genHT            
        ra2_HTMInv             = idevent.ra2_genHTMInv            
        ra2_MHT                = idevent.ra2_genMHT            
        ra2_MET                = idevent.ra2_genMET            
        ra2_dPhi1              = idevent.ra2_genDPhiMHT1
        ra2_dPhi2              = idevent.ra2_genDPhiMHT2
        ra2_dPhi3              = idevent.ra2_genDPhiMHT3
        ra2_Jet1Pt             = idevent.ra2_genJet1Pt         
        ra2_Jet2Pt             = idevent.ra2_genJet2Pt         
        ra2_genBosonPt         = idevent.ra2_genBoson1Pt         
        ra2_genBosonEta        = idevent.ra2_genBoson1Eta         
        ra2_bosonIDPt          = idevent.ra2_boson1Pt         
        ra2_bosonIDEta         = idevent.ra2_boson1Eta         
        ra2_bosonIDPFIsoPt     = pfisoevent.ra2_boson1Pt         
        ra2_bosonIDPFIsoEta    = pfisoevent.ra2_boson1Eta         

        ra2_bosonPt  = 0
        ra2_bosonEta = 0

        ####Fill plots
        for ad,anDir in enumerate(analysisDirs):
            if ad > 2:
                ra2_bosonPt  = ra2_bosonIDPFIsoPt
                ra2_bosonEta = ra2_bosonIDPFIsoEta
                recoPt       = ra2_bosonIDPFIsoPt
                passReco     = pfisoevent.ra2_genPassRecoIso
            else:
                ra2_bosonPt  = ra2_bosonIDPt
                ra2_bosonEta = ra2_bosonIDEta
                recoPt       = ra2_bosonIDPt
                passReco     = idevent.ra2_genPassRecoIso

            for sd,subdir in enumerate(dirs):
                outputFile.cd()
                outputFile.cd(anDir[1]+"/"+subdir[0])

                if (cutF.genCuts(idevent,ad,options.minpt,recoPt,passReco,ra2_HT,ra2_nJetsPt50Eta25,subdir[1])):
                    nJets[ad*len(dirs)+sd]  .Fill(ra2_nJetsPt50Eta25,ra2_EventWt)
                    bJets[ad*len(dirs)+sd]  .Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                    jet1Pt[ad*len(dirs)+sd] .Fill(ra2_Jet1Pt        ,ra2_EventWt)
                    jet2Pt[ad*len(dirs)+sd] .Fill(ra2_Jet2Pt        ,ra2_EventWt)
                    ht[ad*len(dirs)+sd]     .Fill(ra2_HT            ,ra2_EventWt)
                    mht[ad*len(dirs)+sd]    .Fill(ra2_MHT           ,ra2_EventWt)
                    htvsmht[ad*len(dirs)+sd].Fill(ra2_HT,ra2_MHT    ,ra2_EventWt)
                    dPhi1[ad*len(dirs)+sd]  .Fill(ra2_dPhi1         ,ra2_EventWt)
                    dPhi2[ad*len(dirs)+sd]  .Fill(ra2_dPhi2         ,ra2_EventWt)
                    dPhi3[ad*len(dirs)+sd]  .Fill(ra2_dPhi3         ,ra2_EventWt)
                    
                    genBosonPt[ad*len(dirs)+sd] .Fill(ra2_genBosonPt ,ra2_EventWt)
                    genBosonEta[ad*len(dirs)+sd].Fill(ra2_genBosonEta,ra2_EventWt)
                    bosonPt[ad*len(dirs)+sd]    .Fill(ra2_bosonPt    ,ra2_EventWt)
                    bosonEta[ad*len(dirs)+sd]   .Fill(ra2_bosonEta   ,ra2_EventWt)

                if (cutF.genCuts(idevent,ad,options.minpt,recoPt,passReco,ra2_HTMInv,ra2_nJetsPt50Eta25MInv,subdir[1])):
                    nJetsminv[ad*len(dirs)+sd]  .Fill(ra2_nJetsPt50Eta25MInv,ra2_EventWt)
                    bJetsminv[ad*len(dirs)+sd]  .Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                    jet1Ptminv[ad*len(dirs)+sd] .Fill(ra2_Jet1Pt        ,ra2_EventWt)
                    jet2Ptminv[ad*len(dirs)+sd] .Fill(ra2_Jet2Pt        ,ra2_EventWt)
                    htminv[ad*len(dirs)+sd]     .Fill(ra2_HTMInv        ,ra2_EventWt)
                    mhtminv[ad*len(dirs)+sd]    .Fill(ra2_MHT           ,ra2_EventWt)
                    htminvvsmht[ad*len(dirs)+sd].Fill(ra2_HTMInv,ra2_MHT,ra2_EventWt)
                    dPhi1minv[ad*len(dirs)+sd]  .Fill(ra2_dPhi1         ,ra2_EventWt)
                    dPhi2minv[ad*len(dirs)+sd]  .Fill(ra2_dPhi2         ,ra2_EventWt)
                    dPhi3minv[ad*len(dirs)+sd]  .Fill(ra2_dPhi3         ,ra2_EventWt)
                    
                    genBosonPtminv[ad*len(dirs)+sd] .Fill(ra2_genBosonPt ,ra2_EventWt)
                    genBosonEtaminv[ad*len(dirs)+sd].Fill(ra2_genBosonEta,ra2_EventWt)
                    bosonPtminv[ad*len(dirs)+sd]    .Fill(ra2_bosonPt    ,ra2_EventWt)
                    bosonEtaminv[ad*len(dirs)+sd]   .Fill(ra2_bosonEta   ,ra2_EventWt)
        
                
        #########
        i = i + 1

    #####        
    outputFile.cd()

    ####
    for ad,anDir in enumerate(analysisDirs):
        for sd,subdir in enumerate(dirs):
            outputFile.cd()
            outputFile.cd(anDir[1]+"/"+subdir[0])

            nJets[ad*len(dirs)+sd].Write()
            bJets[ad*len(dirs)+sd].Write()
            jet1Pt[ad*len(dirs)+sd].Write()
            jet2Pt[ad*len(dirs)+sd].Write()
            ht[ad*len(dirs)+sd].Write()
            mht[ad*len(dirs)+sd].Write()
            htvsmht[ad*len(dirs)+sd].Write()
            dPhi1[ad*len(dirs)+sd].Write()
            dPhi2[ad*len(dirs)+sd].Write()
            dPhi3[ad*len(dirs)+sd].Write()
            genBosonPt[ad*len(dirs)+sd].Write()
            genBosonEta[ad*len(dirs)+sd].Write()
            bosonPt[ad*len(dirs)+sd].Write()
            bosonEta[ad*len(dirs)+sd].Write()

            nJetsminv[ad*len(dirs)+sd].Write()
            bJetsminv[ad*len(dirs)+sd].Write()
            jet1Ptminv[ad*len(dirs)+sd].Write()
            jet2Ptminv[ad*len(dirs)+sd].Write()
            htminv[ad*len(dirs)+sd].Write()
            mhtminv[ad*len(dirs)+sd].Write()
            htminvvsmht[ad*len(dirs)+sd].Write()
            dPhi1minv[ad*len(dirs)+sd].Write()
            dPhi2minv[ad*len(dirs)+sd].Write()
            dPhi3minv[ad*len(dirs)+sd].Write()
            genBosonPtminv[ad*len(dirs)+sd].Write()
            genBosonEtaminv[ad*len(dirs)+sd].Write()
            bosonPtminv[ad*len(dirs)+sd].Write()
            bosonEtaminv[ad*len(dirs)+sd].Write()

    outputFile.Write()
    outputFile.Close()

##########################
if __name__ == '__main__':
    main()
    
