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
    parser.add_option('-j', action="store",      default=5,     dest="numJets", type="int")
    parser.add_option('-s', action="store",      default="gjets",dest="sample", type="string")
    parser.add_option('-f', action="store",      default="0",    dest="subsec", type="int")
    parser.add_option('-t', action="store",      default="directPhotonsID",dest="treeName", type="string")
#    parser.add_option('-i', action="store",      default="ID",dest="idiso", type="string")
    parser.add_option('-c', action="store",      default=100.0,dest="minpt", type="float")
#    parser.add_option('-p', action="store_true", default=True, dest="dphi")
    options, args = parser.parse_args()

    r.gROOT.SetBatch(True)
    #debug = False
    myWorkingDir = os.getcwd()
    
    subfiles = [
        "*",#0
        "*_?_?_???",#1
        "*_1?_?_???","*_2?_?_???","*_3?_?_???","*_4?_?_???","*_5?_?_???","*_6?_?_???","*_7?_?_???","*_8?_?_???","*_9?_?_???",#10
        "*_11?_?_???","*_12?_?_???","*_13?_?_???","*_14?_?_???","*_15?_?_???","*_16?_?_???","*_17?_?_???","*_18?_?_???","*_19?_?_???",#19
        "*_21?_?_???","*_22?_?_???","*_23?_?_???","*_24?_?_???","*_25?_?_???","*_26?_?_???","*_27?_?_???","*_28?_?_???","*_29?_?_???",#28
        "*_31?_?_???","*_32?_?_???","*_33?_?_???","*_34?_?_???","*_35?_?_???","*_36?_?_???","*_37?_?_???","*_38?_?_???","*_39?_?_???",#37
        "*_41?_?_???","*_42?_?_???","*_43?_?_???","*_44?_?_???","*_45?_?_???","*_46?_?_???","*_47?_?_???","*_48?_?_???","*_49?_?_???",#46
        "*_51?_?_???","*_52?_?_???","*_53?_?_???","*_54?_?_???","*_55?_?_???","*_56?_?_???","*_57?_?_???","*_58?_?_???","*_59?_?_???",#55
        "*_61?_?_???","*_62?_?_???","*_63?_?_???","*_64?_?_???","*_65?_?_???","*_66?_?_???","*_67?_?_???","*_68?_?_???","*_69?_?_???",#64
        "*_71?_?_???","*_72?_?_???","*_73?_?_???","*_74?_?_???","*_75?_?_???","*_76?_?_???","*_77?_?_???","*_78?_?_???","*_79?_?_???",#73
        "*_81?_?_???","*_82?_?_???","*_83?_?_???","*_84?_?_???","*_85?_?_???","*_86?_?_???","*_87?_?_???","*_88?_?_???","*_89?_?_???",#82
        "*_91?_?_???","*_92?_?_???","*_93?_?_???","*_94?_?_???","*_95?_?_???","*_96?_?_???","*_97?_?_???","*_98?_?_???","*_99?_?_???",#91
        "*_101?_?_???","*_102?_?_???","*_103?_?_???","*_104?_?_???","*_105?_?_???","*_106?_?_???","*_107?_?_???","*_108?_?_???","*_109?_?_???",#100
        ]
    outFileName = "phenoPlots_%s_%d.root"%(options.sample,options.subsec)
    if options.debug:
        outFileName = "phenoPlots_%s_test.root"%(options.sample)

    print outFileName
    sys.stdout.flush()
    outputFile  = r.TFile(outFileName,"RECREATE")
    precutsDir  = outputFile.mkdir("precuts")
    dphiDir     = outputFile.mkdir("dphi")
    baselineDir = outputFile.mkdir("baseline")

    analysisDirs = [
        [precutsDir, "precuts"],
        [dphiDir,    "dphi"],
        [baselineDir,"baseline"],
        ]
    for anDir in analysisDirs:
        inclusiveDir  = anDir[0].mkdir("inclusive")
        inclusive3Dir = anDir[0].mkdir("inclusive3")
        jet2Dir       = anDir[0].mkdir("2jets")
        jet3Dir       = anDir[0].mkdir("3jets")
        jet4Dir       = anDir[0].mkdir("4jets")
        jet5Dir       = anDir[0].mkdir("5jets")
        jet6Dir       = anDir[0].mkdir("6jets")
        jet7Dir       = anDir[0].mkdir("7jets")
        jet8Dir       = anDir[0].mkdir("8jets")
        jet3to5Dir    = anDir[0].mkdir("3to5jets")
        jet6to7Dir    = anDir[0].mkdir("6to7jets")

    njbins = [
#        ["precuts"  ,2,100],
        ["inclusive", [2,100]],
        ["inclusive3",[3,100]],
        ["2jets"     ,[2,2]],
        ["3jets"     ,[3,3]],
        ["4jets"     ,[4,4]],
        ["5jets"     ,[5,5]],
        ["6jets"     ,[6,6]],
        ["7jets"     ,[7,7]],
        ["8jets"     ,[8,100]],
        ["3to5jets"  ,[3,5]],
        ["6to7jets"  ,[6,7]],
        
        ]
    htbins = [
        ["none",     [0,10000]],
        ["reduced",  [350,10000]],
        ["inclusive",[500,10000]],
        ["bin1",     [500,900]],
        ["bin2",     [900,1300]],
        ["bin3",     [1300,10000]],
        ]
    mhtbins = [
        ["none",     [0,10000]],
        ["reduced",  [150,10000]],
        ["inclusive",[200,10000]],
        ["bin1",     [200,350]],
        ["bin2",     [350,500]],
        ["bin3",     [500,10000]],
        ]
    myChain = r.TChain('%s/RA2Values'%(options.treeName))

    if options.debug:
        if options.sample=="gjets":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v6/res/*_1_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v6/res/*_1_?_???.root")
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v6/res/*_1_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v6/res/*_1_?_???.root")
        elif options.sample=="zllht":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/dyjetstoll_ht200_gen_tree_v8/res/*_?_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/dyjetstoll_ht400_gen_tree_v8/res/*_?_?_???.root")
        elif options.sample=="zllm50":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/dyjetstoll_m50_gen_tree_v8/res/*_?_?_???.root")
    else:
        if options.sample=="gjets":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zllht":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/dyjetstoll_ht200_gen_tree_v8/res/%s.root"%(subfiles[options.subsec]))
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/dyjetstoll_ht400_gen_tree_v8/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zllm50":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/dyjetstoll_m50_gen_tree_v8/res/%s.root"%(subfiles[options.subsec]))
    
    ###variables
    nJets  = []
    bJets  = []
    ht     = []
    mht    = []
    htvsmht     = []
    genBosonPt = []

    nJetsminv  = []
    bJetsminv  = []
    htminv = []
    mhtminv = []
    htminvvsmht = []
    genBosonPtminv = []

    for ad,anDir in enumerate(analysisDirs):
        if ad > 1:
            for j,jetdir in enumerate(njbins):
                for h,htdir in enumerate(htbins):
                    for m,mhtdir in enumerate(mhtbins):
                        outputFile.cd()
                        outputFile.cd(anDir[1]+"/"+jetdir[0])
                        nJets     .append(r.TH1D(    "h_%s_ht%s_mht%s_nJets_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_nJets_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 15, -0.5, 14.5))
                        nJetsminv .append(r.TH1D("h_%s_ht%s_mht%s_nJetsminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_nJetsminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 15, -0.5, 14.5))

                        bJets     .append(r.TH1D(    "h_%s_ht%s_mht%s_bJets_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_bJets_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 15, -0.5, 14.5))
                        bJetsminv .append(r.TH1D("h_%s_ht%s_mht%s_bJetsminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_bJetsminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 15, -0.5, 14.5))
                        
                        ht     .append(r.TH1D(    "h_%s_ht%s_mht%s_ht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    "h_%s_ht%s_mht%s_ht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    30,  0, 3000))
                        htminv .append(r.TH1D("h_%s_ht%s_mht%s_htminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_htminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),30,  0, 3000))
                        
                        mht    .append(r.TH1D(    "h_%s_ht%s_mht%s_mht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   "h_%s_ht%s_mht%s_mht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   20,  0, 1000))
                        mhtminv.append(r.TH1D("h_%s_ht%s_mht%s_mhtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   "h_%s_ht%s_mht%s_mhtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   20,  0, 1000))
                        
                        htvsmht    .append(r.TH2D(    "h_%s_ht%s_mht%s_htvsmht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    "h_%s_ht%s_mht%s_htvsmht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    30,  0, 3000, 20, 0, 1000))
                        htminvvsmht.append(r.TH2D("h_%s_ht%s_mht%s_htminvvsmht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_htminvvsmht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),30,  0, 3000, 20, 0, 1000))
                        
                        genBosonPt    .append(r.TH1D(    "h_%s_ht%s_mht%s_genBosonPt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),     "h_%s_ht%s_mht%s_genBosonPt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  40,   0, 1000))
                        genBosonPtminv.append(r.TH1D("h_%s_ht%s_mht%s_genBosonPtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_genBosonPtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  40,   0, 1000))
            
            
    ####
    for ad,anDir in enumerate(analysisDirs):
        if ad > 1:
            for j,jetdir in enumerate(njbins):
                for h,htdir in enumerate(htbins):
                    for m,mhtdir in enumerate(mhtbins):
                        #for sd,subdir in enumerate(dirs):
                        print str(ad)+"*"+str(len(njbins))+"*"+str(len(htbins))+"*"+str(len(mhtbins))+"+"+str(j)+"*"+str(len(htbins))+"*"+str(len(mhtbins))+"+"+str(h)+"*"+str(len(mhtbins))+"+"+str(m)
                        print nJets[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]
                        # (0*sdsize)+sd
                        nJets[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        bJets[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        ht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        mht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        htvsmht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        genBosonPt[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        
                        nJetsminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        bJetsminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        htminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        mhtminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        htminvvsmht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        genBosonPtminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()



    ##################
    fChain = myChain

    ###Timing information
    decade  = 0
    century = 0
    tsw = r.TStopwatch()
    tenpcount = 1
    onepcount = 1

    countPT50  = 0
    countPT100 = 0
    countPT150 = 0
    countPT175 = 0

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
        ra2_nJetsPt50Eta25      = event.ra2_nJetsGenPt50Eta25 
        ra2_bJetsPt30Eta24      = event.ra2_bJetsGenPt30Eta24
        ra2_nJetsPt50Eta25MInv  = event.ra2_nJetsGenPt50Eta25MInv
        ra2_PUWt            = event.ra2_PUWt           
        ra2_EventWt         = event.ra2_EventWt        
        ra2_HT              = event.ra2_genHT            
        ra2_HTMInv          = event.ra2_genHTMInv            
        ra2_MHT             = event.ra2_genMHT            
        ra2_genBosonPt      = event.ra2_genBoson1Pt         

        if ra2_nJetsPt50Eta25>1 or ra2_nJetsPt50Eta25MInv>1:
            if ra2_genBosonPt>50:
                countPT50 = countPT50 + 1
            if ra2_genBosonPt>100:
                countPT100 = countPT100 + 1
            if ra2_genBosonPt>150:
                countPT150 = countPT150 + 1
            if ra2_genBosonPt>175:
                countPT175 = countPT175 + 1

        ####Fill plots
        countJetsMInv = event.ra2_nJetsGenPt50Eta25MInv
        countJets = event.ra2_nJetsGenPt50Eta25
        #if ra2_MHT > 200 and event.ra2_genDPhiMHT1 > 0.5 and event.ra2_genDPhiMHT2 > 0.5 and event.ra2_genDPhiMHT3 > 0.3:
        if event.ra2_genDPhiMHT1 > 0.5 and event.ra2_genDPhiMHT2 > 0.5 and event.ra2_genDPhiMHT3 > 0.3:
            for ad,anDir in enumerate(analysisDirs):                
                #if ad > 1:
                for j,jetdir in enumerate(njbins):
                    outputFile.cd()
                    outputFile.cd(anDir[1]+"/"+jetdir[0])
                    for h,htdir in enumerate(htbins):
                        for m,mhtdir in enumerate(mhtbins):
                            if ra2_MHT > mhtdir[1][0] and ra2_MHT < mhtdir[1][1]:
                                # print "%s %s: #J:%d #JI:%d pT:%2.2f"%(anDir[1],jetdir[0],countJets,countJetsMInv,ra2_genBosonPt)
                                if (cutF.phenoCuts(event,ad,options.minpt,countJetsMInv,jetdir[1],ra2_HTMInv,htdir[1],mhtdir[1])):
                                    nJetsminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]     .Fill(ra2_nJetsPt50Eta25MInv,ra2_EventWt)
                                    bJetsminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]     .Fill(ra2_bJetsPt30Eta24    ,ra2_EventWt)
                                    htminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]        .Fill(ra2_HTMInv            ,ra2_EventWt)
                                    mhtminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]       .Fill(ra2_MHT               ,ra2_EventWt)
                                    htminvvsmht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]   .Fill(ra2_HTMInv,ra2_MHT    ,ra2_EventWt)
                                    genBosonPtminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_genBosonPt        ,ra2_EventWt)
                                    
                                if (cutF.phenoCuts(event,ad,options.minpt,countJets,jetdir[1],ra2_HT,htdir[1],mhtdir[1])):
                                    nJets[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]     .Fill(ra2_nJetsPt50Eta25,ra2_EventWt)
                                    bJets[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]     .Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                                    ht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]        .Fill(ra2_HT            ,ra2_EventWt)
                                    mht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]       .Fill(ra2_MHT           ,ra2_EventWt)
                                    htvsmht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]   .Fill(ra2_HT,ra2_MHT    ,ra2_EventWt)
                                    genBosonPt[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_genBosonPt    ,ra2_EventWt)
                                    
                            
        #########
        i = i + 1

    #####        

    print "50:%d 100:%d 150:%d 175:%d"%(countPT50,countPT100,countPT150,countPT175)

    outputFile.cd()
    
    ####
    for ad,anDir in enumerate(analysisDirs):
        if ad > 1:
            for j,jetdir in enumerate(njbins):
                for h,htdir in enumerate(htbins):
                    for m,mhtdir in enumerate(mhtbins):
                        outputFile.cd()
                        outputFile.cd(anDir[1]+"/"+jetdir[0])
                        nJets[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                        bJets[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                        nJetsminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                        bJetsminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                        ht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                        htminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                        mht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                        mhtminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                        htvsmht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                        htminvvsmht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                        genBosonPt[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                        genBosonPtminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()

    outputFile.Write()
    outputFile.Close()

##########################
if __name__ == '__main__':
    main()
    
