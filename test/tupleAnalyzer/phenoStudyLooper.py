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
    myChain = r.TChain('%s/RA2Values'%(options.treeName))

    if options.debug:
        if options.sample=="gjets":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v6/res/*_1_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v6/res/*_1_?_???.root")
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v6/res/*_1_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v6/res/*_1_?_???.root")
    else:
        if options.sample=="gjets":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
    
    ###variables
    nJets  = []
    ht     = []
    mht    = []
    htvsmht     = []
    genBosonPt = []

    nJetsminv  = []
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
                        ht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        mht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        htvsmht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        genBosonPt[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                        
                        nJetsminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
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
        if ra2_MHT > 200 and event.ra2_genDPhiMHT1 > 0.5 and event.ra2_genDPhiMHT2 > 0.5 and event.ra2_genDPhiMHT3 > 0.3:
            for ad,anDir in enumerate(analysisDirs):                
                if ad > 1:
                    for j,jetdir in enumerate(njbins):
                        outputFile.cd()
                        outputFile.cd(anDir[1]+"/"+jetdir[0])
                        for h,htdir in enumerate(htbins):
                            for m,mhtdir in enumerate(mhtbins):
                                if ra2_MHT > mhtdir[1][0] and ra2_MHT < mhtdir[1][1]:
                                    
                                    # print "%s %s: #J:%d #JI:%d pT:%2.2f"%(anDir[1],jetdir[0],countJets,countJetsMInv,ra2_genBosonPt)
                                    if (cutF.phenoCuts(event,ad,options.minpt,countJetsMInv,jetdir[1],ra2_HTMInv,htdir[1],mhtdir[1])):
                                        nJetsminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]  .Fill(ra2_nJetsPt50Eta25MInv,ra2_EventWt)
                                        htminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]     .Fill(ra2_HTMInv            ,ra2_EventWt)
                                        mhtminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_MHT           ,ra2_EventWt)
                                        htminvvsmht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_HTMInv,ra2_MHT    ,ra2_EventWt)
                                        genBosonPtminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_genBosonPt ,ra2_EventWt)
                                        
                                    if (cutF.phenoCuts(event,ad,options.minpt,countJets,jetdir[1],ra2_HT,htdir[1],mhtdir[1])):
                                        nJets[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]  .Fill(ra2_nJetsPt50Eta25,ra2_EventWt)
                                        ht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]     .Fill(ra2_HT            ,ra2_EventWt)
                                        mht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_MHT           ,ra2_EventWt)
                                        htvsmht[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_HT,ra2_MHT,ra2_EventWt)
                                        genBosonPt[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_genBosonPt ,ra2_EventWt)
                                        
                            
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
                        nJetsminv[(0*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
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
    
