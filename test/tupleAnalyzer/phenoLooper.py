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

    dirs = [
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
    myChain = r.TChain('%s/RA2Values'%(options.treeName))

    if options.debug:
        if options.sample=="gjets":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v6/res/*_?_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v6/res/*_?_?_???.root")
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v6/res/*_?_?_???.root")
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v6/res/*_?_?_???.root")
    else:
        if options.sample=="gjets":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v6/res/%s.root"%(subfiles[options.subsec]))
    
    ###variables
    nVertices     = []
    nVerticesReWt = []

    nJets  = []
    bJets  = []
    jet1Pt = []
    jet2Pt = []
    ht     = []
    met    = []
    mht    = []

    htvsmht     = []
    htvsmet     = []

    dPhi1  = []
    dPhi2  = []
    dPhi3  = []

    genBosonPt = []
    genBosonEta = []

    nJetsminv  = []
    htminv = []
    htminvvsmht = []
    htminvvsmet = []

    for ad,anDir in enumerate(analysisDirs):
        for sd,subdir in enumerate(dirs):
            outputFile.cd()
            outputFile.cd(anDir[1]+"/"+subdir[0])
            nVertices    .append(r.TH1D("h_%s_nVertices_%s"%(subdir[0],anDir[1]),    "h_%s_nVertices_%s"%(subdir[0],anDir[1]),    50, -0.5, 49.5))
            nVerticesReWt.append(r.TH1D("h_%s_nVerticesReWt_%s"%(subdir[0],anDir[1]),"h_%s_nVerticesReWt_%s"%(subdir[0],anDir[1]),50, -0.5, 49.5))
            
            nJets .append(r.TH1D("h_%s_nJets_%s"%(subdir[0],anDir[1]), "h_%s_nJets_%s"%(subdir[0],anDir[1]), 15, -0.5, 14.5))
            nJetsminv .append(r.TH1D("h_%s_nJetsminv_%s"%(subdir[0],anDir[1]), "h_%s_nJetsminv_%s"%(subdir[0],anDir[1]), 15, -0.5, 14.5))
            bJets .append(r.TH1D("h_%s_bJets_%s"%(subdir[0],anDir[1]), "h_%s_bJets_%s"%(subdir[0],anDir[1]), 10,  -0.5, 9.5))
            jet1Pt.append(r.TH1D("h_%s_jet1Pt_%s"%(subdir[0],anDir[1]),"h_%s_jet1Pt_%s"%(subdir[0],anDir[1]),40,  0, 1000))
            jet2Pt.append(r.TH1D("h_%s_jet2Pt_%s"%(subdir[0],anDir[1]),"h_%s_jet2Pt_%s"%(subdir[0],anDir[1]),40,  0, 1000))
            ht    .append(r.TH1D("h_%s_ht_%s"%(subdir[0],anDir[1]),    "h_%s_ht_%s"%(subdir[0],anDir[1]),    30,  0, 3000))
            htminv.append(r.TH1D("h_%s_htminv_%s"%(subdir[0],anDir[1]),"h_%s_htminv_%s"%(subdir[0],anDir[1]),30,  0, 3000))
            met   .append(r.TH1D("h_%s_met_%s"%(subdir[0],anDir[1]),   "h_%s_met_%s"%(subdir[0],anDir[1]),   20,  0, 1000))
            mht   .append(r.TH1D("h_%s_mht_%s"%(subdir[0],anDir[1]),   "h_%s_mht_%s"%(subdir[0],anDir[1]),   20,  0, 1000))
            
            htvsmet    .append(r.TH2D("h_%s_htvsmet_%s"%(subdir[0],anDir[1]),    "h_%s_htvsmet_%s"%(subdir[0],anDir[1]),    30,  0, 3000, 20, 0, 1000))
            htminvvsmet.append(r.TH2D("h_%s_htminvvsmet_%s"%(subdir[0],anDir[1]),"h_%s_htminvvsmet_%s"%(subdir[0],anDir[1]),30,  0, 3000, 20, 0, 1000))
            htvsmht    .append(r.TH2D("h_%s_htvsmht_%s"%(subdir[0],anDir[1]),    "h_%s_htvsmht_%s"%(subdir[0],anDir[1]),    30,  0, 3000, 20, 0, 1000))
            htminvvsmht.append(r.TH2D("h_%s_htminvvsmht_%s"%(subdir[0],anDir[1]),"h_%s_htminvvsmht_%s"%(subdir[0],anDir[1]),30,  0, 3000, 20, 0, 1000))
            
            dPhi1 .append(r.TH1D("h_%s_dPhi1_%s"%(subdir[0],anDir[1]),"h_%s_dPhi1_%s"%(subdir[0],anDir[1]),50,  0, 3.2))
            dPhi2 .append(r.TH1D("h_%s_dPhi2_%s"%(subdir[0],anDir[1]),"h_%s_dPhi2_%s"%(subdir[0],anDir[1]),50,  0, 3.2))
            dPhi3 .append(r.TH1D("h_%s_dPhi3_%s"%(subdir[0],anDir[1]),"h_%s_dPhi3_%s"%(subdir[0],anDir[1]),50,  0, 3.2))
            
            genBosonPt .append(r.TH1D("h_%s_genBosonPt_%s"%(subdir[0],anDir[1]), "h_%s_genBosonPt_%s"%(subdir[0],anDir[1]),  40,   0, 1000))
            genBosonEta.append(r.TH1D("h_%s_genBosonEta_%s"%(subdir[0],anDir[1]),"h_%s_genBosonEta_%s"%(subdir[0],anDir[1]),2*50,  -5, 5))
            
            
    ####
    for ad,anDir in enumerate(analysisDirs):
        for sd,subdir in enumerate(dirs):
            print str(ad)+"*"+str(len(dirs))+"+"+str(sd)
            print nJets[ad*len(dirs)+sd]
            #(ad*sdsize)+sd
            nVertices[ad*len(dirs)+sd].Sumw2()
            nVerticesReWt[ad*len(dirs)+sd].Sumw2()
            nJets[ad*len(dirs)+sd].Sumw2()
            nJetsminv[ad*len(dirs)+sd].Sumw2()
            bJets[ad*len(dirs)+sd].Sumw2()
            jet1Pt[ad*len(dirs)+sd].Sumw2()
            jet2Pt[ad*len(dirs)+sd].Sumw2()
            ht[ad*len(dirs)+sd].Sumw2()
            htminv[ad*len(dirs)+sd].Sumw2()
            mht[ad*len(dirs)+sd].Sumw2()
            met[ad*len(dirs)+sd].Sumw2()
            htvsmht[ad*len(dirs)+sd].Sumw2()
            htminvvsmht[ad*len(dirs)+sd].Sumw2()
            htvsmet[ad*len(dirs)+sd].Sumw2()
            htminvvsmet[ad*len(dirs)+sd].Sumw2()
            dPhi1[ad*len(dirs)+sd].Sumw2()
            dPhi2[ad*len(dirs)+sd].Sumw2()
            dPhi3[ad*len(dirs)+sd].Sumw2()
            genBosonPt[ad*len(dirs)+sd].Sumw2()
            genBosonEta[ad*len(dirs)+sd].Sumw2()



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
        ra2_Vertices            = event.ra2_Vertices      
        ra2_nJetsPt50Eta25      = event.ra2_nJetsGenPt50Eta25 
        ra2_nJetsPt50Eta25MInv  = event.ra2_nJetsGenPt50Eta25MInv
        ra2_bJetsPt30Eta24  = event.ra2_bJetsGenPt30Eta24 
        ra2_PUWt            = event.ra2_PUWt           
        ra2_EventWt         = event.ra2_EventWt        
        ra2_HT              = event.ra2_genHT            
        ra2_HTMInv          = event.ra2_genHTMInv            
        ra2_MHT             = event.ra2_genMHT            
        ra2_MET             = event.ra2_genMET            
        ra2_dPhi1           = event.ra2_genDPhiMHT1
        ra2_dPhi2           = event.ra2_genDPhiMHT2
        ra2_dPhi3           = event.ra2_genDPhiMHT3
        ra2_Jet1Pt          = event.ra2_genJet1Pt         
        ra2_Jet2Pt          = event.ra2_genJet2Pt         
        ra2_genBosonPt      = event.ra2_genBoson1Pt         
        ra2_genBosonEta     = event.ra2_genBoson1Eta         

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
        for ad,anDir in enumerate(analysisDirs):                
            for sd,subdir in enumerate(dirs):
                outputFile.cd()
                outputFile.cd(anDir[1]+"/"+subdir[0])
                countJetsMInv = event.ra2_nJetsGenPt50Eta25MInv
                countJets = event.ra2_nJetsGenPt50Eta25

                #print "%s %s: #J:%d #JI:%d pT:%2.2f"%(anDir[1],subdir[0],countJets,countJetsMInv,ra2_genBosonPt)
                if (cutF.phenoCuts(event,ad,options.minpt,countJetsMInv,ra2_HTMInv,subdir[1])):
                    nJetsminv[ad*len(dirs)+sd]  .Fill(ra2_nJetsPt50Eta25MInv,ra2_EventWt)
                    htminv[ad*len(dirs)+sd]     .Fill(ra2_HTMInv            ,ra2_EventWt)
                    htminvvsmht[ad*len(dirs)+sd].Fill(ra2_HTMInv,ra2_MHT    ,ra2_EventWt)
                    htminvvsmet[ad*len(dirs)+sd].Fill(ra2_HTMInv,ra2_MET    ,ra2_EventWt)

                if (cutF.phenoCuts(event,ad,options.minpt,countJets,ra2_HT,subdir[1])):
                    nVertices[ad*len(dirs)+sd]    .Fill(ra2_Vertices,ra2_EventWt)
                    nVerticesReWt[ad*len(dirs)+sd].Fill(ra2_Vertices,ra2_EventWt*ra2_PUWt)
                    
                    nJets[ad*len(dirs)+sd] .Fill(ra2_nJetsPt50Eta25,ra2_EventWt)
                    bJets[ad*len(dirs)+sd] .Fill(ra2_bJetsPt30Eta24,ra2_EventWt)
                    jet1Pt[ad*len(dirs)+sd].Fill(ra2_Jet1Pt        ,ra2_EventWt)
                    jet2Pt[ad*len(dirs)+sd].Fill(ra2_Jet2Pt        ,ra2_EventWt)
                    ht[ad*len(dirs)+sd]    .Fill(ra2_HT            ,ra2_EventWt)
                    mht[ad*len(dirs)+sd]   .Fill(ra2_MHT           ,ra2_EventWt)
                    met[ad*len(dirs)+sd]   .Fill(ra2_MET           ,ra2_EventWt)
                    
                    htvsmht[ad*len(dirs)+sd].Fill(ra2_HT,ra2_MHT,ra2_EventWt)
                    htvsmet[ad*len(dirs)+sd].Fill(ra2_HT,ra2_MET,ra2_EventWt)
                    
                    dPhi1[ad*len(dirs)+sd].Fill(ra2_dPhi1,ra2_EventWt)
                    dPhi2[ad*len(dirs)+sd].Fill(ra2_dPhi2,ra2_EventWt)
                    dPhi3[ad*len(dirs)+sd].Fill(ra2_dPhi3,ra2_EventWt)
                    
                    genBosonPt[ad*len(dirs)+sd] .Fill(ra2_genBosonPt ,ra2_EventWt)
                    genBosonEta[ad*len(dirs)+sd].Fill(ra2_genBosonEta,ra2_EventWt)
        
                
        #########
        i = i + 1

    #####        

    print "50:%d 100:%d 150:%d 175:%d"%(countPT50,countPT100,countPT150,countPT175)

    outputFile.cd()
    
    ####
    for ad,anDir in enumerate(analysisDirs):
        for sd,subdir in enumerate(dirs):
            outputFile.cd()
            outputFile.cd(anDir[1]+"/"+subdir[0])
            nVertices[ad*len(dirs)+sd].Write()
            nVerticesReWt[ad*len(dirs)+sd].Write()
            nJets[ad*len(dirs)+sd].Write()
            nJetsminv[ad*len(dirs)+sd].Write()
            bJets[ad*len(dirs)+sd].Write()
            jet1Pt[ad*len(dirs)+sd].Write()
            jet2Pt[ad*len(dirs)+sd].Write()
            ht[ad*len(dirs)+sd].Write()
            htminv[ad*len(dirs)+sd].Write()
            mht[ad*len(dirs)+sd].Write()
            met[ad*len(dirs)+sd].Write()
            htvsmht[ad*len(dirs)+sd].Write()
            htminvvsmht[ad*len(dirs)+sd].Write()
            htvsmet[ad*len(dirs)+sd].Write()
            htminvvsmet[ad*len(dirs)+sd].Write()
            dPhi1[ad*len(dirs)+sd].Write()
            dPhi2[ad*len(dirs)+sd].Write()
            dPhi3[ad*len(dirs)+sd].Write()
            genBosonPt[ad*len(dirs)+sd].Write()
            genBosonEta[ad*len(dirs)+sd].Write()

    outputFile.Write()
    outputFile.Close()

##########################
if __name__ == '__main__':
    main()
    
