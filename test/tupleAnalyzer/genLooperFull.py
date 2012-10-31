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
    parser.add_option('-s', action="store",      default="gjets",dest="sample", type="string")
    parser.add_option('-f', action="store",      default="0",    dest="subsec", type="int")
    parser.add_option('-p', action="store",      default=100.0,dest="minpt",    type="float")
    parser.add_option('-r', action="store",      default=0.5,dest="cutDR",      type="float")
    options, args = parser.parse_args()

    r.gROOT.SetBatch(True)

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
    outFileName = "genLevelPlots_%s_%d.root"%(options.sample,options.subsec)
    if options.debug:
        outFileName = "genLevelPlots_%s_test.root"%(options.sample)

    print outFileName
    sys.stdout.flush()
    outputFile = r.TFile(outFileName,"RECREATE")
    precutsDir    = outputFile.mkdir("precuts")
    bosonptcutDir = outputFile.mkdir("bosonptcut")
    acceptanceDir = outputFile.mkdir("acceptance")
    recoidDir     = outputFile.mkdir("recoid")
    recopfisoDir  = outputFile.mkdir("recopfiso")

    analysisDirs = [
        [precutsDir   ,"precuts"],
        [bosonptcutDir,"bosonptcut"],
        [acceptanceDir,"acceptance"],
        [recoidDir    ,"recoid"],
        [recopfisoDir ,"recopfiso"],

        ]
    for anDir in analysisDirs:
        inclusiveDir  = anDir[0].mkdir("inclusive")
        inclusivesDir = anDir[0].mkdir("inclusive3")
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
        ["inclusive" ,[2,100]],
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

    idChain    = r.TChain('directPhotonsID/RA2Values')
    pfisoChain = r.TChain('directPhotonsIDPFIso/RA2Values')

    if options.debug:
        if options.sample=="gjets":
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v9/res/*_1_?_???.root")
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v9/res/*_1_?_???.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v9/res/*_1_?_???.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v9/res/*_1_?_???.root")
        elif options.sample=="zinv":
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v9/res/*_?_?_???.root")
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v9/res/*_?_?_???.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v9/res/*_?_?_???.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v9/res/*_?_?_???.root")
    else:
        if options.sample=="gjets":
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v9/res/%s.root"%(subfiles[options.subsec]))
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v9/res/%s.root"%(subfiles[options.subsec]))
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht200_gen_tree_v9/res/%s.root"%(subfiles[options.subsec]))
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/gjetsht400_gen_tree_v9/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv":
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v9/res/*.root")
            idChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v9/res/*.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht200_gen_tree_v9/res/*.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw525p1/genStudy/zinvht400_gen_tree_v9/res/*.root")
    
    ###variables
    nJets   = []
    jet1Pt  = []
    jet2Pt  = []
    ht      = []
    mht     = []
    #htvsmht = []
    genBosonPt  = []
    genBosonEta = []
    genBosonDR  = []
    bosonPt     = []
    bosonEta    = []
    bosonDR     = []

    nJetsminv   = []
    jet1Ptminv  = []
    jet2Ptminv  = []
    htminv      = []
    mhtminv     = []
    #htvsmhtminv = []
    genBosonPtminv  = []
    genBosonEtaminv = []
    genBosonDRminv  = []
    bosonPtminv     = []
    bosonEtaminv    = []
    bosonDRminv     = []

###cut on DR > cutDR
    nJetsCutDR   = []
    jet1PtCutDR  = []
    jet2PtCutDR  = []
    htCutDR      = []
    mhtCutDR     = []
    #htvsmhtCutDR = []
    genBosonPtCutDR  = []
    genBosonEtaCutDR = []
    genBosonDRCutDR  = []
    bosonPtCutDR     = []
    bosonEtaCutDR    = []
    bosonDRCutDR     = []

    nJetsminvCutDR   = []
    jet1PtminvCutDR  = []
    jet2PtminvCutDR  = []
    htminvCutDR      = []
    mhtminvCutDR     = []
    #htvsmhtminvCutDR = []
    genBosonPtminvCutDR  = []
    genBosonEtaminvCutDR = []
    genBosonDRminvCutDR  = []
    bosonPtminvCutDR     = []
    bosonEtaminvCutDR    = []
    bosonDRminvCutDR     = []

    for ad,anDir in enumerate(analysisDirs):
        for j,jetdir in enumerate(njbins):
            for h,htdir in enumerate(htbins):
                for m,mhtdir in enumerate(mhtbins):
                    outputFile.cd()
                    outputFile.cd(anDir[1]+"/"+jetdir[0])
                    workingDir = outputFile.Get(anDir[1]+"/"+jetdir[0])
                    workingDir.mkdir("ht"+htdir[0]+"_mht"+mhtdir[0])
                    outputFile.cd(anDir[1]+"/"+jetdir[0]+"/"+"ht"+htdir[0]+"_mht"+mhtdir[0])
                    
                    nJets .append(r.TH1D("h_%s_ht%s_mht%s_nJets_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_nJets_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 15, -0.5, 14.5))
                    jet1Pt.append(r.TH1D("h_%s_ht%s_mht%s_jet1Pt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet1Pt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),40,  0, 1000))
                    jet2Pt.append(r.TH1D("h_%s_ht%s_mht%s_jet2Pt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet2Pt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),40,  0, 1000))
                    ht    .append(r.TH1D("h_%s_ht%s_mht%s_ht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    "h_%s_ht%s_mht%s_ht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    30,  0, 3000))
                    mht   .append(r.TH1D("h_%s_ht%s_mht%s_mht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   "h_%s_ht%s_mht%s_mht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   20,  0, 1000))
                    #htvsmht.append(r.TH2D("h_%s_ht%s_mht%s_htvsmht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_htvsmht_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),30,  0, 3000, 20, 0, 1000))
                    genBosonPt .append(r.TH1D("h_%s_ht%s_mht%s_genBosonPt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_genBosonPt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  40,   0, 1000))
                    genBosonEta.append(r.TH1D("h_%s_ht%s_mht%s_genBosonEta_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_genBosonEta_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  -5, 5))
                    genBosonDR .append(r.TH1D("h_%s_ht%s_mht%s_genBosonDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_genBosonDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 2.5))
                    bosonPt    .append(r.TH1D("h_%s_ht%s_mht%s_bosonPt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_bosonPt_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  40,   0, 1000))
                    bosonEta   .append(r.TH1D("h_%s_ht%s_mht%s_bosonEta_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_bosonEta_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  -5, 5))
                    bosonDR    .append(r.TH1D("h_%s_ht%s_mht%s_bosonDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_bosonDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 2.5))
                    
                    nJetsminv .append(r.TH1D("h_%s_ht%s_mht%s_nJetsminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_nJetsminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 15, -0.5, 14.5))
                    jet1Ptminv.append(r.TH1D("h_%s_ht%s_mht%s_jet1Ptminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet1Ptminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),40,  0, 1000))
                    jet2Ptminv.append(r.TH1D("h_%s_ht%s_mht%s_jet2Ptminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet2Ptminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),40,  0, 500))
                    htminv    .append(r.TH1D("h_%s_ht%s_mht%s_htminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_htminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1])            ,30,  0, 3000))
                    mhtminv   .append(r.TH1D("h_%s_ht%s_mht%s_mhtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_mhtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1])         ,20,  0, 1000))
                    #htvsmhtminv.append(r.TH2D("h_%s_ht%s_mht%s_htvsmhtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_htvsmhtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),30,  0, 3000, 20, 0, 1000))
                    genBosonPtminv .append(r.TH1D("h_%s_ht%s_mht%s_genBosonPtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_genBosonPtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  40,   0, 1000))
                    genBosonEtaminv.append(r.TH1D("h_%s_ht%s_mht%s_genBosonEtaminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_genBosonEtaminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  -5, 5))
                    genBosonDRminv .append(r.TH1D("h_%s_ht%s_mht%s_genBosonDRminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_genBosonDRminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 2.5))
                    bosonPtminv    .append(r.TH1D("h_%s_ht%s_mht%s_bosonPtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_bosonPtminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  40,   0, 1000))
                    bosonEtaminv   .append(r.TH1D("h_%s_ht%s_mht%s_bosonEtaminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_bosonEtaminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  -5, 5))
                    bosonDRminv    .append(r.TH1D("h_%s_ht%s_mht%s_bosonDRminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_bosonDRminv_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 2.5))

                    ###CutDR
                    nJetsCutDR .append(r.TH1D("h_%s_ht%s_mht%s_nJetsCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_nJetsCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 15, -0.5, 14.5))
                    jet1PtCutDR.append(r.TH1D("h_%s_ht%s_mht%s_jet1PtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet1PtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),40,  0, 1000))
                    jet2PtCutDR.append(r.TH1D("h_%s_ht%s_mht%s_jet2PtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet2PtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),40,  0, 1000))
                    htCutDR    .append(r.TH1D("h_%s_ht%s_mht%s_htCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    "h_%s_ht%s_mht%s_htCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    30,  0, 3000))
                    mhtCutDR   .append(r.TH1D("h_%s_ht%s_mht%s_mhtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   "h_%s_ht%s_mht%s_mhtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   20,  0, 1000))
                    #htvsmhtCutDR.append(r.TH2D("h_%s_ht%s_mht%s_htvsmhtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_htvsmhtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),30,  0, 3000, 20, 0, 1000))
                    genBosonPtCutDR .append(r.TH1D("h_%s_ht%s_mht%s_genBosonPtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_genBosonPtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  40,   0, 1000))
                    genBosonEtaCutDR.append(r.TH1D("h_%s_ht%s_mht%s_genBosonEtaCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_genBosonEtaCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  -5, 5))
                    genBosonDRCutDR .append(r.TH1D("h_%s_ht%s_mht%s_genBosonDRCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_genBosonDRCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 2.5))
                    bosonPtCutDR    .append(r.TH1D("h_%s_ht%s_mht%s_bosonPtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_bosonPtCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  40,   0, 1000))
                    bosonEtaCutDR   .append(r.TH1D("h_%s_ht%s_mht%s_bosonEtaCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_bosonEtaCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  -5, 5))
                    bosonDRCutDR    .append(r.TH1D("h_%s_ht%s_mht%s_bosonDRCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_bosonDRCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 2.5))
                    
                    nJetsminvCutDR .append(r.TH1D("h_%s_ht%s_mht%s_nJetsminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_nJetsminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), 15, -0.5, 14.5))
                    jet1PtminvCutDR.append(r.TH1D("h_%s_ht%s_mht%s_jet1PtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet1PtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),40,  0, 1000))
                    jet2PtminvCutDR.append(r.TH1D("h_%s_ht%s_mht%s_jet2PtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_jet2PtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),40,  0, 1000))
                    htminvCutDR    .append(r.TH1D("h_%s_ht%s_mht%s_htminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    "h_%s_ht%s_mht%s_htminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),    30,  0, 3000))
                    mhtminvCutDR   .append(r.TH1D("h_%s_ht%s_mht%s_mhtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   "h_%s_ht%s_mht%s_mhtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),   20,  0, 1000))
                    #htvsmhtminvCutDR.append(r.TH2D("h_%s_ht%s_mht%s_htvsmhtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_htvsmhtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),30,  0, 3000, 20, 0, 1000))
                    genBosonPtminvCutDR .append(r.TH1D("h_%s_ht%s_mht%s_genBosonPtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_genBosonPtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  40,   0, 1000))
                    genBosonEtaminvCutDR.append(r.TH1D("h_%s_ht%s_mht%s_genBosonEtaminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_genBosonEtaminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  -5, 5))
                    genBosonDRminvCutDR .append(r.TH1D("h_%s_ht%s_mht%s_genBosonDRminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_genBosonDRminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 2.5))
                    bosonPtminvCutDR    .append(r.TH1D("h_%s_ht%s_mht%s_bosonPtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]), "h_%s_ht%s_mht%s_bosonPtminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),  40,   0, 1000))
                    bosonEtaminvCutDR   .append(r.TH1D("h_%s_ht%s_mht%s_bosonEtaminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_bosonEtaminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  -5, 5))
                    bosonDRminvCutDR    .append(r.TH1D("h_%s_ht%s_mht%s_bosonDRminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),"h_%s_ht%s_mht%s_bosonDRminvCutDR_%s"%(jetdir[0],htdir[0],mhtdir[0],anDir[1]),50,  0, 2.5))
                    
    ####
                    nJets[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet1Pt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet2Pt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    ht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    mht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    #htvsmht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonPt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonEta[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonPt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonEta[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    
                    nJetsminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet1Ptminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet2Ptminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    htminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    mhtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    #htvsmhtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonPtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonEtaminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonDRminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonPtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonEtaminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonDRminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()

                    nJetsCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet1PtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet2PtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    htCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    mhtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    #htvsmhtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonPtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonEtaCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonDRCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonPtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonEtaCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonDRCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()

                    nJetsminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet1PtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    jet2PtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    htminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    mhtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    #htvsmhtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonPtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonEtaminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    genBosonDRminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonPtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonEtaminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()
                    bosonDRminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Sumw2()



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
        ra2_PUWt               = idevent.ra2_PUWt           
        ra2_EventWt            = idevent.ra2_EventWt        
        ra2_HT                 = idevent.ra2_genHT            
        ra2_HTMInv             = idevent.ra2_genHTMInv            
        ra2_MHT                = idevent.ra2_genMHT            
        ra2_Jet1Pt             = idevent.ra2_genJet1Pt         
        ra2_Jet2Pt             = idevent.ra2_genJet2Pt         
        ra2_genBosonPt         = idevent.ra2_genBoson1Pt         
        ra2_genBosonEta        = idevent.ra2_genBoson1Eta         
        ra2_genBosonDR         = idevent.ra2_genBoson1MinDR         
        ra2_bosonIDPt          = idevent.ra2_boson1Pt         
        ra2_bosonIDEta         = idevent.ra2_boson1Eta         
        ra2_bosonIDDR          = idevent.ra2_boson1MinDR         
        ra2_bosonIDPFIsoPt     = pfisoevent.ra2_boson1Pt         
        ra2_bosonIDPFIsoEta    = pfisoevent.ra2_boson1Eta         
        ra2_bosonIDPFIsoDR     = pfisoevent.ra2_boson1MinDR         

        ra2_bosonPt  = 0
        ra2_bosonEta = 0
        ra2_bosonDR  = 0

        ####Fill plots
        for ad,anDir in enumerate(analysisDirs):
            if ad > 3:
                ra2_bosonPt  = ra2_bosonIDPFIsoPt
                ra2_bosonEta = ra2_bosonIDPFIsoEta
                ra2_bosonDR  = ra2_bosonIDPFIsoDR 
                recoPt       = ra2_bosonIDPFIsoPt
                passReco     = pfisoevent.ra2_genPassRecoIso
            else:
                ra2_bosonPt  = ra2_bosonIDPt
                ra2_bosonEta = ra2_bosonIDEta
                ra2_bosonDR  = ra2_bosonIDDR 
                recoPt       = ra2_bosonIDPt
                passReco     = idevent.ra2_genPassRecoIso
                
            for j,jetdir in enumerate(njbins):
                if cutF.ra2JetCuts(ra2_nJetsPt50Eta25,jetdir[0],jetdir[1]) or cutF.ra2JetCuts(ra2_nJetsPt50Eta25MInv,jetdir[0],jetdir[1]):
                    outputFile.cd()
                    outputFile.cd(anDir[1]+"/"+jetdir[0])
                    outputFile.cd(anDir[1]+"/"+jetdir[0]+"/"+"ht"+htdir[0]+"_mht"+mhtdir[0])
                    for m,mhtdir in enumerate(mhtbins):
                        if cutF.ra2MHTCuts(ra2_MHT,mhtdir[0],mhtdir[1]):
                            for h,htdir in enumerate(htbins):
                                if cutF.ra2HTCuts(HT,htdir[0],htdir[1]) or cutF.ra2HTCuts(HTMInv,htdir[0],htdir[1]):
                                    if (cutF.genCutsFull(idevent,ad,options.minpt,recoPt,passReco,ra2_nJetsPt50Eta25,jetdir[1],ra2_HT,htdir[1],mhtdir[1],-1)):
                                        nJets[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]  .Fill(ra2_nJetsPt50Eta25,ra2_EventWt)
                                        jet1Pt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_Jet1Pt        ,ra2_EventWt)
                                        jet2Pt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_Jet2Pt        ,ra2_EventWt)
                                        ht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]     .Fill(ra2_HT            ,ra2_EventWt)
                                        mht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_MHT           ,ra2_EventWt)
                                        # htvsmht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_HT,ra2_MHT    ,ra2_EventWt)
                                        
                                        genBosonPt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_genBosonPt ,ra2_EventWt)
                                        genBosonEta[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_genBosonEta,ra2_EventWt)
                                        genBosonDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_genBosonDR ,ra2_EventWt)
                                        bosonPt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_bosonPt    ,ra2_EventWt)
                                        bosonEta[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]   .Fill(ra2_bosonEta   ,ra2_EventWt)
                                        bosonDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_bosonDR    ,ra2_EventWt)
                                        
                                    if (cutF.genCutsFull(idevent,ad,options.minpt,recoPt,passReco,ra2_nJetsPt50Eta25MInv,jetdir[1],ra2_HTMInv,htdir[1],mhtdir[1],-1)):
                                        nJetsminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]  .Fill(ra2_nJetsPt50Eta25MInv,ra2_EventWt)
                                        jet1Ptminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_Jet1Pt        ,ra2_EventWt)
                                        jet2Ptminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_Jet2Pt        ,ra2_EventWt)
                                        htminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]     .Fill(ra2_HTMInv        ,ra2_EventWt)
                                        mhtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_MHT           ,ra2_EventWt)
                                        # htvsmhtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_HTMInv,ra2_MHT,ra2_EventWt)
                                        
                                        genBosonPtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_genBosonPt ,ra2_EventWt)
                                        genBosonEtaminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_genBosonEta,ra2_EventWt)
                                        genBosonDRminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_genBosonDR ,ra2_EventWt)
                                        bosonPtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_bosonPt    ,ra2_EventWt)
                                        bosonEtaminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]   .Fill(ra2_bosonEta   ,ra2_EventWt)
                                        bosonDRminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_bosonDR    ,ra2_EventWt)
                                        
                                        
                                    if (ra2_genBosonDR > options.cutDR):
                                        if (cutF.genCutsFull(idevent,ad,options.minpt,recoPt,passReco,ra2_nJetsPt50Eta25MInv,jetdir[1],ra2_HTMInv,htdir[1],mhtdir[1],options.cutDR)):
                                            nJetsminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]  .Fill(ra2_nJetsPt50Eta25MInv,ra2_EventWt)
                                            jet1Ptminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_Jet1Pt        ,ra2_EventWt)
                                            jet2Ptminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_Jet2Pt        ,ra2_EventWt)
                                            htminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]     .Fill(ra2_HTMInv        ,ra2_EventWt)
                                            mhtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_MHT           ,ra2_EventWt)
                                            # htvsmhtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_HTMInv,ra2_MHT,ra2_EventWt)
                                            
                                            genBosonPtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_genBosonPt ,ra2_EventWt)
                                            genBosonEtaminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_genBosonEta,ra2_EventWt)
                                            genBosonDRminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_genBosonDR ,ra2_EventWt)
                                            bosonPtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_bosonPt    ,ra2_EventWt)
                                            bosonEtaminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]   .Fill(ra2_bosonEta   ,ra2_EventWt)
                                            bosonDRminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_bosonDR    ,ra2_EventWt)
                                            
                                            
                                        if (cutF.genCutsFull(idevent,ad,options.minpt,recoPt,passReco,ra2_nJetsPt50Eta25MInv,jetdir[1],ra2_HTMInv,htdir[1],mhtdir[1],options.cutDR)):
                                            nJetsminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]  .Fill(ra2_nJetsPt50Eta25MInv,ra2_EventWt)
                                            jet1PtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_Jet1Pt        ,ra2_EventWt)
                                            jet2PtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_Jet2Pt        ,ra2_EventWt)
                                            htminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]     .Fill(ra2_HTMInv        ,ra2_EventWt)
                                            mhtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_MHT           ,ra2_EventWt)
                                            # htvsmhtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_HTMInv,ra2_MHT,ra2_EventWt)
                                            
                                            genBosonPtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_genBosonPt ,ra2_EventWt)
                                            genBosonEtaminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Fill(ra2_genBosonEta,ra2_EventWt)
                                            genBosonDRminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m] .Fill(ra2_genBosonDR ,ra2_EventWt)
                                            bosonPtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_bosonPt    ,ra2_EventWt)
                                            bosonEtaminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]   .Fill(ra2_bosonEta   ,ra2_EventWt)
                                            bosonDRminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m]    .Fill(ra2_bosonDR    ,ra2_EventWt)
                                            
                
        #########
        i = i + 1

    #####        
    outputFile.cd()

    ####
    for ad,anDir in enumerate(analysisDirs):
        for j,jetdir in enumerate(njbins):
            for h,htdir in enumerate(htbins):
                for m,mhtdir in enumerate(mhtbins):
                    #outputFile.cd()
                    #outputFile.cd(anDir[1]+"/"+jetdir[0])
                    outputFile.cd(anDir[1]+"/"+jetdir[0]+"/"+"ht"+htdir[0]+"_mht"+mhtdir[0])
                    print("saving histograms in::"+anDir[1]+"/"+jetdir[0]+"/"+"ht"+htdir[0]+"_mht"+mhtdir[0])
                    sys.stdout.flush()
                    
                    nJets[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    jet1Pt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    jet2Pt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    ht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    mht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    # htvsmht[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonPt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonEta[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonPt[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonEta[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    
                    nJetsminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    jet1Ptminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    jet2Ptminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    htminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    mhtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    # htvsmhtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonPtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonEtaminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonDRminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonPtminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonEtaminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonDRminv[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()

                    nJetsCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    jet1PtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    jet2PtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    htCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    mhtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    # htvsmhtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonPtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonEtaCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonDRCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonPtCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonEtaCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonDRCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()

                    nJetsminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    jet1PtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    jet2PtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    htminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    mhtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    # htvsmhtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonPtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonEtaminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    genBosonDRminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonPtminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonEtaminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()
                    bosonDRminvCutDR[(ad*len(njbins)*len(htbins)*len(mhtbins))+(j*len(htbins)*len(mhtbins))+(h*len(mhtbins))+m].Write()

    print "writing output file"
    sys.stdout.flush()
    outputFile.Write()
    ##print "closing output file"
    ##sys.stdout.flush()
    ##outputFile.Close()
    print "done!"
    sys.stdout.flush()

##########################
if __name__ == '__main__':
    main()
    
