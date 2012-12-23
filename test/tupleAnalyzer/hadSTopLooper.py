import sys,os
import ROOT as r
from array import array
import math
import string
import re
import itertools
#import argparse
import optparse

import analysisCuts as cutF

def main() :
    parser = optparse.OptionParser(description="Switch for data/MC running")
    parser.add_option('-m', action="store_true", default=False, dest="isMC")
    parser.add_option('-d', action="store_true", default=False, dest="debug")
    parser.add_option('-b', action="store_true", default=False, dest="tightCSV")
    parser.add_option('-j', action="store",      default=5,     dest="numJets",type="int")
    parser.add_option('-f', action="store",      default=0,     dest="subsec", type="int")
    parser.add_option('-s', action="store",      default="zinv",dest="sample", type="string")
    parser.add_option('-t', action="store",      default="analysis",dest="treeName", type="string")
    options, args = parser.parse_args()

    r.gROOT.SetBatch(True)
    #debug = False
    myWorkingDir = os.getcwd()
    
    outFileName = "topTaggingPlots_%s_tight%d_min%d_job%d.root"%(options.sample,options.tightCSV,options.numJets,options.subsec)
    if options.debug:
        outFileName = "topTaggingPlots_%s_tight%d_min%d_test.root"%(options.sample,options.tightCSV,options.numJets)
        
    print outFileName
    outputFile = r.TFile(outFileName,"RECREATE")
    
    print ('%s%dLoose/RA2Values'%(options.treeName,options.numJets))
    chLoose   = r.TChain('%s%dLoose/RA2Values'%(options.treeName,options.numJets))
    print ('%s%dM/RA2Values'    %(options.treeName,options.numJets))
    chNominal = r.TChain('%s%dM/RA2Values'    %(options.treeName,options.numJets))

    sfCorr = 1.
    subfiles = [
        "*",#0
        "*_?"  ,"*_1?" ,"*_2?" ,"*_3?" ,"*_4?" ,"*_5?" ,"*_6?" ,"*_7?" ,"*_8?" ,"*_9?",#10
        "*_10?","*_11?","*_12?","*_13?","*_14?","*_15?","*_16?","*_17?","*_18?","*_19?",#20
        "*_20?","*_21?","*_22?","*_23?","*_24?","*_25?","*_26?","*_27?","*_28?","*_29?",#30
        "*_30?","*_31?","*_32?","*_33?","*_34?","*_35?","*_36?","*_37?","*_38?","*_39?",#40
        "*_40?","*_41?","*_42?","*_43?","*_44?","*_45?","*_46?","*_47?","*_48?","*_49?",#50
        "*_50?","*_51?","*_52?","*_53?","*_54?","*_55?","*_56?","*_57?","*_58?","*_59?",#60
        "*_60?","*_61?","*_62?","*_63?","*_64?","*_65?","*_66?","*_67?","*_68?","*_69?",#70
        "*_70?","*_71?","*_72?","*_73?","*_74?","*_75?","*_76?","*_77?","*_78?","*_79?",#80
        "*_80?","*_81?","*_82?","*_83?","*_84?","*_85?","*_86?","*_87?","*_88?","*_89?",#90
        "*_90?","*_91?","*_92?","*_93?","*_94?","*_95?","*_96?","*_97?","*_98?","*_99?",#100
        "*_100?","*_101?","*_102?","*_103?","*_104?","*_105?","*_106?","*_107?","*_108?","*_109?",#110
        "*_110?","*_111?","*_112?","*_113?","*_114?","*_115?","*_116?","*_117?","*_118?","*_119?",#120
        "*_120?","*_121?","*_122?","*_123?","*_124?","*_125?","*_126?","*_127?","*_128?","*_129?",#130
        ]
    
    if options.debug:
        if options.sample=="ttbar":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ttjets_tunez2star_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ttjets_tunez2star_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="zinv50":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht50to100_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht50to100_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            sfCorr = 4040980/24063998
        elif options.sample=="zinv100":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht100to200_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht100to200_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            sfCorr = 4416646./4416646.
        elif options.sample=="zinv200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht200to400_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht200to400_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            sfCorr = 5055885./9745619.
        elif options.sample=="zinv400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            sfCorr = 1006928./5095710.
        elif options.sample=="gjets200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht200to400_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht200to400_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            sfCorr = 10494617./10494617.
        elif options.sample=="gjets400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            sfCorr = 1611963./9539562.
        elif options.sample=="zmumu200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht200to400_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht200to400_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            sfCorr = 3789889./3789889.
        elif options.sample=="zmumu400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht400toinfv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht400toinfv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            sfCorr = 1703863./2727789.
            
        elif options.sample=="photon2012a":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="photon2012b":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="photon2012ab":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="photon2012c":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="photonall":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            
        elif options.sample=="muon2012a":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="muon2012b":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="muon2012ab":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="muon2012c":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="muonall":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            
    else:
        if options.sample=="ttbar":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ttjets_tunez2star_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ttjets_tunez2star_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv50":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht50to100_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht50to100_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 4040980/24063998
        elif options.sample=="zinv100":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht100to200_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht100to200_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 4416646./4416646.
        elif options.sample=="zinv200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht200to400_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht200to400_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 5055885./9745619.
        elif options.sample=="zinv400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 1006928./5095710.
        elif options.sample=="gjets200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht200to400_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht200to400_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 10494617./10494617.
        elif options.sample=="gjets400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 1611963./9539562.
        elif options.sample=="zmumu200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht200to400_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht200to400_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 3789889./3789889.
        elif options.sample=="zmumu400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht400toinfv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht400toinfv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 1703863./2727789.
            
        elif options.sample=="photon2012a":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="photon2012b":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="photon2012ab":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="photon2012c":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="photonall":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            
        elif options.sample=="muon2012a":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="muon2012b":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="muon2012ab":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="muon2012c":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="muonall":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chNominal.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            

    csNoMETDir        = outputFile .mkdir("csNoMET")
    csNoMETNJets      = csNoMETDir .mkdir("nJets")
    csNoMETLooseTT    = csNoMETDir .mkdir("looseTT")
    csNoMETNominalTT  = csNoMETDir .mkdir("nominalTT")
    
    csMETDir          = outputFile.mkdir("csMET")
    csMETNJets        = csMETDir  .mkdir("nJets")
    csMETLooseTT      = csMETDir  .mkdir("looseTT")
    csMETNominalTT    = csMETDir  .mkdir("nominalTT")

    invertedDir       = outputFile .mkdir("inverted")
    invertedNJets     = invertedDir.mkdir("nJets")
    invertedLooseTT   = invertedDir.mkdir("looseTT")
    invertedNominalTT = invertedDir.mkdir("nominalTT")

    baselineDir       = outputFile .mkdir("baseline")
    baselineNJets     = baselineDir.mkdir("nJets")
    baselineLooseTT   = baselineDir.mkdir("looseTT")
    baselineNominalTT = baselineDir.mkdir("nominalTT")

    analysisDirs = [
        [csNoMETDir,       "csNoMET" ], #0
        [csNoMETNJets,     "nJets"],    #1
        [csNoMETLooseTT,   "looseTT"],  #2
        [csNoMETNominalTT, "nominalTT"],#3
        [csMETDir,         "csMET" ],   #4
        [csMETNJets,       "nJets"],    #5
        [csMETLooseTT,     "looseTT"],  #6
        [csMETNominalTT,   "nominalTT"],#7
        [invertedDir,      "inverted" ],#8
        [invertedNJets,    "nJets"],    #9
        [invertedLooseTT,  "looseTT"],  #10
        [invertedNominalTT,"nominalTT"],#11
        [baselineDir,      "baseline" ],#12
        [baselineNJets,    "nJets"],    #13
        [baselineLooseTT,  "looseTT"],  #14
        [baselineNominalTT,"nominalTT"],#15
        ]
    
    nVtx      = []
    nJets     = []
    nJetsCSVM = []
    nJetsCSVT = []
    jet1Pt    = []
    jet2Pt    = []
    jet3Pt    = []
    jet4Pt    = []
    met       = []
    
    loose_bestTopJetMass = []
    loose_MTbestTopJet   = []
    loose_MTbJet         = []
    loose_MT2            = []
    
    nominal_bestTopJetMass = []
    nominal_MTbestTopJet   = []
    nominal_MTbJet         = []
    nominal_MT2            = []
    
    dPhi1    = []
    dPhi2    = []
    dPhi3    = []
    dPhi4    = []
    dPhiMinBCSVM = []
    dPhiMinBCSVT = []

    ###variables
    for d,dir in enumerate(analysisDirs):
        outputFile.cd()
        dir[0].cd()

        nVtx     .append(r.TH1D("h_nVtx",     "h_nVtx",     50,  -0.5, 49.5))
        nJets    .append(r.TH1D("h_nJets",    "h_nJets",    15,  -0.5, 14.5))
        nJetsCSVM.append(r.TH1D("h_nJetsCSVM","h_nJetsCSVM",10,  -0.5, 9.5))
        nJetsCSVT.append(r.TH1D("h_nJetsCSVT","h_nJetsCSVT",10,  -0.5, 9.5))
        jet1Pt   .append(r.TH1D("h_jet1Pt",   "h_jet1Pt",   2*50, 0,   1000))
        jet2Pt   .append(r.TH1D("h_jet2Pt",   "h_jet2Pt",   2*50, 0,   1000))
        jet3Pt   .append(r.TH1D("h_jet3Pt",   "h_jet3Pt",   2*50, 0,   1000))
        jet4Pt   .append(r.TH1D("h_jet4Pt",   "h_jet4Pt",   2*50, 0,   1000))
        met      .append(r.TH1D("h_met",      "h_met",      2*50, 0,   1000))
        
        loose_bestTopJetMass.append(r.TH1D("h_loose_bestTopJetMass", "h_loose_bestTopJetMass", 2*50, 0, 1000))
        loose_MTbestTopJet  .append(r.TH1D("h_loose_MTbestTopJet"  , "h_loose_MTbestTopJet"  , 2*75, 0, 1500))
        loose_MTbJet        .append(r.TH1D("h_loose_MTbJet"        , "h_loose_MTbJet"        , 2*75, 0, 1500))
        loose_MT2           .append(r.TH1D("h_loose_MT2"           , "h_loose_MT2"           , 2*50, 0, 1000))
    
        nominal_bestTopJetMass.append(r.TH1D("h_nominal_bestTopJetMass", "h_nominal_bestTopJetMass", 2*50, 0, 1000))
        nominal_MTbestTopJet  .append(r.TH1D("h_nominal_MTbestTopJet"  , "h_nominal_MTbestTopJet"  , 2*75, 0, 1500))
        nominal_MTbJet        .append(r.TH1D("h_nominal_MTbJet"        , "h_nominal_MTbJet"        , 2*75, 0, 1500))
        nominal_MT2           .append(r.TH1D("h_nominal_MT2"           , "h_nominal_MT2"           , 2*50, 0, 1000))
    
        dPhi1   .append(r.TH1D("h_dPhi1","h_dPhi1",50,  0, 3.2))
        dPhi2   .append(r.TH1D("h_dPhi2","h_dPhi2",50,  0, 3.2))
        dPhi3   .append(r.TH1D("h_dPhi3","h_dPhi3",50,  0, 3.2))
        dPhi4   .append(r.TH1D("h_dPhi4","h_dPhi4",50,  0, 3.2))
        dPhiMinBCSVM.append(r.TH1D("h_dPhiMinBCSVM","h_dPhiMinBCSVM",50, 0, 3.2))
        dPhiMinBCSVT.append(r.TH1D("h_dPhiMinBCSVT","h_dPhiMinBCSVT",50, 0, 3.2))

        
    ####
        nVtx[d]     .Sumw2()
        nJets[d]    .Sumw2()
        nJetsCSVM[d].Sumw2()
        nJetsCSVT[d].Sumw2()
        jet1Pt[d]   .Sumw2()
        jet2Pt[d]   .Sumw2()
        jet3Pt[d]   .Sumw2()
        jet4Pt[d]   .Sumw2()
        met[d]      .Sumw2()

        loose_bestTopJetMass[d].Sumw2()
        loose_MTbestTopJet[d]  .Sumw2()
        loose_MTbJet[d]        .Sumw2()
        loose_MT2[d]           .Sumw2()

        nominal_bestTopJetMass[d].Sumw2()
        nominal_MTbestTopJet[d]  .Sumw2()
        nominal_MTbJet[d]        .Sumw2()
        nominal_MT2[d]           .Sumw2()

        dPhi1[d]   .Sumw2()
        dPhi2[d]   .Sumw2()
        dPhi3[d]   .Sumw2()
        dPhi4[d]   .Sumw2()
        dPhiMinBCSVM[d].Sumw2()
        dPhiMinBCSVT[d].Sumw2()

    ##################
    fChain = chNominal

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
    for evLoose,evNominal in itertools.izip(chLoose,chNominal):
        # ==============print number of events done == == == == == == == =
        if ( i==0):
            tsw.Start()
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
                sys.stdout.write("t="+str(time))
                sys.stdout.write(" projected finish=%7d s("%(finTime))
                sys.stdout.write("%2.2f min).   "%(finMin))
                sys.stdout.write("\n")
                sys.stdout.flush()
                tenpcount = tenpcount + 1
        
        elif ( (i*100)/nentries == onepcount ) :
            sys.stdout.write('.')
            sys.stdout.flush()
            onepcount = onepcount + 1

        #sys.stdout.flush()
        #print "PU_et %2.2f"%(evNominal.ra2_PUWt)
        ra2_nVtx            = evNominal.ra2_Vertices 
        ra2_nJetsPt30Eta24  = evNominal.ra2_nJetsPt30Eta24 
        ra2_nJetsCSVM       = evNominal.ra2_nJetsCSVM
        ra2_nJetsCSVT       = evNominal.ra2_nJetsCSVT
        ra2_PUWt            = evNominal.ra2_PUWt           
        ra2_EventWt         = evNominal.ra2_EventWt*sfCorr
        ra2_MET             = evNominal.ra2_MET            
        ra2_dPhi1           = evNominal.ra2_dPhiMET1         
        ra2_dPhi2           = evNominal.ra2_dPhiMET2         
        ra2_dPhi3           = evNominal.ra2_dPhiMET3         
        ra2_dPhi4           = evNominal.ra2_dPhiMET4         
        ra2_dPhiMinBCSVM    = evNominal.ra2_dPhiMETMinBCSVM         
        ra2_dPhiMinBCSVT    = evNominal.ra2_dPhiMETMinBCSVT         
        ra2_Jet1Pt          = evNominal.ra2_Jet1Pt         
        ra2_Jet1Eta         = evNominal.ra2_Jet1Eta        
        ra2_Jet2Pt          = evNominal.ra2_Jet2Pt         
        ra2_Jet2Eta         = evNominal.ra2_Jet2Eta        
        ra2_Jet3Pt          = evNominal.ra2_Jet3Pt         
        ra2_Jet3Eta         = evNominal.ra2_Jet3Eta        
        ra2_Jet4Pt          = evNominal.ra2_Jet4Pt         
        ra2_Jet4Eta         = evNominal.ra2_Jet4Eta        

        passElectronVeto    = evNominal.ra2_passElVeto
        passIndirectTauVeto = evNominal.ra2_passTauVeto
        passTrackVeto       = evNominal.ra2_passIsoTrkVeto

        passLeptonVetos = passElectronVeto and passIndirectTauVeto and passTrackVeto
        if options.sample=="zmumu200" or options.sample=="zmumu400" or options.sample=="muon2012a" or options.sample=="muon2012b" or options.sample=="muon2012ab" or options.sample=="muon2012c" or options.sample=="muonall":
            passLeptonVetos = passElectronVeto and passIndirectTauVeto and passTrackVeto
        else:
            passLeptonVetos = passLeptonVetos and evNominal.ra2_passMuVeto
        ##top tagger variables 
        ra2_loose_bestTopJetMass  = evLoose.ra2_bestTopJetMass 
        ra2_loose_MTbestTopJet    = evLoose.ra2_TbestTopJet    
        ra2_loose_MTbJet          = evLoose.ra2_TbJet          
        ra2_loose_MT2             = evLoose.ra2_MT2            
        passLooseTopTagger = cutF.topTaggerCuts(evLoose)

        ra2_nominal_bestTopJetMass  = evNominal.ra2_bestTopJetMass 
        ra2_nominal_MTbestTopJet    = evNominal.ra2_TbestTopJet    
        ra2_nominal_MTbJet          = evNominal.ra2_TbJet          
        ra2_nominal_MT2             = evNominal.ra2_MT2            
        passNominalTopTagger = cutF.topTaggerCuts(evNominal)

        ####Fill plots
        for d,dir in enumerate(analysisDirs):
            outputFile.cd()
            dir[0].cd()
            fillPlots = False
            if d < 4 and cutF.hadStopControlSample(evNominal,False):
                ###Control sample no MET cut
                if d == 3 and cutF.topTaggerCuts(evNominal):
                    ###Nominal top tagging
                    fillPlots = True
                if d == 2 and cutF.topTaggerCuts(evLoose):
                    ###Loose top tagging
                    fillPlots = True
                if d == 1 and ra2_nJetsPt30Eta24 > (options.numJets-1):
                    ###nJets > minNJets
                    fillPlots = True
                if d == 0:
                    ###No top tagging
                    fillPlots = True
            elif d < 8 and cutF.hadStopControlSample(evNominal,True):
                ###Control sample with MET cut
                if d == 7 and cutF.topTaggerCuts(evNominal):
                    ###Nominal top tagging
                    fillPlots = True
                if d == 6 and cutF.topTaggerCuts(evLoose):
                    ###Loose top tagging
                    fillPlots = True
                if d == 5 and ra2_nJetsPt30Eta24 > (options.numJets-1):
                    ###nJets > minNJets
                    fillPlots = True
                if d == 4:
                    ###No top tagging
                    fillPlots = True

            elif d < 12 and cutF.hadStopBaseline(evNominal,1,options.tightCSV,2,True):
                ###Nominal selections
                if d == 11 and cutF.topTaggerCuts(evNominal):
                    ###Nominal top tagging
                    fillPlots = True
                if d == 10 and cutF.topTaggerCuts(evLoose):
                    ###Loose top tagging
                    fillPlots = True
                if d == 9 and ra2_nJetsPt30Eta24 > (options.numJets-1):
                    ###No top tagging
                    fillPlots = True
                if d == 8:
                    ###No top tagging
                    fillPlots = True

            elif d < 16 and cutF.hadStopBaseline(evNominal,1,options.tightCSV,2,False):
                ###Nominal selections
                if d == 15 and cutF.topTaggerCuts(evNominal):
                    ###Nominal top tagging
                    fillPlots = True
                if d == 14 and cutF.topTaggerCuts(evLoose):
                    ###Loose top tagging
                    fillPlots = True
                if d == 13 and ra2_nJetsPt30Eta24 > (options.numJets-1):
                    ###nJets > minNJets
                    fillPlots = True
                if d == 12:
                    ###No top tagging
                    fillPlots = True

            if fillPlots and passLeptonVetos:
                nVtx[d]     .Fill(ra2_nVtx          ,ra2_EventWt)
                nJets[d]    .Fill(ra2_nJetsPt30Eta24,ra2_EventWt)
                nJetsCSVM[d].Fill(ra2_nJetsCSVM     ,ra2_EventWt)
                nJetsCSVT[d].Fill(ra2_nJetsCSVT     ,ra2_EventWt)
                jet1Pt[d]   .Fill(ra2_Jet1Pt        ,ra2_EventWt)
                jet2Pt[d]   .Fill(ra2_Jet2Pt        ,ra2_EventWt)
                jet3Pt[d]   .Fill(ra2_Jet3Pt        ,ra2_EventWt)
                jet4Pt[d]   .Fill(ra2_Jet4Pt        ,ra2_EventWt)
                met[d]      .Fill(ra2_MET           ,ra2_EventWt)
                
                loose_bestTopJetMass[d].Fill(ra2_loose_bestTopJetMass,ra2_EventWt)
                loose_MTbestTopJet[d]  .Fill(ra2_loose_MTbestTopJet  ,ra2_EventWt)
                loose_MTbJet[d]        .Fill(ra2_loose_MTbJet        ,ra2_EventWt)
                loose_MT2[d]           .Fill(ra2_loose_MT2           ,ra2_EventWt)
                
                nominal_bestTopJetMass[d].Fill(ra2_nominal_bestTopJetMass,ra2_EventWt)
                nominal_MTbestTopJet[d]  .Fill(ra2_nominal_MTbestTopJet  ,ra2_EventWt)
                nominal_MTbJet[d]        .Fill(ra2_nominal_MTbJet        ,ra2_EventWt)
                nominal_MT2[d]           .Fill(ra2_nominal_MT2           ,ra2_EventWt)
                
                dPhi1[d]   .Fill(ra2_dPhi1   ,ra2_EventWt)
                dPhi2[d]   .Fill(ra2_dPhi2   ,ra2_EventWt)
                dPhi3[d]   .Fill(ra2_dPhi3   ,ra2_EventWt)
                dPhi4[d]   .Fill(ra2_dPhi4   ,ra2_EventWt)
                dPhiMinBCSVM[d].Fill(ra2_dPhiMinBCSVM,ra2_EventWt)
                dPhiMinBCSVT[d].Fill(ra2_dPhiMinBCSVT,ra2_EventWt)
                
        ### at least one b-jet
                
        #########
        i = i + 1

    #####        
    outputFile.cd()

    ####
    for d,dir in enumerate(analysisDirs):
        outputFile.cd()
        dir[0].cd()

        nVtx[d]     .Write()
        nJets[d]    .Write()
        nJetsCSVM[d].Write()
        nJetsCSVT[d].Write()
        jet1Pt[d]   .Write()
        jet2Pt[d]   .Write()
        jet3Pt[d]   .Write()
        jet4Pt[d]   .Write()
        met[d]      .Write()

        loose_bestTopJetMass[d].Write()
        loose_MTbestTopJet[d]  .Write()
        loose_MTbJet[d]        .Write()
        loose_MT2[d]           .Write()

        nominal_bestTopJetMass[d].Write()
        nominal_MTbestTopJet[d]  .Write()
        nominal_MTbJet[d]        .Write()
        nominal_MT2[d]           .Write()

        dPhi1[d]   .Write()
        dPhi2[d]   .Write()
        dPhi3[d]   .Write()
        dPhi4[d]   .Write()
        dPhiMinBCSVM[d].Write()
        dPhiMinBCSVT[d].Write()


    outputFile.Write()
    outputFile.Close()

##########################
if __name__ == '__main__':
    main()
    

#  LocalWords:  elif
