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
    
    outFileName = "hadSTopTree_%s_job%d.root"%(options.sample,options.subsec)
    if options.debug:
        outFileName = "hadSTopTree_%s_test.root"%(options.sample)
        
    print outFileName
    outputFile = r.TFile(outFileName,"RECREATE")
    tree     = r.TTree( 'tree', 'tree for sample ' )
    
    print ('%s/RA2Values'%(options.treeName))
    chLoose   = r.TChain('%s/RA2Values'%(options.treeName))

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
        elif options.sample=="zinv50":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht50to100_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="zinv100":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht100to200_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="zinv200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht200to400_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="zinv400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="gjets200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht200to400_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="gjets400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="zmumu200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht200to400_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
        elif options.sample=="zmumu400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht400toinfv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            
        elif options.sample=="photon2012a":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

        elif options.sample=="photon2012b":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

        elif options.sample=="photon2012ab":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

        elif options.sample=="photon2012c":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

        elif options.sample=="photonall":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

        elif options.sample=="muon2012a":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

        elif options.sample=="muon2012b":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

        elif options.sample=="muon2012ab":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

        elif options.sample=="muon2012c":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")

        elif options.sample=="muonall":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/*_sTop_?.root")


    else:
        if options.sample=="ttbar":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ttjets_tunez2star_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv50":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht50to100_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 2844517./24063998.
        elif options.sample=="zinv100":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht100to200_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 2574000./4416646.
        elif options.sample=="zinv200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht200to400_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 3837885./9745619.
        elif options.sample=="zinv400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zinvjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 982928./5095710.#fixed
        elif options.sample=="gjets200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht200to400_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 7188617./10494617.
        elif options.sample=="gjets400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/gjetsht400toinf_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 1599963./9539562.
        elif options.sample=="zmumu200":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht200to400_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 3599062./3789889.
        elif options.sample=="zmumu400":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/zlljetsht400toinfv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            sfCorr = 1673863./2727789.
            
        elif options.sample=="photon2012a":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

        elif options.sample=="photon2012b":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

        elif options.sample=="photon2012ab":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

        elif options.sample=="photon2012c":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

        elif options.sample=="photonall":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/photondata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/singlephotondata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

        elif options.sample=="muon2012a":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

        elif options.sample=="muon2012b":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

        elif options.sample=="muon2012ab":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

        elif options.sample=="muon2012c":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

        elif options.sample=="muonall":
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Av1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Arecoverv1_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Bv4_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))
            chLoose.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/doublemudata2012Cv2_reco_tree_topTagged_09Nov/condor_output/%s.root"%(subfiles[options.subsec]))

            

    run         = array( 'i', [ 0 ] )
    event       = array( 'i', [ 0 ] )
    lumi        = array( 'i', [ 0 ] )
    nVtx        = array( 'i', [ 0 ] )
    nJets       = array( 'i', [ 0 ] )
    nJetsCSVM   = array( 'i', [ 0 ] )
    nJetsCSVT   = array( 'i', [ 0 ] )
    metVal      = array( 'd', [ 0. ] )
    dphi1       = array( 'd', [ 0. ] )
    dphi2       = array( 'd', [ 0. ] )
    dphi3       = array( 'd', [ 0. ] )
    dphi4       = array( 'd', [ 0. ] )
    dphiMinCSVM = array( 'd', [ 0. ] )
    dphiMinCSVT = array( 'd', [ 0. ] )

    ra2_nJets  = array( 'i', [ 0 ] )
    ra2_htVal  = array( 'd', [ 0. ] )
    ra2_mhtVal = array( 'd', [ 0. ] )
    ra2_dphi1  = array( 'd', [ 0. ] )
    ra2_dphi2  = array( 'd', [ 0. ] )
    ra2_dphi3  = array( 'd', [ 0. ] )

    bosonPt     = array( 'd', [ 0. ] )
    
    photonPt     = array( 'd', [ 0. ] )
    photonEta    = array( 'd', [ 0. ] )
    photonMinDR  = array( 'd', [ 0. ] )
    photonJet1DR = array( 'd', [ 0. ] )

    muon1Pt     = array( 'd', [ 0. ] )
    muon1Eta    = array( 'd', [ 0. ] )
    muon1MinDR  = array( 'd', [ 0. ] )
    muon1Jet1DR = array( 'd', [ 0. ] )
    muon2Pt     = array( 'd', [ 0. ] )
    muon2Eta    = array( 'd', [ 0. ] )
    muon2MinDR  = array( 'd', [ 0. ] )
    muon2Jet1DR = array( 'd', [ 0. ] )
    dimuonPt    = array( 'd', [ 0. ] )
    dimuonEta   = array( 'd', [ 0. ] )
    dimuonMinDR = array( 'd', [ 0. ] )
    dimuonMinDR = array( 'd', [ 0. ] )
    dimuonJet1DR= array( 'd', [ 0. ] )
    dimuonM     = array( 'd', [ 0. ] )

    jet1Pt    = array( 'd', [ 0. ] )
    jet1Eta   = array( 'd', [ 0. ] )
    jet2Pt    = array( 'd', [ 0. ] )
    jet2Eta   = array( 'd', [ 0. ] )
    jet3Pt    = array( 'd', [ 0. ] )
    jet3Eta   = array( 'd', [ 0. ] )
    jet4Pt    = array( 'd', [ 0. ] )
    jet4Eta   = array( 'd', [ 0. ] )

    csvmet1Pt    = array( 'd', [ 0. ] )
    csvmet1Eta   = array( 'd', [ 0. ] )
    csvmet2Pt    = array( 'd', [ 0. ] )
    csvmet2Eta   = array( 'd', [ 0. ] )

    csvtet1Pt    = array( 'd', [ 0. ] )
    csvtet1Eta   = array( 'd', [ 0. ] )
    csvtet2Pt    = array( 'd', [ 0. ] )
    csvtet2Eta   = array( 'd', [ 0. ] )

    loose_bestTopJetMass  = array( 'd', [ 0. ] )
    loose_MTbestTopJet    = array( 'd', [ 0. ] )
    loose_MTbJet          = array( 'd', [ 0. ] )
    loose_MT2             = array( 'd', [ 0. ] )
    
    nominal_bestTopJetMass  = array( 'd', [ 0. ] )
    nominal_MTbestTopJet    = array( 'd', [ 0. ] )
    nominal_MTbJet          = array( 'd', [ 0. ] )
    nominal_MT2             = array( 'd', [ 0. ] )

    eventWt   = array( 'd', [ 0. ] )
    puWt      = array( 'd', [ 0. ] )

    passElVeto     = array( 'b', [ 0 ] )
    passMuVeto     = array( 'b', [ 0 ] )
    passTauVeto    = array( 'b', [ 0 ] )
    passIsoTrkVeto = array( 'b', [ 0 ] )
    passLeptonVeto = array( 'b', [ 0 ] )
    passTopTaggerL = array( 'b', [ 0 ] )
    passTopTaggerN = array( 'b', [ 0 ] )


    tree.Branch( 'nJets',    nJets,    'nJets/I' )
    tree.Branch( 'nJetsCSVM',nJetsCSVM,'nJetsCSVM/I' )
    tree.Branch( 'nJetsCSVT',nJetsCSVT,'nJetsCSVT/I' )

    tree.Branch( 'metVal',     metVal,     'metVal/D' )
    tree.Branch( 'dphi1',      dphi1,      'dphi1/D' )
    tree.Branch( 'dphi2',      dphi2,      'dphi2/D' )
    tree.Branch( 'dphi3',      dphi3,      'dphi3/D' )
    tree.Branch( 'dphi4',      dphi4,      'dphi4/D' )
    tree.Branch( 'dphiMinCSVM',dphiMinCSVM,'dphiMinCSVM/D' )
    tree.Branch( 'dphiMinCSVT',dphiMinCSVT,'dphiMinCSVT/D' )

    tree.Branch( 'ra2_nJets',    ra2_nJets,    'ra2_nJets/I' )
    tree.Branch( 'ra2_htVal',    ra2_htVal,    'ra2_htVal/D' )
    tree.Branch( 'ra2_mhtVal',   ra2_mhtVal,   'ra2_mhtVal/D' )
    tree.Branch( 'ra2_dphi1',    ra2_dphi1,    'ra2_dphi1/D' )
    tree.Branch( 'ra2_dphi2',    ra2_dphi2,    'ra2_dphi2/D' )
    tree.Branch( 'ra2_dphi3',    ra2_dphi3,    'ra2_dphi3/D' )

    tree.Branch( 'bosonPt',    bosonPt,   'bosonPt/D' )

    tree.Branch( 'photonPt',   photonPt,   'photonPt/D' )
    tree.Branch( 'photonEta',  photonEta,  'photonEta/D' )
    tree.Branch( 'photonMinDR',photonMinDR,'photonMinDR/D' )
    
    tree.Branch( 'muon1Pt',    muon1Pt,    'muon1Pt/D' )
    tree.Branch( 'muon1Eta',   muon1Eta,   'muon1Eta/D' )
    tree.Branch( 'muon1MinDR', muon1MinDR, 'muon1MinDR/D' )
    tree.Branch( 'muon2Pt',    muon2Pt,    'muon2Pt/D' )
    tree.Branch( 'muon2Eta',   muon2Eta,   'muon2Eta/D' )
    tree.Branch( 'muon2MinDR', muon2MinDR, 'muon2MinDR/D' )
    tree.Branch( 'dimuonPt',   dimuonPt,   'dimuonPt/D' )
    tree.Branch( 'dimuonEta',  dimuonEta,  'dimuonEta/D' )
    tree.Branch( 'dimuonMinDR',dimuonMinDR,'dimuonMinDR/D' )
    tree.Branch( 'dimuonM',    dimuonM,    'dimuonM/D' )

    tree.Branch( 'jet1Pt', jet1Pt, 'jet1Pt/D' )
    tree.Branch( 'jet1Eta',jet1Eta,'jet1Eta/D' )
    tree.Branch( 'jet2Pt', jet2Pt, 'jet2Pt/D' )
    tree.Branch( 'jet2Eta',jet2Eta,'jet2Eta/D' )
    tree.Branch( 'jet3Pt', jet3Pt, 'jet3Pt/D' )
    tree.Branch( 'jet3Eta',jet3Eta,'jet3Eta/D' )
    tree.Branch( 'jet4Pt', jet4Pt, 'jet4Pt/D' )
    tree.Branch( 'jet4Eta',jet4Eta,'jet4Eta/D' )

    tree.Branch( 'loose_bestTopJetMass',loose_bestTopJetMass,'loose_bestTopJetMass/D')
    tree.Branch( 'loose_MTbestTopJet',  loose_MTbestTopJet  ,'loose_MTbestTopJet/D')
    tree.Branch( 'loose_MTbJet',        loose_MTbJet        ,'loose_MTbJet/D')
    tree.Branch( 'loose_MT2',           loose_MT2           ,'loose_MT2/D')
                                                                                        
    tree.Branch( 'nominal_bestTopJetMass',nominal_bestTopJetMass,'nominal_bestTopJetMass/D')
    tree.Branch( 'nominal_MTbestTopJet',  nominal_MTbestTopJet  ,'nominal_MTbestTopJet/D')
    tree.Branch( 'nominal_MTbJet',        nominal_MTbJet        ,'nominal_MTbJet/D')
    tree.Branch( 'nominal_MT2',           nominal_MT2           ,'nominal_MT2/D')

    tree.Branch( 'passElVeto'    ,passElVeto    ,   'passElVeto/O' )
    tree.Branch( 'passMuVeto'    ,passMuVeto    ,   'passMuVeto/O' )
    tree.Branch( 'passTauVeto'   ,passTauVeto   ,   'passTauVeto/O' )
    tree.Branch( 'passIsoTrkVeto',passIsoTrkVeto,   'passIsoTrkVeto/O' )
    tree.Branch( 'passLeptonVeto',passLeptonVeto,   'passLeptonVeto/O' )
    tree.Branch( 'passTopTaggerL',passTopTaggerL,   'passTopTaggerL/O' )
    tree.Branch( 'passTopTaggerN',passTopTaggerN,   'passTopTaggerN/O' )

    tree.Branch( 'eventWt',  eventWt,   'eventWt/D' )
    tree.Branch( 'puWt',     puWt,      'puWt/D' )

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
        nVtx[0]         = evNominal.ra2_Vertices 
        nJets[0]        = evNominal.ra2_nJetsPt30Eta24
        nJetsCSVM[0]    = evNominal.ra2_nJetsCSVM
        nJetsCSVT[0]    = evNominal.ra2_nJetsCSVT
        puWt[0]         = evNominal.ra2_PUWt
        eventWt[0]      = evNominal.ra2_EventWt/sfCorr
        metVal[0]       = evNominal.ra2_MET
        dphi1[0]        = evNominal.ra2_dPhiMET1
        dphi2[0]        = evNominal.ra2_dPhiMET2
        dphi3[0]        = evNominal.ra2_dPhiMET3
        dphi4[0]        = evNominal.ra2_dPhiMET4
        dphiMinCSVM[0]  = evNominal.ra2_dPhiMETMinBCSVM
        dphiMinCSVT[0]  = evNominal.ra2_dPhiMETMinBCSVT
        jet1Pt[0]       = evNominal.ra2_Jet1Pt
        jet1Eta[0]      = evNominal.ra2_Jet1Eta
        jet2Pt[0]       = evNominal.ra2_Jet2Pt
        jet2Eta[0]      = evNominal.ra2_Jet2Eta
        jet3Pt[0]       = evNominal.ra2_Jet3Pt
        jet3Eta[0]      = evNominal.ra2_Jet3Eta
        jet4Pt[0]       = evNominal.ra2_Jet4Pt
        jet4Eta[0]      = evNominal.ra2_Jet4Eta

        csvmJet1Pt[0]      = event.ra2_CSVMJet1Pt
        csvmJet1Eta[0]     = event.ra2_CSVMJet1Eta
        csvmJet2Pt[0]      = event.ra2_CSVMJet2Pt
        csvmJet2Eta[0]     = event.ra2_CSVMJet2Eta
        csvtJet1Pt[0]      = event.ra2_CSVTJet1Pt
        csvtJet1Eta[0]     = event.ra2_CSVTJet1Eta
        csvtJet2Pt[0]      = event.ra2_CSVTJet2Pt
        csvtJet2Eta[0]     = event.ra2_CSVTJet2Eta

        ra2_nJets[0]  = evNominal.ra2_nJetsPt50Eta25
        ra2_htVal[0]  = evNominal.ra2_HT
        ra2_mhtVal[0] = evNominal.ra2_MHT
        ra2_dphi1[0]  = evNominal.ra2_dPhiMHT1
        ra2_dphi2[0]  = evNominal.ra2_dPhiMHT2
        ra2_dphi3[0]  = evNominal.ra2_dPhiMHT3

        bosonPt[0]     = 1000
        photonPt[0]    = -10
        photonEta[0]   =  10
        photonMinDR[0] = -1
        
        muon1Pt[0]     = -10
        muon1Eta[0]    =  10
        muon1MinDR[0]  = -1
        muon2Pt[0]     = -10
        muon2Eta[0]    =  10
        muon2MinDR[0]  = -1
        dimuonPt[0]    = -10
        dimuonEta[0]   =  10
        dimuonMinDR[0] = -1
        dimuonM[0]     = -10

        passElVeto[0]     = evNominal.ra2_passElVeto
        passTauVeto[0]    = evNominal.ra2_passTauVeto
        passIsoTrkVeto[0] = evNominal.ra2_passIsoTrkVeto

        if options.sample=="zmumu200" or options.sample=="zmumu400" or options.sample=="muon2012a" or options.sample=="muon2012b" or options.sample=="muon2012ab" or options.sample=="muon2012c" or options.sample=="muonall":
            muon1Pt[0]     = evNominal.ra2_Muon1Pt
            muon1Eta[0]    = evNominal.ra2_Muon1Eta
            muon1MinDR[0]  = evNominal.ra2_Muon1MinDR
            muon2Pt[0]     = evNominal.ra2_Muon2Pt
            muon2Eta[0]    = evNominal.ra2_Muon2Eta
            muon2MinDR[0]  = evNominal.ra2_Muon2MinDR
            dimuonPt[0]    = evNominal.ra2_DiMuonPt
            dimuonEta[0]   = evNominal.ra2_DiMuonEta
            dimuonMinDR[0] = evNominal.ra2_DiMuonMinDR
            dimuonM[0]     = evNominal.ra2_DiMuonInvM
            bosonPt[0]     = dimuonPt[0]
            passMuVeto[0]  = True

        elif options.sample=="zinv50" or options.sample=="zinv100" or options.sample=="zinv200" or options.sample=="zinv400":
            bosonPt[0]     = 1000
            passMuVeto[0] = evNominal.ra2_passMuVeto

        else:
            photonPt[0]    = evNominal.ra2_Photon1Pt
            photonEta[0]   = evNominal.ra2_Photon1Eta
            photonMinDR[0] = evNominal.ra2_Photon1MinDR
            
            bosonPt[0]    = photonPt[0]
            passMuVeto[0] = evNominal.ra2_passMuVeto

        ##########
        passLeptonVeto[0] = passElVeto[0] and passTauVeto[0] and passIsoTrkVeto[0] and passMuVeto[0]

        ##top tagger variables 
        loose_bestTopJetMass[0] = evLoose.ra2_bestTopJetMass 
        loose_MTbestTopJet[0]   = evLoose.ra2_TbestTopJet    
        loose_MTbJet[0]         = evLoose.ra2_TbJet          
        loose_MT2[0]            = evLoose.ra2_MT2            
        passTopTaggerL[0]       = cutF.topTaggerCuts(evLoose)

        nominal_bestTopJetMass[0] = evNominal.ra2_bestTopJetMass 
        nominal_MTbestTopJet[0]   = evNominal.ra2_TbestTopJet    
        nominal_MTbJet[0]         = evNominal.ra2_TbJet          
        nominal_MT2[0]            = evNominal.ra2_MT2            
        passTopTaggerN[0]         = cutF.topTaggerCuts(evNominal)

        sTopSelection = jet1Pt[0] > 70 and abs(jet1Eta[0]) < 2.4 and jet2Pt[0] > 70 and abs(jet2Eta[0]) < 2.4 and metVal[0] > 75 and bosonPt[0] > 75
        ra2Selection = ra2_mhtVal[0] > 75 and bosonPt[0] > 75 and ra2_htVal[0] > 300 and ra2_dphi1[0] > 0.5 and ra2_dphi2[0] > 0.5 and ra2_dphi3[0] > 0.3
        passSelection = sTopSelection or ra2Selection
        if (passSelection):
            tree.Fill()
        #########
        i = i + 1

    #####        
    outputFile.cd()

    ####
    outputFile.Write()
    outputFile.Close()

##########################
if __name__ == '__main__':
    main()
    

#  LocalWords:  elif
