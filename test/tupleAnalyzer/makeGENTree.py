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
    parser.add_option('-s', action="store",      default="gjets400v1",dest="sample", type="string")
    parser.add_option('-f', action="store",      default="0",    dest="subsec", type="int")
    parser.add_option('-t', action="store",      default="analysisGEN",dest="treeName", type="string")
    parser.add_option('-x', action="store",      default=0,dest="minht",    type="float")
    parser.add_option('-y', action="store",      default=0,dest="minmht",   type="float")
    parser.add_option('-p', action="store",      default=100,dest="minpt",    type="float")
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

    gjetsSamples = ["gjets200","gjets400","gjets400v1","gjets400v2"]
    zinvSamples = ["zinv200","zinv400","zinv400v1","zinv400v2"]
    zmumuSamples = ["zmumu200","zmumu400","zmumu400v1","zmumu400v2"]
    outFileName = "genTree_mht%d_ht%d_%s_%s_%d.root"%(options.minmht,options.minht,options.treeName,options.sample,options.subsec)
    if options.debug:
        outFileName = "genTree_mht%d_ht%d_%s_%s_test.root"%(options.minmht,options.minht,options.treeName,options.sample)

    print outFileName
    sys.stdout.flush()
    outFile  = r.TFile(outFileName,"RECREATE")
    tree     = r.TTree( 'gen', 'tree for gen info ' )
    
    nVtx      = array( 'i', [ 0 ] )
    nJetsMHT  = array( 'i', [ 0 ] )
    nJetsHT   = array( 'i', [ 0 ] )
    nJetsHTMInv = array( 'i', [ 0 ] )
    htVal       = array( 'd', [ 0. ] )
    htMInvVal   = array( 'd', [ 0. ] )
    mhtVal    = array( 'd', [ 0. ] )

    dphi1     = array( 'd', [ 0. ] )
    dphi2     = array( 'd', [ 0. ] )
    dphi3     = array( 'd', [ 0. ] )

    nGenJetsMHT  = array( 'i', [ 0 ] )
    nGenJetsHT   = array( 'i', [ 0 ] )
    htGenVal  = array( 'd', [ 0. ] )
    mhtGenVal = array( 'd', [ 0. ] )
    gendphi1     = array( 'd', [ 0. ] )
    gendphi2     = array( 'd', [ 0. ] )
    gendphi3     = array( 'd', [ 0. ] )

    genBosonM     = array( 'd', [ 0. ] )
    genBosonPt    = array( 'd', [ 0. ] )
    genBosonEta   = array( 'd', [ 0. ] )
    genBosonMinDR = array( 'd', [ 0. ] )

    photonPt  = array( 'd', [ 0. ] )
    photonEta = array( 'd', [ 0. ] )
    photonMinDR = array( 'd', [ 0. ] )
    photonPixV = array( 'b', [ 0 ] )
    photonCSEV = array( 'b', [ 0 ] )
    photonIso  = array( 'b', [ 0 ] )

    ##daughter1M     = array( 'd', [ 0. ] )
    ##daughter1Pt    = array( 'd', [ 0. ] )
    ##daughter1Eta   = array( 'd', [ 0. ] )
    ##
    ##daughter2M     = array( 'd', [ 0. ] )
    ##daughter2Pt    = array( 'd', [ 0. ] )
    ##daughter2Eta   = array( 'd', [ 0. ] )

    #jet1Pt    = array( 'd', [ 0. ] )
    #jet1Eta   = array( 'd', [ 0. ] )
    #jet2Pt    = array( 'd', [ 0. ] )
    #jet2Eta   = array( 'd', [ 0. ] )
    #jet3Pt    = array( 'd', [ 0. ] )
    #jet3Eta   = array( 'd', [ 0. ] )
    #jet4Pt    = array( 'd', [ 0. ] )
    #jet4Eta   = array( 'd', [ 0. ] )
    eventWt   = array( 'd', [ 0. ] )
    puWt      = array( 'd', [ 0. ] )

    passAcceptance  = array( 'b', [ 0 ] )
    passRecoID      = array( 'b', [ 0 ] )
    passRecoTightID = array( 'b', [ 0 ] )
    passRecoIDPixV  = array( 'b', [ 0 ] )
    passRecoIDCSEV  = array( 'b', [ 0 ] )
    passRecoIDPFIso = array( 'b', [ 0 ] )

    passRA2ElVeto    = array( 'b', [ 0 ] )
    passRA2MuVeto    = array( 'b', [ 0 ] )
    passDirIsoElVeto = array( 'b', [ 0 ] )
    passDirIsoMuVeto = array( 'b', [ 0 ] )
    passIsoTrkVeto   = array( 'b', [ 0 ] )


    tree.Branch( 'nVtx',       nVtx,        'nVtx/I' )
    tree.Branch( 'nJetsMHT',   nJetsMHT,    'nJetsMHT/I' )
    tree.Branch( 'nJetsHT',    nJetsHT,     'nJetsHT/I' )
    tree.Branch( 'nJetsHTMInv',nJetsHTMInv, 'nJetsHTMInv/I' )
    tree.Branch( 'mhtVal',     mhtVal,      'mhtVal/D' )
    tree.Branch( 'htVal',      htVal,       'htVal/D' )
    tree.Branch( 'htMInvVal',  htMInvVal,   'htMInvVal/D' )
    tree.Branch( 'dphi1',      dphi1,       'dphi1/D' )
    tree.Branch( 'dphi2',      dphi2,       'dphi2/D' )
    tree.Branch( 'dphi3',      dphi3,       'dphi3/D' )

    tree.Branch( 'nGenJetsMHT',   nGenJetsMHT,    'nGenJetsMHT/I' )
    tree.Branch( 'nGenJetsHT',    nGenJetsHT,     'nGenJetsHT/I' )
    tree.Branch( 'mhtGenVal',     mhtGenVal,      'mhtGenVal/D' )
    tree.Branch( 'htGenVal',      htGenVal,       'htGenVal/D' )
    tree.Branch( 'gendphi1',      gendphi1,       'gendphi1/D' )
    tree.Branch( 'gendphi2',      gendphi2,       'gendphi2/D' )
    tree.Branch( 'gendphi3',      gendphi3,       'gendphi3/D' )

    tree.Branch( 'genBosonM',    genBosonM,    'genBosonM/D' )
    tree.Branch( 'genBosonPt',   genBosonPt,   'genBosonPt/D' )
    tree.Branch( 'genBosonEta',  genBosonEta,  'genBosonEta/D' )
    tree.Branch( 'genBosonMinDR',genBosonMinDR,'genBosonMinDR/D' )

    tree.Branch( 'photonPt',   photonPt,   'photonPt/D' )
    tree.Branch( 'photonEta',  photonEta,  'photonEta/D' )
    tree.Branch( 'photonMinDR',photonMinDR,'photonMinDR/D' )
    tree.Branch( 'photonPixV', photonPixV, 'photonPixV/O' )
    tree.Branch( 'photonCSEV', photonCSEV, 'photonCSEV/O' )
    tree.Branch( 'photonIso',  photonIso,  'photonIso/O' )
    
    #tree.Branch( 'daughter1M',    daughter1M,    'daughter1M/D' )
    #tree.Branch( 'daughter1Pt',   daughter1Pt,   'daughter1Pt/D' )
    #tree.Branch( 'daughter1Eta',  daughter1Eta,  'daughter1Eta/D' )
    #
    #tree.Branch( 'daughter2M',    daughter2M,    'daughter2M/D' )
    #tree.Branch( 'daughter2Pt',   daughter2Pt,   'daughter2Pt/D' )
    #tree.Branch( 'daughter2Eta',  daughter2Eta,  'daughter2Eta/D' )

    #tree.Branch( 'jet1Pt',   jet1Pt,    'jet1Pt/D' )
    #tree.Branch( 'jet1Eta',  jet1Eta,   'jet1Eta/D' )
    #tree.Branch( 'jet2Pt',   jet2Pt,    'jet2Pt/D' )
    #tree.Branch( 'jet2Eta',  jet2Eta,   'jet2Eta/D' )
    #tree.Branch( 'jet3Pt',   jet3Pt,    'jet3Pt/D' )
    #tree.Branch( 'jet3Eta',  jet3Eta,   'jet3Eta/D' )
    #tree.Branch( 'jet4Pt',   jet4Pt,    'jet4Pt/D' )
    #tree.Branch( 'jet4Eta',  jet4Eta,   'jet4Eta/D' )
    tree.Branch( 'eventWt',  eventWt,   'eventWt/D' )
    tree.Branch( 'puWt',     puWt,      'puWt/D' )

    tree.Branch( 'passAcceptance', passAcceptance, 'passAcceptance/O' )
    tree.Branch( 'passRecoID',     passRecoID,     'passRecoID/O' )
    tree.Branch( 'passRecoTightID',passRecoTightID,'passRecoTightID/O' )
    tree.Branch( 'passRecoIDPixV', passRecoIDPixV, 'passRecoIDPixV/O' )
    tree.Branch( 'passRecoIDCSEV', passRecoIDCSEV, 'passRecoIDCSEV/O' )
    tree.Branch( 'passRecoIDPFIso',passRecoIDPFIso,'passRecoIDPFIso/O' )

    tree.Branch( 'passRA2ElVeto'   ,passRA2ElVeto   ,'passRA2ElVeto/O' )
    tree.Branch( 'passRA2MuVeto'   ,passRA2MuVeto   ,'passRA2MuVeto/O' )
    tree.Branch( 'passDirIsoElVeto',passDirIsoElVeto,'passDirIsoElVeto/O' )
    tree.Branch( 'passDirIsoMuVeto',passDirIsoMuVeto,'passDirIsoMuVeto/O' )
    tree.Branch( 'passIsoTrkVeto'  ,passIsoTrkVeto  ,'passIsoTrkVeto/O' )

    myChain = r.TChain('%s/RA2Values'%(options.treeName))
    ##################

    if options.debug:
        if options.sample=="gjets200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets200_trees_mar7/*_?_?_???.root")
        elif options.sample=="gjets400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_v2_trees_mar7/*_?_?_???.root")
        elif options.sample=="gjets400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_trees_mar7/*_?_?_???.root")
        elif options.sample=="gjets400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_v2_trees_mar7/*_?_?_???.root")
        elif options.sample=="zinv200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets200_trees_mar7/*_?_?_???.root")
        elif options.sample=="zinv400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_ext_trees_mar7/*_?_?_???.root")
        elif options.sample=="zinv400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_trees_mar7/*_?_?_???.root")
        elif options.sample=="zinv400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_ext_trees_mar7/*_?_?_???.root")
        elif options.sample=="zmumu200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu200_trees_mar7/*_job?.root")
        elif options.sample=="zmumu400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu400_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu400_ext_trees_mar7/*_?_?_???.root")
        elif options.sample=="zmumu400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu400_trees_mar7/*_?_?_???.root")
        elif options.sample=="zmumu400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu400_ext_trees_mar7/*_?_?_???.root")

    else:
        if options.sample=="gjets200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets200_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="gjets400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_v2_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="gjets400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="gjets400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_v2_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets200_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_ext_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zinvjets400_ext_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zmumu200":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu200_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zmumu400":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu400_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu400_ext_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zmumu400v1":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu400_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zmumu400v2":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/zmumu400_ext_trees_mar7/%s.root"%(subfiles[options.subsec]))

    fChain   = myChain
    idChain  = myChain

    print '%s/RA2Values'%(options.treeName)
    print fChain
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
    for event in fChain:
        # ==============print number of events done == == == == == == == =
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
                # sys.stdout.write("t=7.2f"%(time))
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

        nVtx[0]        = event.ra2_Vertices
        nJetsMHT[0]    = event.ra2_nJetsPt30Eta50
        nJetsHT[0]     = event.ra2_nJetsPt50Eta25
        htVal[0]       = event.ra2_HT
        htMInvVal[0]   = event.ra2_HTMInv
        mhtVal[0]      = event.ra2_MHT
        dphi1[0]       = event.ra2_dPhiMHT1
        dphi2[0]       = event.ra2_dPhiMHT2
        dphi3[0]       = event.ra2_dPhiMHT3

        if options.sample in gjetsSamples:
            htVal[0]          = event.ra2_gen_HT
            mhtVal[0]         = event.ra2_gen_MHT
            dphi1[0]          = event.ra2_gen_dPhiMHT1
            dphi2[0]          = event.ra2_gen_dPhiMHT2
            dphi3[0]          = event.ra2_gen_dPhiMHT3
            nJetsMHT[0]       = event.ra2_gen_nJetsPt30Eta50
            nJetsHT[0]        = event.ra2_gen_nJetsPt50Eta25
            nGenJetsMHT[0]    = event.ra2_gen_nGenJetsPt30Eta50
            nGenJetsHT[0]     = event.ra2_gen_nGenJetsPt50Eta25
            htGenVal[0]       = event.ra2_gen_GenHT
            mhtGenVal[0]      = event.ra2_gen_GenMHT

        if options.sample in zinvSamples:
            htVal[0]          = event.ra2_HT
            mhtVal[0]         = event.ra2_MHT
            dphi1[0]          = event.ra2_dPhiMHT1
            dphi2[0]          = event.ra2_dPhiMHT2
            dphi3[0]          = event.ra2_dPhiMHT3
            nJetsMHT[0]       = event.ra2_nJetsPt30Eta50
            nJetsHT[0]        = event.ra2_nJetsPt50Eta25
            nGenJetsMHT[0]    = event.ra2_nGenJetsPt30Eta50
            nGenJetsHT[0]     = event.ra2_nGenJetsPt50Eta25
            htGenVal[0]       = event.ra2_genHT
            mhtGenVal[0]      = event.ra2_genMHT

        puWt[0]        = event.ra2_PUWt
        eventWt[0]     = event.ra2_EventWt

        genBosonM[0]     = event.ra2_genBoson1M
        genBosonPt[0]    = event.ra2_genBoson1Pt
        genBosonEta[0]   = event.ra2_genBoson1Eta
        genBosonMinDR[0] = event.ra2_genBoson1MinDR

        passRA2ElVeto[0]     = event.ra2_passRA2ElVeto
        passRA2MuVeto[0]     = event.ra2_passRA2MuVeto
        passDirIsoElVeto[0]  = event.ra2_passDirIsoElVeto
        passDirIsoMuVeto[0]  = event.ra2_passDirIsoMuVeto
        passIsoTrkVeto[0]    = event.ra2_passIsoTrkVeto

        passAcceptance[0]  = True
        passRecoID[0]      = True
        passRecoIDPixV[0]  = True
        passRecoIDCSEV[0]  = True
        passRecoIDPFIso[0] = True

        photonPt[0]    = -10.
        photonEta[0]   = -10.
        photonMinDR[0] = 10.
        photonPixV[0] = False
        photonCSEV[0] = False
        photonIso[0]  = False

        #daughter1M[0]     = event.ra2_daughter1M
        #daughter1Pt[0]    = event.ra2_daughter1Pt
        #daughter1Eta[0]   = event.ra2_daughter1Eta
        #
        #daughter2M[0]     = event.ra2_daughter2M
        #daughter2Pt[0]    = event.ra2_daughter2Pt
        #daughter2Eta[0]   = event.ra2_daughter2Eta

        if options.sample in gjetsSamples or options.sample == "qcd":
            passAcceptance[0]  = abs(genBosonEta[0])<1.4442 or (abs(genBosonEta[0])>1.566
                                                                and abs(genBosonEta[0])<2.5)
            passRecoID[0]      = event.ra2_gen_genMatchRecoID
            passRecoIDPixV[0]  = event.ra2_gen_genMatchRecoIDPixV
            passRecoIDCSEV[0]  = event.ra2_gen_genMatchRecoIDCSEV
            passRecoIDPFIso[0] = event.ra2_gen_genMatchRecoIDIso
            photonPt[0]    = event.ra2_Photon1Pt
            photonEta[0]   = event.ra2_Photon1Eta
            photonMinDR[0] = event.ra2_Photon1MinDR
            photonPixV[0]  = event.ra2_Photon1PixelVeto
            photonCSEV[0]  = event.ra2_Photon1EConvVeto
            photonIso[0]   = event.ra2_Photon1IsTightPFIso

        elif options.sample == "zmumu200" or \
                 options.sample == "zmumu400" or \
                 options.sample == "zmumu400v1" or \
                 options.sample == "zmumu400v2":
            passAcceptance[0]  = abs(genBosonEta[0])<2.1

        extra = True
        #if options.sample == "gjets":
        #    extra = cutF.ra2PhotonSelection(event,options.minpt, options.cutDR)
        #if extra and genBosonPt > 100 and mhtVal > 100 and htVal > 300:
        if extra and genBosonPt > 100 and mhtVal > options.minmht and htVal > options.minht:
            tree.Fill()

        #########
        i = i + 1
    #tree.Write()
    outFile.Write()
    outFile.Close()

    ##########################
if __name__ == '__main__':
    main()
    print "done with __main__!"
    sys.stdout.flush()
    


#  LocalWords:  genPassRecoID
