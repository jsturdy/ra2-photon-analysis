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
    parser.add_option('-l', action="store_true", default=False, dest="applyVetos")
    parser.add_option('-s', action="store",      default="gjets",dest="sample", type="string")
    parser.add_option('-f', action="store",      default="0",    dest="subsec", type="int")
    parser.add_option('-t', action="store",      default="directPhotonsID",dest="treeName", type="string")
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

    outFileName = "genTreeDR%2.1f_%s_%s_%d.root"%(options.cutDR,options.treeName,options.sample,options.subsec)
    if options.debug:
        outFileName = "genTreeDR%2.1f_%s_%s_test.root"%(options.cutDR,options.treeName,options.sample)

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
    genBosonM     = array( 'd', [ 0. ] )
    genBosonPt    = array( 'd', [ 0. ] )
    genBosonEta   = array( 'd', [ 0. ] )
    genBosonMinDR = array( 'd', [ 0. ] )
    photonPt  = array( 'd', [ 0. ] )
    photonEta = array( 'd', [ 0. ] )
    photonMinDR = array( 'd', [ 0. ] )
    jet1Pt    = array( 'd', [ 0. ] )
    jet1Eta   = array( 'd', [ 0. ] )
    jet2Pt    = array( 'd', [ 0. ] )
    jet2Eta   = array( 'd', [ 0. ] )
    jet3Pt    = array( 'd', [ 0. ] )
    jet3Eta   = array( 'd', [ 0. ] )
    eventWt   = array( 'd', [ 0. ] )
    puWt      = array( 'd', [ 0. ] )

    passAcceptance  = array( 'b', [ 0 ] )
    passRecoID      = array( 'b', [ 0 ] )
    passRecoIDPFIso = array( 'b', [ 0 ] )


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

    tree.Branch( 'genBosonM',    genBosonM,    'genBosonM/D' )
    tree.Branch( 'genBosonPt',   genBosonPt,   'genBosonPt/D' )
    tree.Branch( 'genBosonEta',  genBosonEta,  'genBosonEta/D' )
    tree.Branch( 'genBosonMinDR',genBosonMinDR,'genBosonMinDR/D' )

    tree.Branch( 'photonPt',   photonPt,   'photonPt/D' )
    tree.Branch( 'photonEta',  photonEta,  'photonEta/D' )
    tree.Branch( 'photonMinDR',photonMinDR,'photonMinDR/D' )
    
    tree.Branch( 'jet1Pt',   jet1Pt,    'jet1Pt/D' )
    tree.Branch( 'jet1Eta',  jet1Eta,   'jet1Eta/D' )
    tree.Branch( 'jet2Pt',   jet2Pt,    'jet2Pt/D' )
    tree.Branch( 'jet2Eta',  jet2Eta,   'jet2Eta/D' )
    tree.Branch( 'jet3Pt',   jet3Pt,    'jet3Pt/D' )
    tree.Branch( 'jet3Eta',  jet3Eta,   'jet3Eta/D' )
    tree.Branch( 'eventWt',  eventWt,   'eventWt/D' )
    tree.Branch( 'puWt',     puWt,      'puWt/D' )

    tree.Branch( 'passAcceptance', passAcceptance, 'passAcceptance/O' )
    tree.Branch( 'passRecoID',     passRecoID,     'passRecoID/O' )
    tree.Branch( 'passRecoIDPFIso',passRecoIDPFIso,'passRecoIDPFIso/O' )

    vetoSuffix = ""
    if not options.applyVetos:
        vetoSuffix = "NoVeto"
    myChain = r.TChain('%s%s/RA2Values'%(options.treeName,vetoSuffix))
    pfisoChain = r.TChain('%sPFIso%s/RA2Values'%(options.treeName,vetoSuffix))
    ##################

    if options.debug:
        if options.sample=="gjets":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ra2Studies/gjetsht400_gen_tree_535_v1/res/*_?_?_???.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ra2Studies/gjetsht400_gen_tree_535_v1/res/*_?_?_???.root")
        elif options.sample=="qcd":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ra2Studies/qcdpt25to3000_gen_tree_535_v1/res/*_?_?_???.root")
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ra2Studies/qcdpt25to3000_gen_tree_535_v1/res/*_?_?_???.root")
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ra2Studies/zinvht400_gen_tree_535_v1/res/*_?_?_???.root")
    else:
        if options.sample=="gjets":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ra2Studies/gjetsht400_gen_tree_535_v1/res/%s.root"%(subfiles[options.subsec]))
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ra2Studies/gjetsht400_gen_tree_535_v1/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="qcd":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ra2Studies/qcdpt25to3000_gen_tree_535_v1/res/%s.root"%(subfiles[options.subsec]))
            pfisoChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ra2Studies/qcdpt25to3000_gen_tree_535_v1/res/%s.root"%(subfiles[options.subsec]))
        elif options.sample=="zinv":
            myChain.Add("/uscmst1b_scratch/lpc1/3DayLifetime/sturdy07/tmp/cmssw535/treeMaker/ra2Studies/zinvht400_gen_tree_535_v1/res/%s.root"%(subfiles[options.subsec]))

    fChain   = myChain
    isoChain = myChain
    if options.sample=="gjets" or options.sample=="qcd":
        idChain  = myChain
        isoChain = pfisoChain

    print '%s%s/RA2Values'%(options.treeName,vetoSuffix)
    print fChain
    print '%sPFIso%s/RA2Values'%(options.treeName,vetoSuffix)
    print isoChain
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
    for event,pfisoevent in itertools.izip(fChain,isoChain):
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
        nJetsMHT[0]    = event.ra2_nJetsGenPt30Eta50
        nJetsHT[0]     = event.ra2_nJetsGenPt50Eta25
        nJetsHTMInv[0] = event.ra2_nJetsGenPt50Eta25MInv
        htVal[0]       = event.ra2_genHT
        htMInvVal[0]   = event.ra2_genHTMInv
        mhtVal[0]      = event.ra2_genMHT
        dphi1[0]       = event.ra2_genDPhiMHT1
        dphi2[0]       = event.ra2_genDPhiMHT2
        dphi3[0]       = event.ra2_genDPhiMHT3
        jet1Pt[0]      = event.ra2_genJet1Pt
        jet1Eta[0]     = event.ra2_genJet1Eta
        jet2Pt[0]      = event.ra2_genJet2Pt
        jet2Eta[0]     = event.ra2_genJet2Eta
        jet3Pt[0]      = event.ra2_genJet3Pt
        jet3Eta[0]     = event.ra2_genJet3Eta
        puWt[0]        = event.ra2_PUWt
        eventWt[0]     = event.ra2_EventWt

        genBosonM[0]     = event.ra2_genBoson1M
        genBosonPt[0]    = event.ra2_genBoson1Pt
        genBosonEta[0]   = event.ra2_genBoson1Eta
        genBosonMinDR[0] = event.ra2_genBoson1MinDR


        passAcceptance[0]  = True
        passRecoID[0]      = True
        passRecoIDPFIso[0] = True
        photonPt[0]    = -10.
        photonEta[0]   = 10.
        photonMinDR[0] = 10.

        if options.sample == "gjets" or options.sample == "qcd":
            passAcceptance[0]  = abs(genBosonEta[0])<1.4442 or (abs(genBosonEta[0])>1.566 and abs(genBosonEta[0])<2.5)
            passRecoID[0]      = event.ra2_genPassRecoIso
            passRecoIDPFIso[0] = pfisoevent.ra2_genPassRecoIso
            photonPt[0]    = event.ra2_boson1Pt
            photonEta[0]   = event.ra2_boson1Eta
            photonMinDR[0] = event.ra2_boson1MinDR

        extra = True
        #if options.sample == "gjets":
        #    extra = cutF.ra2PhotonSelection(event,options.minpt, options.cutDR)
        if extra and genBosonPt > 100 and mhtVal > 100 and (htVal > 300 or htMInvVal >300):
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
    

