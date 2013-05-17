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
    parser.add_option('-s', action="store",      default="signal",dest="template", type="string") ##signal, fit, bkgd
    #parser.add_option('-t', action="store",      default="directPhotons",dest="treeName", type="string")
    parser.add_option('-f', action="store",      default="0",    dest="subsec", type="int")
    options, args = parser.parse_args()

    r.gROOT.SetBatch(True)
    validTemplates = {'signal':"analysisFitTemplate",'fit':"analysisFitTemplate",'bkgd':"analysisFakes"}
    templateName = "%s"%(options.template.lower())
    if templateName not in validTemplates.keys():
        print "You have asked to use template %s"%(templateName)
        print "Please use a valid template from:"
        print validTemplates.keys()
        sys.exit()

    treeName = validTemplates[templateName]
    
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

    outFileName = "purityTree_%s_%d.root"%(templateName,options.subsec)
    if options.debug:
        outFileName = "purityTree_test.root"

    print outFileName
    sys.stdout.flush()
    outFile  = r.TFile(outFileName,"UPDATE")
    tree     = r.TTree( '%s'%(templateName), 'tree for %s template'%(templateName) )
    
    nVtx      = array( 'i', [ 0 ] )
    nJetsMHT  = array( 'i', [ 0 ] )
    nJetsHT   = array( 'i', [ 0 ] )
    htVal     = array( 'd', [ 0. ] )
    mhtVal    = array( 'd', [ 0. ] )
    dphi1     = array( 'd', [ 0. ] )
    dphi2     = array( 'd', [ 0. ] )
    dphi3     = array( 'd', [ 0. ] )
    photonPt        = array( 'd', [ 0. ] )
    photonEta       = array( 'd', [ 0. ] )
    photonMinDR     = array( 'd', [ 0. ] )
    photonPixV      = array( 'b', [ 0 ] )
    photonCSEV      = array( 'b', [ 0 ] )
    photonSigIeIe   = array( 'd', [ 0.0 ] )
    photonIsoPFCH   = array( 'd', [ 0.0 ] )
    photonIsoPFNU   = array( 'd', [ 0.0 ] )
    photonIsoPFGA   = array( 'd', [ 0.0 ] )
    photonIsoPFComb = array( 'd', [ 0.0 ] )

    eventWt   = array( 'd', [ 0. ] )
    puWt      = array( 'd', [ 0. ] )


    tree.Branch( 'nVtx',       nVtx,        'nVtx/I' )
    tree.Branch( 'nJetsMHT',   nJetsMHT,    'nJetsMHT/I' )
    tree.Branch( 'nJetsHT',    nJetsHT,     'nJetsHT/I' )
    tree.Branch( 'mhtVal',     mhtVal,      'mhtVal/D' )
    tree.Branch( 'htVal',      htVal,       'htVal/D' )
    tree.Branch( 'dphi1',      dphi1,       'dphi1/D' )
    tree.Branch( 'dphi2',      dphi2,       'dphi2/D' )
    tree.Branch( 'dphi3',      dphi3,       'dphi3/D' )

    tree.Branch( 'photonPt',       photonPt,       'photonPt/D' )
    tree.Branch( 'photonEta',      photonEta,      'photonEta/D' )
    tree.Branch( 'photonMinDR',    photonMinDR,    'photonMinDR/D' )
    tree.Branch( 'photonPixV',     photonPixV,     'photonPixV/O' )
    tree.Branch( 'photonCSEV',     photonCSEV,     'photonCSEV/O' )
    tree.Branch( 'photonSigIeIe',  photonSigIeIe,  'photonSigIeIe/D' )
    tree.Branch( 'photonIsoPFCH',  photonIsoPFCH,  'photonIsoPFCH/D' )
    tree.Branch( 'photonIsoPFNU',  photonIsoPFNU,  'photonIsoPFNU/D' )
    tree.Branch( 'photonIsoPFGA',  photonIsoPFGA,  'photonIsoPFGA/D' )
    tree.Branch( 'photonIsoPFComb',photonIsoPFComb,'photonIsoPFComb/D' )
    
    tree.Branch( 'eventWt',  eventWt,   'eventWt/D' )
    tree.Branch( 'puWt',     puWt,      'puWt/D' )

    myChain = r.TChain('%s/RA2Values'%(treeName))
    ##################

    if options.debug:
        if templateName=="signal":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets200_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_v2_trees_mar7/*_?_?_???.root")
        elif templateName=="bkgd" or templateName=="fit":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012a_v1_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012a_recoverv1_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012b_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012c_v1_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012c_v2_trees_mar7/*_?_?_???.root")
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012d_v1_trees_mar7/*_?_?_???.root")

    else:
        if templateName=="signal":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets200_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_gjets400_v2_trees_mar7/%s.root"%(subfiles[options.subsec]))
        elif templateName=="bkgd" or templateName=="fit":
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012a_v1_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012a_recoverv1_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012b_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012c_v1_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012c_v2_trees_mar7/%s.root"%(subfiles[options.subsec]))
            myChain.Add("dcache:///pnfs/cms/WAX/11/store/user/sturdy07/RA2_535_Flats/photon_run2012d_v1_trees_mar7/%s.root"%(subfiles[options.subsec]))

    fChain   = myChain

    print '%s/RA2Values'%(treeName)
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
    trueSignalPhotons = 0
    falseSignalPhotons = 0
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
        mhtVal[0]      = event.ra2_MHT
        dphi1[0]       = event.ra2_dPhiMHT1
        dphi2[0]       = event.ra2_dPhiMHT2
        dphi3[0]       = event.ra2_dPhiMHT3
        puWt[0]        = event.ra2_PUWt
        eventWt[0]     = event.ra2_EventWt

        photonPt[0]      = event.ra2_Photon1Pt
        photonEta[0]     = event.ra2_Photon1Eta
        photonMinDR[0]   = event.ra2_Photon1MinDR
        photonSigIeIe[0] = event.ra2_Photon1SigmaIetaIeta
        photonPixV[0]    = event.ra2_Photon1PixelVeto
        photonCSEV[0]    = event.ra2_Photon1EConvVeto
        photonIsoPFCH[0]   = event.ra2_Photon1pfCH
        photonIsoPFNU[0]   = event.ra2_Photon1pfNU
        photonIsoPFGA[0]   = event.ra2_Photon1pfGA
        photonIsoPFComb[0] = event.ra2_Photon1pfCH+event.ra2_Photon1pfNU+event.ra2_Photon1pfGA

        fillTree = False
        if event.ra2_Photon1SigmaIetaIeta > 1:
            print "showerShape variable = %2.4f(%2.4f)"%(event.ra2_Photon1SigmaIetaIeta,photonSigIeIe[0])
        if templateName == "signal":
            if options.debug:
                print event.ra2_Photon1PDGID,event.ra2_Photon1SigmaIetaIeta
            if event.ra2_Photon1PDGID == 22:
                fillTree = True
                trueSignalPhotons = trueSignalPhotons + 1
                
            else:
                #print "signal photon matched to pdgid:%d"%(event.ra2_Photon1PDGID)
                falseSignalPhotons = falseSignalPhotons + 1
        elif templateName == "fit":
            fillTree = True
        elif templateName == "bkgd":
            #print event.ra2_Photon1IsTightPFIso
            if not event.ra2_Photon1IsTightPFIso:
                fillTree = True

        if fillTree:
            tree.Fill()

        #########
        i = i + 1
    print "signalPhotons:%d, fakePhotons:%d"%(trueSignalPhotons,falseSignalPhotons)
    #tree.Write()
    outFile.Write()
    outFile.Close()

    ##########################
if __name__ == '__main__':
    main()
    print "done with __main__!"
    sys.stdout.flush()
    


#  LocalWords:  genPassRecoID
