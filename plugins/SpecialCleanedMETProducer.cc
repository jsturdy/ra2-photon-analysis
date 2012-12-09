// -*- C++ -*-
//
// Package:    Objects
// Class:      Objects
// 
/**\class Objects SpecialCleanedMETProducer.cc ZInvisibleBkgds/Photons/plugins/SpecialCleanedMETProducer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: SpecialCleanedMETProducer.cc,v 1.1 2012/08/30 09:45:14 sturdy Exp $
//
//


#include "ZInvisibleBkgds/Photons/interface/SpecialCleanedMETProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//
// constants, enums and typedefs
//
namespace zinvtools {
  typedef SpecialCleanedMETProducer<pat::Electron, 2> SpecialPATElectronCleanedMETProducer;
  typedef SpecialCleanedMETProducer<pat::Muon, 2>     SpecialPATMuonCleanedMETProducer;
  //  typedef SpecialCleanedMETProducer<pat::Tau, 2>      SpecialPATTauCleanedMETProducer;
  typedef SpecialCleanedMETProducer<pat::Photon, 1>   SpecialPATPhotonCleanedMETProducer;
}  

//define this as a plug-in
using namespace zinvtools;
DEFINE_FWK_MODULE(SpecialPATElectronCleanedMETProducer);
DEFINE_FWK_MODULE(SpecialPATMuonCleanedMETProducer);
//DEFINE_FWK_MODULE(SpecialPATTauCleanedMETProducer);
DEFINE_FWK_MODULE(SpecialPATPhotonCleanedMETProducer);
