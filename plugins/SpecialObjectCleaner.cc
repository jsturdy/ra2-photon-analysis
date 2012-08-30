// -*- C++ -*-
//
// Package:    Objects
// Class:      Objects
// 
/**\class Objects SpecialObjectCleaner.cc ZInvisibleBkgds/Photons/plugins/SpecialObjectCleaner.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: SpecialObjectCleaner.cc,v 1.1 2012/08/19 23:45:25 sturdy Exp $
//
//


#include "ZInvisibleBkgds/Photons/interface/SpecialObjectCleaner.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//
// constants, enums and typedefs
//
namespace zinvtools {
  typedef SpecialObjectCleaner<pat::Electron> SpecialPATElectronCleanedJetCollection;
  typedef SpecialObjectCleaner<pat::Muon>     SpecialPATMuonCleanedJetCollection;
  typedef SpecialObjectCleaner<pat::Tau>      SpecialPATTauCleanedJetCollection;
  typedef SpecialObjectCleaner<pat::Photon>   SpecialPATPhotonCleanedJetCollection;
}  

//define this as a plug-in
using namespace zinvtools;
DEFINE_FWK_MODULE(SpecialPATElectronCleanedJetCollection);
DEFINE_FWK_MODULE(SpecialPATMuonCleanedJetCollection);
DEFINE_FWK_MODULE(SpecialPATPhotonCleanedJetCollection);
DEFINE_FWK_MODULE(SpecialPATTauCleanedJetCollection);
