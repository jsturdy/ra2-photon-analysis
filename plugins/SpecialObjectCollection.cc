// -*- C++ -*-
//
// Package:    Objects
// Class:      Objects
// 
/**\class Objects SpecialObjectCollection.cc ZInvisibleBkgds/Photons/plugins/SpecialObjectCollection.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Jared Sturdy
//         Created:  Wed Apr 18 16:06:24 CDT 2012
// $Id: SpecialObjectCollection.cc,v 1.1 2012/08/19 23:45:25 sturdy Exp $
//
//


#include "ZInvisibleBkgds/Photons/interface/SpecialObjectCollection.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//
// constants, enums and typedefs
//
namespace zinvtools {
  typedef SpecialObjectCollection<pat::Electron> SpecialPATElectronCollection;
  typedef SpecialObjectCollection<pat::Muon>     SpecialPATMuonCollection;
  typedef SpecialObjectCollection<pat::Tau>      SpecialPATTauCollection;
  typedef SpecialObjectCollection<pat::Photon>   SpecialPATPhotonCollection;
  typedef SpecialObjectCollection<pat::Jet>      SpecialPATJetCollection;
}  

//define this as a plug-in
using namespace zinvtools;
DEFINE_FWK_MODULE(SpecialPATElectronCollection);
DEFINE_FWK_MODULE(SpecialPATMuonCollection);
DEFINE_FWK_MODULE(SpecialPATPhotonCollection);
DEFINE_FWK_MODULE(SpecialPATTauCollection);
DEFINE_FWK_MODULE(SpecialPATJetCollection);
