//
// $Id: PATObjectSelector.h,v 1.7 2010/06/17 10:08:39 gpetrucc Exp $
//

#ifndef SandBox_Skims_RA2ObjectSelector_h
#define SandBox_Skims_RA2ObjectSelector_h

#include "DataFormats/Common/interface/RefVector.h"

#include "CommonTools/UtilAlgos/interface/StringCutObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/ObjectSelector.h"
#include "CommonTools/UtilAlgos/interface/SingleElementCollectionSelector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PFParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"

#include "ZInvisibleBkgds/Photons/plugins/RA2BasicJetSelector.h"

#include <vector>


namespace pat {


  /* typedef SingleObjectSelector< */
  /*             std::vector<Jet>, */
  /*             StringCutObjectSelector<Jet> */
  /*         > RA2BasicJetSelector; */
  typedef SingleObjectSelector<
              std::vector<Jet>,
              StringCutObjectSelector<Jet>,
              edm::RefVector<std::vector<Jet> >
          > RA2BasicJetRefSelector;
}

#endif
