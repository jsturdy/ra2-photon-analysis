
#fix something that always gave me problems, probably due to my own misconfiguration
perl -pi -e 's/process.//g' $CMSSW_BASE/src/SandBox/Skims/python/RA2CleaningFilterResults_cfg.py
rm $CMSSW_BASE/src/ZInvisibleBkgds/Photons/plugins/AddElectronUserData.cc 
rm $CMSSW_BASE/src/ZInvisibleBkgds/Photons/plugins/AddMuonUserData.cc 
rm $CMSSW_BASE/src/ZInvisibleBkgds/Photons/plugins/AddPhotonUserData.cc 
rm $CMSSW_BASE/src/ZInvisibleBkgds/Photons/plugins/RA2BasicJetSelector.h 
rm $CMSSW_BASE/src/ZInvisibleBkgds/Photons/plugins/RA2ObjectSelector.cc 
rm $CMSSW_BASE/src/ZInvisibleBkgds/Photons/plugins/RA2ObjectSelector.h 
rm $CMSSW_BASE/src/ZInvisibleBkgds/Photons/interface/AddElectronUserData.h 
rm $CMSSW_BASE/src/ZInvisibleBkgds/Photons/interface/AddMuonUserData.h 
rm $CMSSW_BASE/src/ZInvisibleBkgds/Photons/interface/AddPhotonUserData.h
