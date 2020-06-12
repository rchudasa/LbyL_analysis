#!/bin/tcsh
setenv HOME /afs/cern.ch/work/r/rchudasa/private/CMSSW_7_5_8_patch4/src/HeavyIonsAnalysis/PhotonAnalysis/test/aug_reco_data_check_for_lumi/hiforest/id_eff/hybrid
setenv WORK $PWD
cd ${HOME}
cmsenv
cp check_data_driven_id_eff.C ${WORK}
cp filelist_24.txt ${WORK}
cd ${WORK}
g++ check_data_driven_id_eff.C `root-config --cflags --libs` -O2 -o check_data_driven_id_eff.exe
./check_data_driven_id_eff.exe "filelist_24.txt" 
cp test.root /afs/cern.ch/work/r/rchudasa/private/CMSSW_7_5_8_patch4/src/HeavyIonsAnalysis/PhotonAnalysis/test/aug_reco_data_check_for_lumi/hiforest/id_eff/hybrid/ntuple_24.root
cp check_data_driven_id_eff.C /afs/cern.ch/work/r/rchudasa/private/CMSSW_7_5_8_patch4/src/HeavyIonsAnalysis/PhotonAnalysis/test/aug_reco_data_check_for_lumi/hiforest/id_eff/hybrid/check_data_driven_id_eff.C
