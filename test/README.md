# Run-time !

## Example
cmsRun ../python/Cfg_std.py maxEvents=100

cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_relValZee.root inputFiles_load=samples/RelValZEE_RECO_zeroThreshold.txt

## Initial tests
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_nuGun_v2.root inputFiles=file:/afs/cern.ch/work/m/mratti/private/Generation/CMSSW_10_0_3/src/test_generation/test/mg_test_1000evts/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=1000

cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_relValZee_v2.root maxEvents=1000 inputFiles_load=samples/reduced_RelValZEE_RECO_zeroThreshold.txt

## 10K evts nugun
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_nuGun_v10.root inputFiles=file:/mnt/t3nfs01/data01/shome/mratti/cmssw_workarea/Generation/CMSSW_10_0_3/src/test_generation/mg_test_10000evts/sgejob-9706668/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=10000

## 1K evts nugun with PFrechit collection that we actually want to see
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_nuGun_MOD.root inputFiles=file:/shome/mratti/cmssw_workarea/Generation/CMSSW_10_0_3_mod/src/test_generation/mg_test_1000evts/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=1000 

## Basic plotting monitoring of output - gives output parallel to input file
root -l -b -q "loopdir.C(\"outputfiles/test_relValZee_v3_numEvent1000.root\", \"ecalnoisestudy\")"

## Example to copy to web page
scp test_relValZee_v3_numEvent1000.root.pdf mratti@lxplus.cern.ch:/eos/user/m/mratti/www/EcalDPG/zee/diagnostic_plots/.
