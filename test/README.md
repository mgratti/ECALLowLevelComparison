# Run-time !

## Examples of usage
cmsRun ../python/Cfg_std.py maxEvents=100

cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_relValZee.root inputFiles_load=samples/RelValZEE_RECO_zeroThreshold.txt 

## Tests
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_nuGun_v2.root inputFiles=file:/afs/cern.ch/work/m/mratti/private/Generation/CMSSW_10_0_3/src/test_generation/test/mg_test_1000evts/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=1000

cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_relValZee_v2.root maxEvents=1000 inputFiles_load=samples/reduced_RelValZEE_RECO_zeroThreshold.txt 
