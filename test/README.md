# Run-time !

## Example
cmsRun ../python/Cfg_std.py maxEvents=100

cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_relValZee.root inputFiles_load=samples/RelValZEE_RECO_zeroThreshold.txt


## Rechits & PFRechits

### neutrino gun
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_nuGun_v2.root inputFiles=file:/afs/cern.ch/work/m/mratti/private/Generation/CMSSW_10_0_3/src/test_generation/test/mg_test_1000evts/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=1000

cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_nuGun_v10.root inputFiles=file:/mnt/t3nfs01/data01/shome/mratti/cmssw_workarea/Generation/CMSSW_10_0_3/src/test_generation/mg_test_10000evts/sgejob-9706668/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=10000

### relval Zee
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_relValZee_v2.root maxEvents=1000 inputFiles_load=samples/reduced_RelValZEE_RECO_zeroThreshold.txt

### neutrino gun in full-readout 
cmsRun ../python/Cfg_std.py outputFile=test_nuGun_fullReadout_v1.root maxEvents=1000 inputFiles_load=/shome/mratti/cmssw_workarea/Generation/CMSSW_10_0_3/src/test_generation/mg_test_EcalSRSettingsRcd/SingleNuE10_GEN_SIM_DIGI_RECO_FullReadout.root 

### neutrino gun with correct PFrechit collection saved 
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_nuGun_MOD.root inputFiles=file:/shome/mratti/cmssw_workarea/Generation/CMSSW_10_0_3_mod/src/test_generation/mg_test_1000evts/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=1000 

## PFClusters

### double photon 
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_v4.root inputFiles=file:/shome/mratti/cmssw_workarea/Generation/temp_doublePhoton/test_1000evts/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=1000

cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_v5.root inputFiles=file:/shome/mratti/cmssw_workarea/Generation/CMSSW_10_0_3_mod/src/test_generation/doublePhoton_10000evts/sgejob-2266125/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=10000

### double photon with changed seeding and gathering thresholds
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_seedingGathering_v0.root  inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v0/NEVTS10000_seed1.0_GATHER1.0/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=10000

cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_seed0.5_GATHER0.5_v2.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v0/NEVTS10000_seed0.5_GATHER0.5/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=10000
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_seed0.5_GATHER1.0_v2.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v0/NEVTS10000_seed0.5_GATHER1.0/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=10000
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_seed0.5_GATHER2.0_v2.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v0/NEVTS10000_seed0.5_GATHER2.0/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=10000
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_seed1.0_GATHER0.5_v2.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v0/NEVTS10000_seed1.0_GATHER0.5/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=10000
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_seed1.0_GATHER1.0_v2.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v0/NEVTS10000_seed1.0_GATHER1.0/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=10000
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_seed1.0_GATHER2.0_v2.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v0/NEVTS10000_seed1.0_GATHER2.0/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=10000
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_seed2.0_GATHER0.5_v2.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v0/NEVTS10000_seed2.0_GATHER0.5/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=10000
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_seed2.0_GATHER1.0_v2.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v0/NEVTS10000_seed2.0_GATHER1.0/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=10000
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_photonGun_seed2.0_GATHER2.0_v2.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v0/NEVTS10000_seed2.0_GATHER2.0/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=10000



## Basic plotting / monitoring of output 
root -l -b -q "loopdir.C(\"outputfiles/test_relValZee_v3_numEvent1000.root\", \"ecalnoisestudy\")"

## Example to copy to web page
scp test_relValZee_v3_numEvent1000.root.pdf mratti@lxplus.cern.ch:/eos/user/m/mratti/www/EcalDPG/zee/diagnostic_plots/.
