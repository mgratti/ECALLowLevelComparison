# Run-time !

## Example
```
cmsRun ../python/Cfg_std.py maxEvents=100

cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_relValZee.root inputFiles_load=samples/RelValZEE_RECO_zeroThreshold.txt
```

## Rechits & PFRechits

### neutrino gun
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_nuGun_v2.root inputFiles=file:/afs/cern.ch/work/m/mratti/private/Generation/CMSSW_10_0_3/src/test_generation/test/mg_test_1000evts/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=1000

cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_nuGun_v10.root inputFiles=file:/mnt/t3nfs01/data01/shome/mratti/cmssw_workarea/Generation/CMSSW_10_0_3/src/test_generation/mg_test_10000evts/sgejob-9706668/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=10000
```
### relval Zee
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_relValZee_v2.root maxEvents=1000 inputFiles_load=samples/reduced_RelValZEE_RECO_zeroThreshold.txt
```
### neutrino gun in full-readout
```
cmsRun ../python/Cfg_std.py outputFile=test_nuGun_fullReadout_v1.root maxEvents=1000 inputFiles_load=/shome/mratti/cmssw_workarea/Generation/CMSSW_10_0_3/src/test_generation/mg_test_EcalSRSettingsRcd/SingleNuE10_GEN_SIM_DIGI_RECO_FullReadout.root
```
### neutrino gun with correct PFrechit collection saved
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/test_nuGun_MOD.root inputFiles=file:/shome/mratti/cmssw_workarea/Generation/CMSSW_10_0_3_mod/src/test_generation/mg_test_1000evts/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=1000
```
## PFClusters

### double photon
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/check_std.root inputFiles=file:/shome/mratti/cmssw_workarea/Generation/CMSSW_10_0_3_mod/src/test_generation/doublePhoton_10000evts/sgejob-2295886/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=1000
```
### double photon with changed seeding and gathering thresholds
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/check_v9.root  inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v0/NEVTS10000_seed1.0_GATHER1.0/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=1000
```
### Notes on versions
EcalGen:
```
$PNFS/EcalGen/PROD_SeedingGathering_v1/  --> Run2 conditions, 10K evts, w/o tracker, 100X_upgrade2018_realistic_v7,
$PNFS/EcalGen/PROD_SeedingGathering_v3/  --> Run2 conditions, 50K evts, w/o tracker, 100X_upgrade2018_realistic_v7,
$PNFS/EcalGen/PROD_SeedingGathering_v4/  --> Run3 conditions, 50K evts, w/o tracker, 102X_upgrade2018_realistic_EcalAging_mid2021_235fb_v1
$PNFS/EcalGen/PROD_SeedingGathering_v5/  --> Run3 conditions, 50K evts, w/o tracker, 102X_upgrade2018_realistic_EcalAging_mid2022_315fb_v1
$PNFS/EcalGen/PROD_SeedingGathering_v6/  --> Run3 conditions, 50K evts, w/o tracker, 102X_upgrade2018_realistic_EcalAging_mid2023_400fb_v1
$PNFS/EcalGen/PROD_SeedingGathering_v7/  --> Run3 conditions, 50K evts, w   tracker, 102X_upgrade2018_realistic_EcalAging_mid2021_235fb_v1
```
ECALNoiseStudy:
```
v8 -> with DeltaR=0.05
v9 -> with DeltaR=0.08
```
## Basic plotting / monitoring of output
```
root -l -b -q "loopdir.C(\"outputfiles/test_relValZee_v3_numEvent1000.root\", \"ecalnoisestudy\")"
```
## Example to copy to web page
```
scp test_relValZee_v3_numEvent1000.root.pdf mratti@lxplus.cern.ch:/eos/user/m/mratti/www/EcalDPG/zee/diagnostic_plots/.
```
