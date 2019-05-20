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

### neutrino gun with run-2 vs run-3 tags

```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run2_ecalV9.root inputFiles=/store/user/mratti/EcalGen/GEN_SIM_DIGI/SingleNu/Run2Cond/100X_upgrade2018_realistic_v7/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=10000

cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run2_new_ecalV9.root inputFiles=/store/user/mratti/EcalGen/GEN_SIM_DIGI/SingleNu/Run2Cond/102X_upgrade2018_realistic_v15/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=10000

cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run3_1_ecalV9.root inputFiles=/store/user/mratti/EcalGen/GEN_SIM_DIGI/SingleNu/Run3Cond/102X_upgrade2018_realistic_EcalAging_mid2021_235fb_v1/SingleNuE10_GEN_SIM_DIGI_RECO.root  maxEvents=10000

cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run3_2_ecalV9.root inputFiles=/store/user/mratti/EcalGen/GEN_SIM_DIGI/SingleNu/Run3Cond/102X_upgrade2018_realistic_EcalAging_mid2022_315fb_v1/SingleNuE10_GEN_SIM_DIGI_RECO.root  maxEvents=10000

cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run3_3_ecalV9.root inputFiles=/store/user/mratti/EcalGen/GEN_SIM_DIGI/SingleNu/Run3Cond/102X_upgrade2018_realistic_EcalAging_mid2023_400fb_v1/SingleNuE10_GEN_SIM_DIGI_RECO.root  maxEvents=10000
```

### neutrino gun for UL tests

GT="103X_mc2017_realistic_v2_AB_v01"
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run2_UL_AB_ecalV10.root inputFiles=file:/shome/mratti/cmssw_workarea/Generation/CMSSW_10_4_0_pre2/src/ECALGen/private_production/reco/nugun_reco_103X_AB/RelValNuGun_103X_mc2017_realistic_v2_AB_v01_HS-v1.root maxEvents=15000
```

GT="103X_mc2017_realistic_v2_AC_v01" 
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run2_UL_AC_ecalV10.root inputFiles=file:/shome/mratti/cmssw_workarea/Generation/CMSSW_10_4_0_pre2/src/ECALGen/private_production/reco/nugun_reco_103X_AC/RelValNuGun_103X_mc2017_realistic_v2_AC_v01_HS-v1.root maxEvents=15000
```

GT="94X_mc2017_realistic_v11"
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/MinBias_Run2_Fall17_ecalV10.root inputFiles=file:root://cms-xrd-global.cern.ch//store/mc/RunIIFall17FSPremix/MinBias_TuneCP2_13TeV-pythia8/GEN-SIM-RECO/PUMoriond17_94X_mc2017_realistic_v11-v1/80000/7CB15E49-7BF1-E811-8096-008CFAC919F8.root maxEvents=10000
```

Fall17 conditions
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run2_Fall17_central_ecalV10.root inputFiles=file:root://cms-xrd-global.cern.ch//store/relval/CMSSW_11_4_0_pre2/RelValNuGun/GEN-SIM-RECO/103X_mc2017_realistic_v2-v1/20000/E3B48B18-B84F-4346-BEF4-CA8508B2675F.root     maxEvents=15000
```

THIS SAMPLE IS BUGGED
GT="94X_mc2017_realistic_v10"
```
cmsRun ../python/Cfg_std.py  outputFile=outputfiles/SingleNu_Run2_Fall17_ecalV10.root inputFiles=file:/shome/mratti/cmssw_workarea/Generation/CMSSW_9_4_0_patch1/src/ECALGen/private_production/gen_sim_digi_reco/nugun_Fall17/SingleNuE10_GEN_SIM_DIGI_RECO.root maxEvents=15000
```

### neutrino gun for Run-2 vs Run-3 after March19
edit in Cfg GT="105X_upgrade2018_realistic_v3"

end-of-Run2 50K evts :
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_ecalV13.root inputFiles=/store/user/mratti/EcalGen/GEN_SIM_DIGI/SingleNu/Run2Cond/105X_upgrade2018_realistic_v3/SingleNuE10_GEN_SIM_DIGI_RECO_ULPFrecHits.root maxEvents=50000

```
end-of-Run3 50K evts :
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run3_105X_upgrade2018_realistic_v3_450ifb_ecalV13.root inputFiles=/store/user/mratti/EcalGen/GEN_SIM_DIGI/SingleNu/Run3Cond/105X_upgrade2018_realistic_v3/SingleNuE10_GEN_SIM_DIGI_RECO_50K.root maxEvents=50000
```

Samples with full readout:

end-of-Run2 minimal pfrechit thresholds:
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_FR_ecalV13.root inputFiles=/store/user/mratti/EcalGen/GEN_SIM_DIGI/SingleNu/Run2Cond/105X_upgrade2018_realistic_v3/SingleNuE10_GEN_SIM_DIGI_RECO_FullReadout.root maxEvents=5000
```

end-of-Run2 UL 2018 pfrechit thresholds:
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_FR_ULPFrecHits_ecalV13.root inputFiles=/store/user/mratti/EcalGen/GEN_SIM_DIGI/SingleNu/Run2Cond/105X_upgrade2018_realistic_v3/SingleNuE10_GEN_SIM_DIGI_RECO_FullReadout_ULPFrecHits.root maxEvents=5000
```

end-of-Run3:
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/SingleNu_Run3_105X_upgrade2018_realistic_v3_450ifb_FR_ecalV13.root inputFiles=/store/user/mratti/EcalGen/GEN_SIM_DIGI/SingleNu/Run3Cond/105X_upgrade2018_realistic_v3/SingleNuE10_GEN_SIM_DIGI_RECO_FullReadout.root maxEvents=5000
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
v10 -> as v9 but change eta bin end-caps
```

### double electron
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/check_doubleEle.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v7/NEVTS50000_seed1.0_GATHER1.0/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root maxEvents=1000 anaName="DoubleElectron"
```

$PNFS/EcalGen/PROD_SeedingGathering_v7/  --> Run3 conditions, 50K evts, w   tracker, 102X_upgrade2018_realistic_EcalAging_mid2021_235fb_v1


##  PFClusters for Run-2 vs Run-3 after March19
GT=105X_upgrade2018_realistic_v3_180ifb

Double photon 180/fb:
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/DoublePhoton_180ifb_nominal_ecalV14.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v26/NEVTS150000_seed1.0_GATHER1.0/EGM_GEN_SIM_DIGI_RECO.root maxEvents=150000 anaName="DoublePhoton"
```
Double photon 450/fb:
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/DoublePhoton_450ifb_nominal_ecalV14.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v25/NEVTS150000_seed1.0_GATHER1.0/EGM_GEN_SIM_DIGI_RECO.root maxEvents=150000 anaName="DoublePhoton"
```

Double electron 180/fb:
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/DoubleElectron_180ifb_nominal_ecalV14.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v27/NEVTS150000_seed1.0_GATHER1.0/EGM_GEN_SIM_DIGI_RECO.root maxEvents=150000 anaName="DoubleElectron"
```

Double electron 450/fb:
```
cmsRun ../python/Cfg_std.py outputFile=outputfiles/DoubleElectron_450ifb_nominal_ecalV14.root inputFiles=/store/user/mratti/EcalGen/PROD_SeedingGathering_v28/NEVTS150000_seed1.0_GATHER1.0/EGM_GEN_SIM_DIGI_RECO.root maxEvents=150000 anaName="DoubleElectron"
```


## Basic plotting / monitoring of output
```
root -l -b -q "loopdir.C(\"outputfiles/test_relValZee_v3_numEvent1000.root\", \"ecalnoisestudy\")"
```
## Example to copy to web page
```
scp test_relValZee_v3_numEvent1000.root.pdf mratti@lxplus.cern.ch:/eos/user/m/mratti/www/EcalDPG/zee/diagnostic_plots/.
```
