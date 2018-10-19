VERSION="prodV4_ecalV9"
PL="PROD_SeedingGathering_v4"

for iSEED in 0.5 1.0 2.0 5.0 10.0
do
  for iGATHER in 1.0 2.0 5.0 10.0
  do
    cmsRun ../python/Cfg_std.py outputFile="outputfiles/test_photonGun_seed"$iSEED"_gather"$iGATHER"_v"$VERSION".root" inputFiles="/store/user/mratti/EcalGen/${PL}/NEVTS50000_seed"$iSEED"_GATHER"$iGATHER"/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root" maxEvents=50000
  done
done
