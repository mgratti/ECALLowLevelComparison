#dir='anaRechits_SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_FR_ecalV13'
#dir='anaRechits_SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_FR_ULPFrecHits_ecalV13'
#dir='anaRechits_SingleNu_Run3_105X_upgrade2018_realistic_v3_450ifb_FR_ecalV13'
#dir='anaRechits_SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_ULPFrecHits_ecalV13'
#dir='anaRechits_SingleNu_Run3_105X_upgrade2018_realistic_v3_450ifb_ecalV13'
#dir='anaEoverEtrue_DoublePhoton_450ifb_nominal_ecalV14'
#dir='anaEoverEtrue_DoublePhoton_180ifb_NewSeed4_ecalV14DR5'
#dir='anaRechits_SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_FR_07_06_19_ecalV15'
#dir='anaRechits_SingleNu_Run3_105X_upgrade2018_realistic_v3_450ifb_FR_07_06_19_ecalV15'
#dir='anaRechits_SingleNu_Run3_105X_upgrade2018_realistic_v3_550ifb_FR_07_06_19_ecalV15'
dir='anaRechits_SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_13_08_19_3sSeed_noPU_ecalV15'
#dir='anaRechits_SingleNu_Run2_105X_upgrade2018_realistic_v3_180ifb_13_08_19_3sSeed_wAutumn18PU_ecalV15'
path='./plots'

cd $path
cp ../HTACCESS ${dir}/.htaccess
tar -czvf ${dir}.tar.gz ${dir}

scp ${dir}.tar.gz mratti@lxplus.cern.ch:/eos/user/m/mratti/www/EcalDPG/RecHits/.
#scp ${dir}.tar.gz mratti@lxplus.cern.ch:/eos/user/m/mratti/www/EcalDPG/PFclusters/.

cd ..

