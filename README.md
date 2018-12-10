
# ECALNoiseStudy
Package to perform analysis on ECAL reconstruciton objects, by M.G. Ratti.

Started from a slim version of M.Bharti's Validation/EcalValidation code, cloned from C.Rovelli  

Currently set-up to do analysis of 

| Analysis Type  | Neutrinos | Photons | Electrons | Zee |
| -------------  | --------- | --- |      --- |      --- |
| Rechits        | yes  | yes | yes | yes |
| PFrechits      | yes  | yes | yes | yes |
| PFClusters     | no | yes | no  | no |
| SuperClusters  | no | no | yes | TODO |

Studies so far do not include pile-up.

### Installation 
CMSSW release depends on the GT of the samples you want to run on

```
cmsrel CMSSW_10_2_4
cd CMSSW_10_2_4/src
cmsenv

git clone git@github.com:mgratti/ECALNoiseStudy.git ECAL/ECALLowLevelComparison
scram b -j 8
```

### Run - go to test and see examples in README.md
