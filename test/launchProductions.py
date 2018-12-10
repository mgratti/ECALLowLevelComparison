# Script to submit several productions of ecal noise study

import sys
print(sys.version)
import os
import time
import math
import itertools
import subprocess


ecalVersion='ecalV10'

prodVersions_list = [

{ 'name':'PROD_SeedingGathering_v4', 'anaName': 'DoublePhoton'},
{ 'name':'PROD_SeedingGathering_v5', 'anaName': 'DoublePhoton'},
{ 'name':'PROD_SeedingGathering_v6', 'anaName': 'DoublePhoton'},
{ 'name':'PROD_SeedingGathering_v12','anaName': 'DoublePhoton'},

{ 'name':'PROD_SeedingGathering_v7', 'anaName': 'DoubleElectron'},
{ 'name':'PROD_SeedingGathering_v8', 'anaName': 'DoubleElectron'},
{ 'name':'PROD_SeedingGathering_v9', 'anaName': 'DoubleElectron'},
{ 'name':'PROD_SeedingGathering_v11','anaName': 'DoubleElectron'},

]

params={}
params["nevts"] =     [50000]
params["gathering"] = [10.0, 5.0, 1.0, 2.0, 0.5] # multiplier
params["seeding"] =   [10.0, 5.0, 1.0, 2.0, 0.5] # multiplier


##### Submission
for iEl in prodVersions_list:
  prodVersion = iEl['name']
  anaName = iEl['anaName']

  logsDir = 'logs_{p}_{e}'.format(p=prodVersion,e=ecalVersion)

  # create logs dir
  command = 'mkdir -p {l}'.format(l=logsDir)
  if not os.path.isdir(logsDir):
    subprocess.check_output(command, shell=True)
  else: raise RuntimeError('logsDir already present please check')
  
  parameters_set = list(itertools.product(params["nevts"],params["gathering"],params["seeding"] ))
  for iset in parameters_set:
    inevts = iset[0]
    igathering = iset[1]
    iseeding = iset[2]



    command = "qsub -o {l} -e {l} job_template.sh {a} {s} {g} {e} {p} {n}".format(l=logsDir,a=anaName, s=iseeding, g=igathering, e=ecalVersion, p=prodVersion, n=inevts)
    print command
    print "Going to submit job_template.sh for {a} {s} {g} {e} {p} {n}".format(l=logsDir,a=anaName, s=iseeding, g=igathering, e=ecalVersion, p=prodVersion, n=inevts)
    os.system(command)
