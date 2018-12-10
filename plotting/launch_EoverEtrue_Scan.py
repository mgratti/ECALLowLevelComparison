
import sys
import os
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

for iEl in prodVersions_list:

  prodVersion = iEl['name']
  anaName = iEl['anaName']

  command = "python EoverEtrue_Scan.py -e {e} -p {p} -a {a}".format(e=ecalVersion, p=prodVersion, a=anaName)
  print "\n Going to run EoverEtrue analysis for ecalVersion={e}, prodVersion={p}, anaName={a}".format(e=ecalVersion, p=prodVersion, a=anaName)
  os.system(command)

  # copy the htaccess file
  command = 'cp HTACCESS plots/anaEoverEtrue_{p}_{e}/.htaccess'.format(e=ecalVersion, p=prodVersion)
  os.system(command)
