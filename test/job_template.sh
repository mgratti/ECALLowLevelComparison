#!/bin/bash

echo $#;
if [ $# -le 5 ]; then
    echo "USAGE: ${0} analysisName iseeding igathering ecalVersion prodVersion nevents ";
    exit 1;
fi

ANANAME=$1
iSEED=$2
iGATHER=$3
EL=$4
PL=$5
NEVENTS=$6

############ BATCH QUEUE DIRECTIVES ##############################
# Lines beginning with #$ are used to set options for the SGE
# queueing system (same as specifying these options to the qsub
# command

# Job name (defines name seen in monitoring by qstat and the
#     job script's stderr/stdout names)
#$ -N example_job

### Specify the queue on which to run
#$ -q all.q

# Change to the current working directory from which the job got
# submitted (will also result in the job report stdout/stderr being
# written to this directory)
#$ -cwd



# here you could change location of the job report stdout/stderr files
#  if you did not want them in the submission directory
#  #$ -o /shome/username/mydir/
#  #$ -e /shome/username/mydir/

##################################################################



source $VO_CMS_SW_DIR/cmsset_default.sh
shopt -s expand_aliases


cmsenv 
echo "going to run"

cmsRun ../python/Cfg_std.py outputFile="outputfiles/test_${ANANAME}_seed"${iSEED}"_gather"${iGATHER}"_"${PL}"_"${EL}".root" inputFiles="/store/user/mratti/EcalGen/${PL}/NEVTS50000_seed"${iSEED}"_GATHER"${iGATHER}"/EGM-RunIISpring18_GEN_SIM_DIGI_RECO.root" maxEvents=${NEVENTS} anaName=${ANANAME}

echo "finished running"

