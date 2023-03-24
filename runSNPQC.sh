#!/bin/sh
#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --job-name=SNPQC


## print start date and time
echo Job started on:
date -u

########

## NOTE: Do not store confidential information in this file use the config file

######


source $1

module purge
module load Python

echo "#######################################################"
echo "#### Step1 main QC"
echo "#######################################################"

sh ${SCRIPTDIR}/1_QC.sh

module load R

echo "#######################################################"
echo "#### Step2 Checking ethnicity"
echo "#######################################################"

sh ${SCRIPTDIR}/2_CheckEthnicity.sh

echo "#######################################################"
echo "#### Step3 Checking relatedness"
echo "#######################################################"
  
sh ${SCRIPTDIR}/3_CheckRelatedness.sh


echo "#######################################################"
echo "#### Step4 Preparing Imputation input files"
echo "#######################################################"
module purge
module load VCFtools
sh ${SCRIPTDIR}/4_formatForImputation.sh ALL ${RefDir}/1000GP_Phase3_combined.legend

## print finish date and time
echo Job finished on:
date -u

