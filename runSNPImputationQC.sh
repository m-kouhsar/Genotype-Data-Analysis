#!/bin/sh
#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --job-name=SNPImp

### print start date and time
echo Job started on:
date -u

set -e
####### 

### NOTE: Do not store confidenial information in this file use the config file

######

source $1

mkdir -p ${IMPUTEDIR}/PostImpQCoutput_${FILEPREFIX}

mkdir -p ${IMPUTEDIR}/PostImpQCoutput_${FILEPREFIX}/Summarize

module load R

sh ${SCRIPTDIR}/1_unzip.sh ${Password} ${IMPUTEDIR} ${IMPUTEDIR}
Rscript ${SCRIPTDIR}/2_summarizeImputation.r ${IMPUTEDIR}  ${IMPUTEDIR}/PostImpQCoutput_${FILEPREFIX}/Summarize  ${RefDir}/1000GP_Phase3_combined.legend ALL
#Rscript ${SCRIPTDIR}/6_summarizeImputation.r ${IMPUTEDIR}  ${IMPUTEDIR}/PostImpQCoutput_${FILEPREFIX}/Summarize  ${RefDir}/1000GP_Phase3_combined.legend EUR
####Rscript ${SCRIPTDIR}/preprocessing/6_summarizeImputation.r ${DATADIR}/Imputation_Results/EUR ${RefDir}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab EUR

module purge
module load VCFtools
module load BCFtools
module load Python
#### combine imputation output separately for ALL and EUR versions
sh ${SCRIPTDIR}/3_combineImputationOutput.sh  ${IMPUTEDIR}  PostImpQCoutput_${FILEPREFIX}

sh ${SCRIPTDIR}/4_QC.sh  ${IMPUTEDIR}/PostImpQCoutput_${FILEPREFIX}

sh ${SCRIPTDIR}/5_CheckEthnicity.sh   ${IMPUTEDIR}/PostImpQCoutput_${FILEPREFIX} 
  
sh ${SCRIPTDIR}/6_CheckRelatedness.sh   ${IMPUTEDIR}/PostImpQCoutput_${FILEPREFIX}

#### reformat for use with verifyBamID
####sh preprocessing/reformatForVerifyBamID.sh ${IMPUTEDIR}/ImputationOutput/All/hg38

#### print end date and time
echo Job finished:
date -u


