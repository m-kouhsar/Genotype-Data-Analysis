#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq
#SBATCH --time=25:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --job-name=SNPImp

### print start date and time
echo Job started on:
date -u
###

module purge 2>/dev/null
module load VCFtools 2>/dev/null
module load BCFtools 2>/dev/null
module load Python 2>/dev/null
module load R 2>/dev/null

source $1

echo "#############################################################################"
echo " "
echo "Input parameters:"
echo "    Files prefix: $FilePrefix"
echo "    Imputed data directory: $InDir"
echo "    Output data directory: $OutDir"
echo "    Ziped files password: $Password"
echo "    Minor allel frequency: $MAF"
echo "    SNP missingness threshold: $GENO"
echo "    Samples missingness threshold: $MIND"
echo "    Hardy-Weinberg equilibrium exact test p-value threshold: $HWE"
echo "    Is ethnicity of the samples checked? $CheckEthnicity"
echo "    Is reletness of the samples checked? $CheckReletedness"
echo "    Is Sex of the samples checked? $CheckSex"
echo "    Input data format: $InputFormat"
echo "    Chromosomal Separated inputs: $CombinedInputs"
echo "    Genrating summary plots for imputed data? $sumplots"
echo "    Reference genome binary files prefix: $RefGenome_binary"
echo "    Reference genome legend file: $RefGenome_legend"
echo "    Reference genome samples file: $RefGenome_samples"
echo "    Reference genome ped file: $RefGenome_ped"
echo "    Reference genome info file: $RefGenome_info"
echo "    Scripts directory: $ScriptDir"
echo " "
echo "#############################################################################"

mkdir -p ${OutDir}/QCoutput_${FilePrefix}

echo "###############################################################"
echo "Preparing Input data..."
echo "###############################################################"
sh ${ScriptDir}/1_PreparingInputs.sh

if [ "$sumplots" == yes ]
then
  mkdir -p ${OutDir}/QCoutput_${FilePrefix}/Summarize
  echo "###############################################################"
  echo "Generating summerize imputation plots..."
  echo "###############################################################"
  Rscript ${ScriptDir}/2_summarizeImputation.r ${InDir}  ${OutDir}/QCoutput_${FilePrefix}/Summarize  ${RefGenome_legend}
fi

echo "###############################################################"
echo "Running general QC..."
echo "###############################################################"
sh ${ScriptDir}/3_QC.sh

if [ $CheckReletedness = "yes" ]
then
  echo "#######################################################"
  echo "Checking relatedness relatedness..."
  echo "#######################################################"
  sh ${ScriptDir}/4_CheckRelatedness.sh
fi

if [ $CheckEthnicity = "yes" ]
then
  echo "#######################################################"
  echo "Checkinging ethnicity..."
  echo "#######################################################"
  sh ${ScriptDir}/5_CheckEthnicity.sh
fi

if [ $FormatForImputation = "yes" ]
	echo "#######################################################"
	echo "Preparing Imputation input files..."
	echo "#######################################################"
	sh ${ScriptDir}/6_FormatForImputation.sh
fi
#### print end date and time
echo Job finished:
date -u


