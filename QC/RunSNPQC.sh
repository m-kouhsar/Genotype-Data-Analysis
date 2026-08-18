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
#SBATCH --job-name=SNPQC
#SBATCH --output=SNPQC.%j.out

##################################################################################
# Running QC on genotype data in plink binary format
# Requiered tools:
#            -plink
#            -king
#            -VCFtools
#            -BCFtools
#            -Python
#            -perl
#            -R
# Required R packages:
#            -data.table
#            -ggplot2
#            -ggpubr
#            -dplyr
##################################################################################

# Exit immediately if a command exits with a non-zero status (-e),
# treat unset variables as an error (-u)
set -eu

### print start date and time
echo Job started on:
date -u
###

if [ -z "$1" ]; then
    echo "Error: No configuration file supplied."
    echo "Usage: sbatch $0 /path/to/your/config.file"
    exit 1
fi

set -a 
source "$1"
set +a

echo "#############################################################################"
echo " "
echo "Input parameters:"
echo "    Input files prefix: $InputPrefix"
echo "    Ziped files password: $Password"
echo "    Minor allel frequency: $MAF"
echo "    SNP missingness threshold: $GENO"
echo "    Samples missingness threshold: $MIND"
echo "    Hardy-Weinberg equilibrium exact test p-value threshold: $HWE"
echo "    Check ethnicity? $CheckEthnicity"
echo "    Check reletness? $CheckReletedness"
echo "    Check sex? $CheckSex"
echo "    Preparing the data for imputation? $FormatForImputation"
echo "    Adding "chr" to chromosomes name in VCF file (see Michigan Imputation Server website)? $Addchr"
echo "    Chromosomal Separated inputs (from imputation server)? $CombinedInputs"
echo "    Genrating summary plots for imputed data? $sumplots"
echo "    Reference genome binary files prefix: $RefGenome_binary"
echo "    Reference genome legend file: $RefGenome_legend"
echo "    Reference genome samples file: $RefGenome_samples"
echo "    Scripts directory: $ScriptDir"
echo " "
echo "#############################################################################"

export InDir=$(dirname -- "$InputPrefix")
export FilePrefix=$(basename -- "$InputPrefix")
export OutDir="${InDir}/${FilePrefix}_QC"

mkdir -p "${OutDir}"

echo "###############################################################"
echo "Preparing Input data..."
echo "###############################################################"
bash "${ScriptDir}/1_PreparingInputs.sh"

if [ "$sumplots" == "yes" ]
then
  mkdir -p "${OutDir}/Summarize"
  echo "###############################################################"
  echo "Generating summerize imputation plots..."
  echo "###############################################################"
  Rscript "${ScriptDir}/2_summarizeImputation.R" "${InDir}"  "${OutDir}/Summarize"  "${RefGenome_legend}"
fi

if [ "$GeneralQC" = "yes" ]
then
  echo "###############################################################"
  echo "Running general QC..."
  echo "###############################################################"
  bash "${ScriptDir}/3_QC.sh"
fi

if [ "$CheckReletedness" = "yes" ]
then
  echo "#######################################################"
  echo "Checking relatedness..."
  echo "#######################################################"
  bash "${ScriptDir}/4_CheckRelatedness.sh"
fi

if [ "$CheckEthnicity" = "yes" ]
then
  echo "#######################################################"
  echo "Checkinging ethnicity..."
  echo "#######################################################"
  bash "${ScriptDir}/5_CheckEthnicity.sh"
fi

if [ "$FormatForImputation" = "yes" ]
then
	echo "#######################################################"
	echo "Preparing Imputation input files..."
	echo "#######################################################"
	bash "${ScriptDir}/6_FormatForImputation.sh"
fi

rm ${OutDir}/${FilePrefix}.bim
rm ${OutDir}/${FilePrefix}.bed
rm ${OutDir}/${FilePrefix}.fam
#### print end date and time
echo Job finished:
date -u


