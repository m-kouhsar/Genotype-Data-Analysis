
#!/bin/bash

## format files for use with Michegan Imputation Server

## EXECUTION
# sh SNPArray/preprocessing/4_formatForimputation.sh <population> <SNP ref file>
# where 
# <population > is 3 letter code for super population state ALL for no subsetting by population
# <SNP ref file> is an input file of 
# script needs to be executed from <git repo>/array/

## REQUIRES the following variables in config file
# PROCESSDIR, IMPUTEDIR, FILEPREFIX

## REQUIRES the following software
# plink, perl,

## INPUT
#  # binary plink files following prelim QC

## OUTPUT
# vcf files split by chr for upload to michegan imputation server
set -e

cd $OutDir

mkdir -p ImputationInput

in_file="${FilePrefix}_QC_final"

if [[ ! -f "${FilePrefix}_QC_final.bim" ]]; then
  echo "Warning: There is no general QC output file (${in_file}.bim/bed/fam). "
  echo "         Using original input (${FilePrefix}.bim/bed/fam) instead."
  in_file="${FilePrefix}"
fi

plink --bfile ${in_file} --freq --out ImputationInput/${in_file}

cp ${in_file}.bim ImputationInput/${in_file}.bim
cp ${in_file}.bed ImputationInput/${in_file}.bed
cp ${in_file}.fam ImputationInput/${in_file}.fam

cd ImputationInput
perl ${ScriptDir}/HRC-1000G-check-bim.pl -b ${in_file}.bim -f ${in_file}.frq -r $RefGenome_legend -g --1000g

# sed -i  '1 s+'"${FilePrefix}_QC_final"'+'"${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}_QC_final"'+' Run-plink.sh
if [ "$Addchr" == "yes" ]
then
	sed -i '6,$ s/--make-bed/--recode vcf --output-chr chrM/' Run-plink.sh
else
	sed -i '6,$ s/--make-bed/--recode vcf/' Run-plink.sh
fi

sed -i '1i\set e' Run-plink.sh
sed -i '1i\#!/bin/bash' Run-plink.sh
bash Run-plink.sh
#for file in *.vcf; do awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ${file} > with_chr_${file}; vcf-sort with_chr_${file} | bgzip -c > ${file}.gz;done
for file in *.vcf; do vcf-sort ${file} | bgzip -c > ${file}.gz;done
 
rm *.vcf
#rm *.txt *.log
rm *.b* *.f*

