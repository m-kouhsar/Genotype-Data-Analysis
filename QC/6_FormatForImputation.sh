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

cd ${OutDir}/QCoutput_${FilePrefix}

mkdir -p ImputationInput


plink --bfile ${FilePrefix}_QC_final --freq --out ImputationInput/${FilePrefix}_QC_final

cp ${FilePrefix}_QC_final.bim ImputationInput/${FilePrefix}_QC_final.bim
cd ImputationInput
perl ${ScriptDir}/HRC-1000G-check-bim.pl -b ${FilePrefix}_QC_final.bim -f ${FilePrefix}_QC_final.frq -r $RefGenome_legend -g --1000g


sed -i  '1 s+'"${FilePrefix}_QC_final"'+'"${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}_QC_final"'+' Run-plink.sh
sed -i  '6,$ s+--make-bed+--recode vcf+' Run-plink.sh
sh Run-plink.sh
#for file in *.vcf; do awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' ${file} > with_chr_${file}; vcf-sort with_chr_${file} | bgzip -c > ${file}.gz;done
for file in *.vcf; do vcf-sort ${file} | bgzip -c > ${file}.gz;done
 
rm *.vcf
rm *.txt *.log
rm *.b* *.f*

