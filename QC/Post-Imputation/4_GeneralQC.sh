#!/bin/bash


cd ${1}

if [ ! -f ${FILEPREFIX}_PostImp_1-22_BiAllelic.bed ]
then

  plink --vcf  ${FILEPREFIX}_postImpute_1-22_BiAllelic.recode.vcf --double-id --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic

fi

echo "General QC for Imputed data"

## some samples were run twice (two different brain regions) identify these:
king -b ${FILEPREFIX}_PostImp_1-22_BiAllelic.bed --duplicate --prefix ${FILEPREFIX}_king

if [ -f ${FILEPREFIX}_king.con ]
then
  tail -n +2 ${FILEPREFIX}_king.con > ${FILEPREFIX}_king.tmp
	cut -f 1,2 ${FILEPREFIX}_king.tmp > ${FILEPREFIX}_duplicateSamples.tmp
	cut -f 3,4 ${FILEPREFIX}_king.tmp >> ${FILEPREFIX}_duplicateSamples.tmp
	sort ${FILEPREFIX}_duplicateSamples.tmp | uniq > ${FILEPREFIX}_duplicateSamples.txt
	rm ${FILEPREFIX}_duplicateSamples.tmp
fi

##if any duplicate samples calc missingness rate to identify those to exclude
if [ -s ${FILEPREFIX}_duplicateSamples.txt ]
then 
	plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic --missing --out ${FILEPREFIX}_duplicateSamples

	## use python script to identify duplicated with greatest missingness
	python ${SCRIPTDIR}/ExcludeDuplicates.py ${FILEPREFIX}_king.tmp ${FILEPREFIX}_duplicateSamples.imiss ${FILEPREFIX}_dupsToExclude.txt
  rm ${FILEPREFIX}_king.tmp

	## remove duplicates
	plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic --remove ${FILEPREFIX}_dupsToExclude.txt --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_1
else
	cp ${FILEPREFIX}_PostImp_1-22_BiAllelic.bed  ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_1.bed
	cp ${FILEPREFIX}_PostImp_1-22_BiAllelic.bim  ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_1.bim
	cp ${FILEPREFIX}_PostImp_1-22_BiAllelic.fam  ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_1.fam
fi	


## remove variants at the same position (i.e. triallelic)
awk '{if ($1 != 0) print $1":"$4}' ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_1.bim > pos.tmp
sort pos.tmp | uniq -d > dupLocs.txt
awk -F ":" '{print $1,$2-1,$2,"set1", "set2"}' dupLocs.txt > ${FILEPREFIX}_PostImp_1-22_BiAllelic_positionsExclude.txt

plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_1 --exclude range ${FILEPREFIX}_PostImp_1-22_BiAllelic_positionsExclude.txt --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_2

rm pos.tmp
rm dupLocs.txt
 
## check for runs of homozygosity
awk '{if ($1 >= 1 && $1 <= 22) print $2}' ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_2.bim > ${FILEPREFIX}_PostImp_1-22_BiAllelic_autosomalVariants.txt
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_2 --extract ${FILEPREFIX}_PostImp_1-22_BiAllelic_autosomalVariants.txt --maf 0.01 --hwe 0.00001 --mind 0.02 --geno 0.05 --indep-pairwise 5000 1000 0.2 --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_ld.auto

plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_2 --extract ${FILEPREFIX}_PostImp_1-22_BiAllelic_ld.auto.prune.in --het --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_roh
## exclude anyone with |Fhet| > 0.2
awk '{if ($6 > 0.2 || $6 < -0.2) print $1,$2}' ${FILEPREFIX}_PostImp_1-22_BiAllelic_roh.het > ${FILEPREFIX}_PostImp_1-22_BiAllelic_excessHet.txt
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_2 --remove ${FILEPREFIX}_PostImp_1-22_BiAllelic_excessHet.txt --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_3
#rm autosomalVariants.txt

##remove variants with 3+ alleles
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_3 --biallelic-only strict  --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_4


## filter sample and variant missingness, HWE, rare variants and exclude variants with no position
###remove samples with >5% missing values, and SNPs with > 1% missing values, Hardy-Weinburg equilibrium P < 0.00001, and a minor allele frequency of <5%
awk '{if ($1 == 0) print $2}' ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_4.bim > ${FILEPREFIX}_PostImp_1-22_BiAllelic_noLocPos.tmp
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_4 --exclude ${FILEPREFIX}_PostImp_1-22_BiAllelic_noLocPos.tmp --maf 0.05 --hwe 0.00001 --mind 0.05 --geno 0.01 --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd

### 0 variants removed due to missing genotype data (--geno).
### --hwe: 859 variants removed due to Hardy-Weinberg exact test.
### 39572976 variants removed due to minor allele threshold(s)
### 6576073 variants and 217 people pass filters and QC.

## write list of samples that passed QC for CNV calling
#cut -f 1,2 --delimiter=" " QCoutput_${FILEPREFIX}/${FILEPREFIX}_QCd.fam > QCoutput_${FILEPREFIX}/${FILEPREFIX}_ID_Map.txt
cut -f 2 --delimiter=" " ${FILEPREFIX}_PostImp_1-22_BiAllelic.fam > ${FILEPREFIX}_PostImp_1-22_BiAllelic_PastQCSamples.txt

## clean up intermediate files but keep log files
rm ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_*.b*
rm ${FILEPREFIX}_PostImp_1-22_BiAllelic_update_*.fam
 

## calc PCS within sample only
# LD prune
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd --indep 50 5 1.5 --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd.ld
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd --extract ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd.ld.prune.in --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd.ld.prune

gcta64 --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd.ld.prune --make-grm-bin --autosome --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd_GCTA
gcta64 --grm ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd_GCTA --pca --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd.pca

rm ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd.ld.prune*

###### plot PCs to identify outliers
module load R
Rscript ${SCRIPTDIR}/plotPCs.r ${1} ${FILEPREFIX}_PostImp_1-22_BiAllelic 3


## extract SNP probes for comparison with DNAm data
##plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd --extract ${RefDir}/RSprobes.txt --recodeA --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_DNAmSNPs 

