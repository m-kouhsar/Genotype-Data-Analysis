#!/bin/bash


cd ${OutDir}/QCoutput_${FilePrefix}

echo
echo "[INFO] Removing duplicate samples..."
echo
## some samples were run twice (two different brain regions) identify these:
king -b ${FilePrefix}.bed --duplicate --prefix ${FilePrefix}_king

if [ -f ${FilePrefix}_king.con ]
then
  tail -n +2 ${FilePrefix}_king.con > ${FilePrefix}_king.tmp
	cut -f 1,2 ${FilePrefix}_king.tmp > ${FilePrefix}_duplicateSamples.tmp
	cut -f 3,4 ${FilePrefix}_king.tmp >> ${FilePrefix}_duplicateSamples.tmp
	sort ${FilePrefix}_duplicateSamples.tmp | uniq > ${FilePrefix}_duplicateSamples.txt
	rm ${FilePrefix}_duplicateSamples.tmp
  
fi

##if any duplicate samples calc missingness rate to identify those to exclude
if [ -s ${FilePrefix}_duplicateSamples.txt ]
then 
	plink --bfile ${FilePrefix} --missing --out ${FilePrefix}_duplicateSamples

	## use python script to identify duplicated with greatest missingness
	python ${ScriptDir}/ExcludeDuplicates.py ${FilePrefix}_king.tmp ${FilePrefix}_duplicateSamples.imiss ${FilePrefix}_dupsToExclude.txt
	rm ${FILEPREFIX}_king.tmp
	## remove duplicates
	plink --bfile ${FilePrefix} --remove ${FilePrefix}_dupsToExclude.txt --make-bed --out ${FilePrefix}_update_1
else
	cp ${FilePrefix}.bed ${FilePrefix}_update_1.bed
	cp ${FilePrefix}.bim ${FilePrefix}_update_1.bim
	cp ${FilePrefix}.fam ${FilePrefix}_update_1.fam
fi	

echo
echo "[INFO] Removing variants at the same position (i.e. triallelic)..."
echo
## remove variants at the same position (i.e. triallelic)
awk '{if ($1 != 0) print $1":"$4}' ${FilePrefix}_update_1.bim > pos.tmp
sort pos.tmp | uniq -d > dupLocs.txt
awk -F ":" '{print $1,$2-1,$2,"set1", "set2"}' dupLocs.txt > ${FilePrefix}_positionsExclude.txt

plink --bfile ${FilePrefix}_update_1 --exclude range ${FilePrefix}_positionsExclude.txt --make-bed --out ${FilePrefix}_update_2

rm pos.tmp
rm dupLocs.txt

if [ $CheckSex = "yes" ]
then

  echo
  echo "[INFO] Performing sex check on samples with enough data..."
  echo
  ## perform sex check on samples with enough data
  plink --bfile ${FilePrefix}_update_2 --check-sex --out ${FilePrefix}

  ## we keep samples which do not have sex predicted but check their predicted sex with other softwares such as GenomeStudio
  ## exclude mismatched samples
  ## If a samples doesn't have real sex, we can use the predicted sex for it
  awk '{if ( $4 != $3 && $3 != 0 && $4 != 0 ) print $1,$2,$3,$4,$5,$6}'  ${FilePrefix}.sexcheck >> ${FilePrefix}_sexErrors.txt
  plink --bfile ${FilePrefix}_update_2 --remove ${FilePrefix}_sexErrors.txt --make-bed --out ${FilePrefix}_update_3

else
	cp ${FilePrefix}_update_2.bed ${FilePrefix}_update_3.bed
	cp ${FilePrefix}_update_2.bim ${FilePrefix}_update_3.bim
	cp ${FilePrefix}_update_2.fam ${FilePrefix}_update_3.fam
fi

echo
echo "[INFO] Checking the homozygosity..."
echo
## check for runs of homozygosity
awk '{if ($1 >= 1 && $1 <= 22) print $2}' ${FilePrefix}_update_3.bim > autosomalVariants.txt
plink --bfile ${FilePrefix}_update_3 --extract autosomalVariants.txt --maf $MAF --hwe $HWE --mind $MIND --geno $GENO --indep-pairwise 5000 1000 0.2 --out ${FilePrefix}_ld.auto

plink --bfile ${FilePrefix}_update_3 --extract ${FilePrefix}_ld.auto.prune.in --het --out ${FilePrefix}_roh
## exclude anyone with |Fhet| > 0.2
awk '{if ($6 > 0.2 || $6 < -0.2) print $1,$2}' ${FilePrefix}_roh.het > ${FilePrefix}_excessHet.txt
plink --bfile ${FilePrefix}_update_3 --remove ${FilePrefix}_excessHet.txt --make-bed --out ${FilePrefix}_update_4

echo
echo "[INFO] removing variants with 3+ alleles..."
echo
##remove variants with 3+ alleles
plink --bfile ${FilePrefix}_update_4 --biallelic-only strict  --make-bed --out ${FilePrefix}_update_5

echo
echo "[INFO] Final General QC Step: Filtering sample and variant missingness, HWE, rare variants and exclude variants with no position..."
echo
## filter sample and variant missingness, HWE, rare variants and exclude variants with no position
awk '{if ($1 == 0) print $2}' ${FilePrefix}_update_5.bim > ${FilePrefix}_noLocPos.tmp
plink --bfile ${FilePrefix}_update_5 --exclude ${FilePrefix}_noLocPos.tmp --maf $MAF --hwe $HWE --mind $MIND --geno $GENO --make-bed --out ${FilePrefix}_QC_final

## write list of samples that passed QC for CNV calling
#cut -f 1,2 --delimiter=" " QCoutput_${FilePrefix}/${FilePrefix}_QC.fam > QCoutput_${FilePrefix}/${FilePrefix}_ID_Map.txt
cut -f 2 --delimiter=" " ${FilePrefix}_QC_final.fam > ${FilePrefix}_PastQCSamples.txt

echo
echo "[INFO] Cleaning up intermediate files but keep log files..."
echo
## clean up intermediate files but keep log files
rm ${FilePrefix}_update_*.*

## calc PCS within sample only
# LD prune
#plink --bfile ${FilePrefix}_QC_final --indep 50 5 1.5 --out ${FilePrefix}_QC_final.ld
#plink --bfile ${FilePrefix}_QC_final --extract ${FilePrefix}_QC_final.ld.prune.in --make-bed --out ${FilePrefix}_QC_final.ld.prune

#mkdir -p GCTA_${FilePrefix}

#gcta64 --bfile ${FilePrefix}_QC_final.ld.prune --make-grm-bin --autosome --out GCTA_${FilePrefix}/${FilePrefix}_QC_final_GCTA
#gcta64 --grm GCTA_${FilePrefix}/${FilePrefix}_QC_final_GCTA --pca --thread-num 16 --out GCTA_${FilePrefix}/${FilePrefix}_QC_final.pca

#rm ${FilePrefix}_QC_final.ld.prune*

## extract SNP probes for comparison with DNAm data
#plink --bfile ${FilePrefix}_QC_final --extract ${RefProbe} --recodeA --out ${FilePrefix}_QC_final_DNAmSNPs


echo "General QC was done"
