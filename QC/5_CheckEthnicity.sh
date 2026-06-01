#!/bin/bash
set -e

cd ${OutDir}/QCoutput_${FilePrefix}
mkdir -p Ethnicity 
echo
echo "[INFO] change variant ids to chr:bp in input data"
echo
awk '{if ($1 != 0) print $2,"chr"$1":"$4}' ${FilePrefix}_QC_final.bim > Ethnicity/${FilePrefix}_UpdateVariantsFormat.txt
plink --bfile ${FilePrefix}_QC_final --update-name Ethnicity/${FilePrefix}_UpdateVariantsFormat.txt --make-bed --out Ethnicity/${FilePrefix}_QC_final_1kgIDs

echo
echo "[INFO} change variant ids to chr:bp in reference data"
echo
RefGenome_file=$(basename $RefGenome_binary)
awk '{if ($1 != 0) print $2,"chr"$1":"$4}' ${RefGenome_binary}.bim > Ethnicity/${RefGenome_file}_UpdateVariantsFormat.txt
plink --bfile $RefGenome_binary --update-name Ethnicity/${RefGenome_file}_UpdateVariantsFormat.txt --biallelic-only --make-bed --out Ethnicity/${RefGenome_file}_Updated
RefGenome_binary="Ethnicity/${RefGenome_file}_Updated"
# first merge with 1000 genomes and filter variants to those in common
# need to test initially in case of error with triallelic variants

echo
echo "[INFO] Merging the data with ref genome"
echo
plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs --bmerge $RefGenome_binary --maf 0.1 --geno 0.1 --make-bed --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G || true

if [ -e Ethnicity/${FilePrefix}_QC_final_mergedw1000G*.missnp ]
then
  echo
  echo "[INFO} Removing problematic SNPs from the data..."
  echo
  plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs --exclude Ethnicity/${FilePrefix}_QC_final_mergedw1000G*.missnp --make-bed --out Ethnicity/${FilePrefix}_QC_final_1kgIDs.1
  echo
  echo "[INFO] Try to merge again"
  echo
  plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs.1 --bmerge $RefGenome_binary --maf 0.1 --geno 0.1 --make-bed --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G
fi

echo "[INFO] LD prune"
plink --bfile Ethnicity/${FilePrefix}_QC_final_mergedw1000G --indep 50 5 1.5 --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G.ld
plink --bfile Ethnicity/${FilePrefix}_QC_final_mergedw1000G --extract Ethnicity/${FilePrefix}_QC_final_mergedw1000G.ld.prune.in --make-bed --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G.ld.prune

#rm Ethnicity/${FilePrefix}_QC_final_1kgIDs*
#rm Ethnicity/${FilePrefix}_QC_final_mergedw1000G*

echo "[INFO] use GCTA to calc PCs"
gcta64 --bfile Ethnicity/${FilePrefix}_QC_final_mergedw1000G.ld.prune --make-grm-bin --autosome --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G
gcta64 --grm Ethnicity/${FilePrefix}_QC_final_mergedw1000G --pca --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G.pca

#rm Ethnicity/${FilePrefix}_QC_final_mergedw1000G*grm*

echo
echo "[INFO] Generating plots..."
echo 
Rscript ${ScriptDir}/PlotEthnicity.R ./Ethnicity ${FilePrefix}_QC_final ${RefGenome_samples} #$RefGenome_ped $RefGenome_info

#plink --bfile ${FilePrefix}_QC_final --remove ${FilePrefix}_QC_final_EthnicityOutliers.txt --make-bed --out Ethnicity/${FilePrefix}_QC_final
#mv Relatedness/${FilePrefix}_QC_final.b* ./
#mv Relatedness/${FilePrefix}_QC_final.fam ./

echo "[INFO] Check ethnithity is done"
