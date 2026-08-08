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

plink --bfile $RefGenome_binary --allow-extra-chr --chr 1-22 --make-bed --out Ethnicity/${RefGenome_file}_1_22
# 1. Extract duplicate variant IDs (skipping their first occurrence)
awk '{print $2}' ${RefGenome_binary}.bim | awk 'seen[$0]++' > Ethnicity/${RefGenome_file}_dup_variants_to_remove.txt

# 2. Exclude the duplicate instances
plink --bfile Ethnicity/${RefGenome_file}_1_22 \
      --exclude Ethnicity/${RefGenome_file}_dup_variants_to_remove.txt \
      --make-bed \
      --out Ethnicity/${RefGenome_file}_noDup

awk '{if ($1 != 0) print $2,"chr"$1":"$4}' Ethnicity/${RefGenome_file}_noDup.bim > Ethnicity/${RefGenome_file}_UpdateVariantsFormat.txt
plink --bfile Ethnicity/${RefGenome_file}_noDup --update-name Ethnicity/${RefGenome_file}_UpdateVariantsFormat.txt --biallelic-only --make-bed --out Ethnicity/${RefGenome_file}_Updated
RefGenome_binary="Ethnicity/${RefGenome_file}_Updated"
# first merge with 1000 genomes and filter variants to those in common
# need to test initially in case of error with triallelic variants

echo
echo "[INFO] Merging the data with ref genome"
echo
plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs --bmerge $RefGenome_binary --maf 0.1 --geno 0.1 --make-bed --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G || true

if [ -e Ethnicity/${FilePrefix}_QC_final_mergedw1000G-merge.missnp ]
then
  echo
  echo "[INFO] Problematic SNPs detected. Attempting to flip strands..."
  echo
  # Step 1: Try flipping the mismatched SNPs in your input data
  plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs \
        --flip Ethnicity/${FilePrefix}_QC_final_mergedw1000G-merge.missnp \
        --make-bed \
        --out Ethnicity/${FilePrefix}_QC_final_1kgIDs_flipped

  echo
  echo "[INFO] Trying to merge again after flipping..."
  echo
  plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs_flipped \
        --bmerge $RefGenome_binary \
        --maf 0.1 --geno 0.1 \
        --make-bed \
        --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G || true

  # Step 2: If it STILL fails, exclude the remaining stubborn SNPs from BOTH datasets
  if [ -e Ethnicity/${FilePrefix}_QC_final_mergedw1000G-merge.missnp ]
  then
    echo
    echo "[INFO] Flipping didn't resolve all SNPs. Excluding them from both datasets..."
    echo
    
    # Exclude from input data
    plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs_flipped \
          --exclude Ethnicity/${FilePrefix}_QC_final_mergedw1000G-merge.missnp \
          --make-bed \
          --out Ethnicity/${FilePrefix}_QC_final_1kgIDs.1
          
    # Exclude from reference data
    plink --bfile $RefGenome_binary \
          --exclude Ethnicity/${FilePrefix}_QC_final_mergedw1000G-merge.missnp \
          --make-bed \
          --out Ethnicity/${RefGenome_file}_clean
          
    echo
    echo "[INFO] Final attempt to merge..."
    echo
    plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs.1 \
          --bmerge Ethnicity/${RefGenome_file}_clean \
          --maf 0.1 --geno 0.1 \
          --make-bed \
          --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G
  fi
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
