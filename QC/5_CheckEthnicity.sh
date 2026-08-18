#!/bin/bash
set -eu

cd $OutDir
mkdir -p Ethnicity 
InputData=${FilePrefix}_QC_final
if [[ ! -f "${InputData}.bim" ]]; then
  echo "Warning: There is no general QC output file (${InputData}.bim/bed/fam). "
  echo "         Using original input (${FilePrefix}.bim/bed/fam) instead."
  InputData="${FilePrefix}"
fi

echo
echo "[INFO] change variant ids to chr:bp in input data"
echo

awk '{if ($1 != 0) print $2,"chr"$1":"$4}' ${InputData}.bim > Ethnicity/${FilePrefix}_UpdateVariantsFormat.txt
plink --bfile $InputData --update-name Ethnicity/${FilePrefix}_UpdateVariantsFormat.txt --make-bed --out Ethnicity/${InputData}_1kgIDs

echo
echo "[INFO} change variant ids to chr:bp in reference data"
echo

RefGenome_file=$(basename $RefGenome_binary)

plink --bfile $RefGenome_binary --allow-extra-chr --chr 1-22 --make-bed --out Ethnicity/${RefGenome_file}_1_22
# 1. Parse the .bim file. Keep the first instance as-is. 
# Rename subsequent instances and write those new names to an exclusion file.
awk -v excl_file="Ethnicity/${RefGenome_file}_1_22_dup_variants_to_remove.txt" '{
    if (seen[$2] == 0) {
        seen[$2] = 1;
        print $0;
    } else {
        seen[$2]++;
        # Rename the duplicate by appending a suffix
        $2 = $2 "_dup" seen[$2];
        print $0;
        # Output this renamed ID to the exclusion list using the passed variable
        print $2 > excl_file;
    }
}' OFS='\t' Ethnicity/${RefGenome_file}_1_22.bim > Ethnicity/${RefGenome_file}_1_22_temp.bim

# Copy the binary and family files
cp Ethnicity/${RefGenome_file}_1_22.bed Ethnicity/${RefGenome_file}_1_22_temp.bed
cp Ethnicity/${RefGenome_file}_1_22.fam Ethnicity/${RefGenome_file}_1_22_temp.fam

# 2. Exclude the duplicate instances
plink --bfile Ethnicity/${RefGenome_file}_1_22_temp \
      --exclude Ethnicity/${RefGenome_file}_1_22_dup_variants_to_remove.txt \
      --make-bed \
      --out Ethnicity/${RefGenome_file}_1_22_noDup

# 3. Cleanup temporary files
rm Ethnicity/${RefGenome_file}_1_22_temp.bim
rm Ethnicity/${RefGenome_file}_1_22_temp.bed
rm Ethnicity/${RefGenome_file}_1_22_temp.fam

awk '{if ($1 != 0) print $2,"chr"$1":"$4}' Ethnicity/${RefGenome_file}_1_22_noDup.bim > Ethnicity/${RefGenome_file}_1_22_noDup_UpdateVariantsFormat.txt
plink --bfile Ethnicity/${RefGenome_file}_1_22_noDup --update-name Ethnicity/${RefGenome_file}_1_22_noDup_UpdateVariantsFormat.txt --biallelic-only --make-bed --out Ethnicity/${RefGenome_file}_1_22_noDup_1kgIDs

echo
echo "[INFO] Finding shared SNPs betweed input and reference data"
echo

awk 'NR==FNR {seen[$2]=1; next} seen[$2] {print $2}' Ethnicity/${InputData}_1kgIDs.bim Ethnicity/${RefGenome_file}_1_22_noDup_1kgIDs.bim > Ethnicity/shared_snps.txt

echo "---> number of shared SNPs:"
wc -l Ethnicity/shared_snps.txt

echo
echo "[INFO] Filtering the data based on shared SNPs"
echo

plink --bfile Ethnicity/${InputData}_1kgIDs --extract Ethnicity/shared_snps.txt --make-bed --out Ethnicity/${InputData}_1kgIDs_sharedSNP
plink --bfile Ethnicity/${RefGenome_file}_1_22_noDup_1kgIDs --extract Ethnicity/shared_snps.txt --make-bed --out Ethnicity/${RefGenome_file}_1_22_noDup_1kgIDs_sharedSNP

echo
echo "[INFO] Merging the data with ref genome"
echo
plink --bfile Ethnicity/${InputData}_1kgIDs_sharedSNP --bmerge Ethnicity/${RefGenome_file}_1_22_noDup_1kgIDs_sharedSNP --maf 0.1 --geno 0.1 --make-bed --out Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G || true

if [ -e Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G-merge.missnp ]
then
  echo
  echo "[INFO] Problematic SNPs detected. Attempting to flip strands..."
  echo
  # Step 1: Try flipping the mismatched SNPs in your input data
  plink --bfile Ethnicity/${InputData}_1kgIDs_sharedSNP \
        --flip Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G-merge.missnp \
        --make-bed \
        --out Ethnicity/${InputData}_1kgIDs_sharedSNP_flipped
  mv Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G-merge.missnp  Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G-merge_1.missnp

  echo
  echo "[INFO] Trying to merge again after flipping..."
  echo
  plink --bfile Ethnicity/${InputData}_1kgIDs_sharedSNP_flipped \
        --bmerge Ethnicity/${RefGenome_file}_1_22_noDup_1kgIDs_sharedSNP \
        --maf 0.1 --geno 0.1 \
        --make-bed \
        --out Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G || true
  
  # Step 2: If it STILL fails, exclude the remaining stubborn SNPs from BOTH datasets
  if [ -e Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G-merge.missnp ]
  then
    echo
    echo "[INFO] Flipping didn't resolve all SNPs. Excluding them from both datasets..."
    echo
    
    # Exclude from input data
    plink --bfile Ethnicity/${InputData}_1kgIDs_sharedSNP_flipped \
          --exclude Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G-merge.missnp \
          --make-bed \
          --out Ethnicity/${InputData}_1kgIDs_sharedSNP_flipped_clean
          
    # Exclude from reference data
    plink --bfile Ethnicity/${RefGenome_file}_1_22_noDup_1kgIDs_sharedSNP \
          --exclude Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G-merge.missnp \
          --make-bed \
          --out Ethnicity/${RefGenome_file}_1_22_noDup_1kgIDs_sharedSNP_clean
          
    echo
    echo "[INFO] Final attempt to merge..."
    echo
    plink --bfile Ethnicity/${InputData}_1kgIDs_sharedSNP_flipped_clean \
          --bmerge Ethnicity/${RefGenome_file}_1_22_noDup_1kgIDs_sharedSNP_clean \
          --maf 0.1 --geno 0.1 \
          --make-bed \
          --out Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G
  fi
fi

echo "[INFO] LD prune"
plink --bfile Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G --maf $MAF --hwe $HWE --indep-pairwise 50 5 0.2 --out Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G.ld
plink --bfile Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G --extract Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G.ld.prune.in --make-bed --out Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G.ld.prune

echo "[INFO] use GCTA to calc PCs"
gcta64 --bfile Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G.ld.prune --make-grm-bin --autosome --out Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G.ld.prune
gcta64 --grm Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G.ld.prune --pca --out Ethnicity/${InputData}_1kgIDs_sharedSNP_mergedw1000G.ld.prune.pca

echo
echo "[INFO] Generating plots..."
echo 
Rscript ${ScriptDir}/PlotEthnicity.R ./Ethnicity ${InputData}_1kgIDs_sharedSNP_mergedw1000G ${RefGenome_samples} 

rm Ethnicity/*.bim
rm Ethnicity/*.bed
rm Ethnicity/*.fam

echo "[INFO] Check ethnithity is done"
