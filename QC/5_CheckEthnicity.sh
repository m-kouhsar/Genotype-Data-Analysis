
cd ${OutDir}/QCoutput_${FilePrefix}
mkdir -p Ethnicity 
# change variant ids to chr:bp
awk '{if ($1 != 0) print $2,"chr"$1":"$4}' ${FilePrefix}_QC_final.bim > Ethnicity/${FilePrefix}_QC_final_updateTo1KGFormat.txt
plink --bfile ${FilePrefix}_QC_final --update-name Ethnicity/${FilePrefix}_QC_final_updateTo1KGFormat.txt --make-bed --out Ethnicity/${FilePrefix}_QC_final_1kgIDs

# first merge with 1000 genomes and filter variants to those in common
# need to test initially in case of error with triallelic variants

plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs --bmerge $RefGenome_binary --maf 0.1 --geno 0.1 --make-bed --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G_test

if [ -e Ethnicity/${FilePrefix}_QC_final_mergedw1000G_test-merge.missnp ]
then
  plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs --exclude Ethnicity/${FilePrefix}_QC_final_mergedw1000G_test-merge.missnp --make-bed --out Ethnicity/${FilePrefix}_QC_final_1kgIDs.1
  plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs.1 --bmerge $RefGenome_binary --maf 0.1 --geno 0.1 --make-bed --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G_test
fi
## issue with variants at same position but different alleles - exclude these
plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs --exclude Ethnicity/${FilePrefix}_QC_final_mergedw1000G_test-merge.missnp --make-bed --out Ethnicity/${FilePrefix}_QC_final_1kgIDs_forMerge

plink --bfile Ethnicity/${FilePrefix}_QC_final_1kgIDs_forMerge --bmerge $RefGenome_binary --maf 0.2 --geno 0.05 --make-bed --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G
# LD prune
plink --bfile Ethnicity/${FilePrefix}_QC_final_mergedw1000G --indep 50 5 1.5 --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G.ld
plink --bfile Ethnicity/${FilePrefix}_QC_final_mergedw1000G --extract Ethnicity/${FilePrefix}_QC_final_mergedw1000G.ld.prune.in --make-bed --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G.ld.prune

rm Ethnicity/${FilePrefix}_QC_final_1kgIDs_forMerge*
rm Ethnicity/${FilePrefix}_QC_final_mergedw1000G_test*

# use GCTA to calc PCs
gcta64 --bfile Ethnicity/${FilePrefix}_QC_final_mergedw1000G.ld.prune --make-grm-bin --autosome --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G
gcta64 --grm Ethnicity/${FilePrefix}_QC_final_mergedw1000G --pca --out Ethnicity/${FilePrefix}_QC_final_mergedw1000G.pca

rm Ethnicity/${FilePrefix}_QC_final_mergedw1000G*grm*

Rscript ${ScriptDir}/PlotEthnicity.r ${OutDir}/QCoutput_${FilePrefix}/Ethnicity ${FilePrefix}_QC_final ${RefGenome_samples} $RefGenome_ped $RefGenome_info

#plink --bfile ${FilePrefix}_QC_final --remove ${FilePrefix}_QC_final_EthnicityOutliers.txt --make-bed --out Ethnicity/${FilePrefix}_QC_final
#mv Relatedness/${FilePrefix}_QC_final.b* ./
#mv Relatedness/${FilePrefix}_QC_final.fam ./

echo "Check ethnithity is done"
