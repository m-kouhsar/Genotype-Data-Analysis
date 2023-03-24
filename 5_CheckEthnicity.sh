
cd ${1}

# change variant ids to chr:bp
awk '{if ($1 != 0) print $2,"chr"$1":"$4}' ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd.bim > ${FILEPREFIX}_PostImp_1-22_BiAllelic_updateTo1KGFormat.txt
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd --update-name ${FILEPREFIX}_PostImp_1-22_BiAllelic_updateTo1KGFormat.txt --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd_1kgIDs

# first merge with 1000 genomes and filter variants to those in common
# need to test initially in case of error with triallelic variants
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd_1kgIDs --bmerge ${RefDir}/1000g_GRCh37_HG19_Referance1.bed ${RefDir}/1000g_GRCh37_HG19_Referance1.bim ${RefDir}/1000g_GRCh37_HG19_Referance1.fam --maf 0.1 --geno 0.1 --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G_test

## issue with variants at same position but different alleles - exclude these
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd_1kgIDs --exclude ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G_test-merge.missnp --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_1kgIDs_forMerge

plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_1kgIDs_forMerge --bmerge ${RefDir}/1000g_GRCh37_HG19_Referance1.bed ${RefDir}/1000g_GRCh37_HG19_Referance1.bim ${RefDir}/1000g_GRCh37_HG19_Referance1.fam --maf 0.2 --geno 0.05 --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G
# LD prune
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G --indep 50 5 1.5 --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G.ld
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G --extract ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G.ld.prune.in --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G.ld.prune

rm ${FILEPREFIX}_PostImp_1-22_BiAllelic_1kgIDs_forMerge*
rm ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G_test*

# use GCTA to calc PCs
gcta64 --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G.ld.prune --make-grm-bin --autosome --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G
gcta64 --grm ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G --pca --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G.pca

rm ${FILEPREFIX}_PostImp_1-22_BiAllelic_mergedw1000G*grm*

module load R

Rscript ${SCRIPTDIR}/plotEthnicity.r ${1} ${FILEPREFIX}_PostImp_1-22_BiAllelic ${RefDir} 

#plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd --remove ${FILEPREFIX}_PostImp_1-22_BiAllelic_EthnicityOutliers.txt --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd

echo "Check ethnithity is done"
