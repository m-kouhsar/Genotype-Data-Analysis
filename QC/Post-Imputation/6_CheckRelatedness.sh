
cd ${1}
echo "Checking Relatedness"
## check for relatedness with other samples with KING
king -b ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd.bed --kinship --prefix ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd_king

## check for relatedness with other samples with plink
plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd --genome --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd_ibd

module load R

Rscript ${SCRIPTDIR}/plotrelatedness.r ${1} ${FILEPREFIX}_PostImp_1-22_BiAllelic

plink --bfile ${FILEPREFIX}_PostImp_1-22_BiAllelic_QCd --remove ${FILEPREFIX}_PostImp_1-22_BiAllelic_RelatednessOutliers.txt  --make-bed --out ${FILEPREFIX}_PostImp_1-22_BiAllelic_QC_final




