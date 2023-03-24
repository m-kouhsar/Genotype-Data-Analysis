
cd $1

echo "Combine Imputated files (chr 1 .. 22)"
#### COMBINE VCF FILES
bcftools concat chr{1..22}.dose.vcf.gz -o ${2}/${FILEPREFIX}_postImpute_1-22.dose.vcf.gz

###Remove Tri-allellic variants
vcftools --gzvcf ${2}/${FILEPREFIX}_postImpute_1-22.dose.vcf.gz --min-alleles 2 --max-alleles 2 --recode --out ${2}/${FILEPREFIX}_postImpute_1-22_BiAllelic

###Create plink files  
plink --vcf  ${2}/${FILEPREFIX}_postImpute_1-22_BiAllelic.recode.vcf --double-id --make-bed --out ${2}/${FILEPREFIX}_PostImp_1-22_BiAllelic

