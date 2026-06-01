#!/bin/bash
set -e 

cd ${OutDir}/QCoutput_${FilePrefix}
mkdir -p Relatedness
## check for relatedness with other samples with KING
king -b ${FilePrefix}_QC_final.bed --kinship --prefix Relatedness/${FilePrefix}_QC_final_king

## check for relatedness with other samples with plink
plink --bfile ${FilePrefix}_QC_final --genome --out Relatedness/${FilePrefix}_QC_final_ibd

echo "[INFO] Genrating the relatedness plots..."

Rscript ${ScriptDir}/PlotRelatedness.R ./Relatedness $Kinship ${FilePrefix}_QC_final


plink --bfile ${FilePrefix}_QC_final --remove Relatedness/${FilePrefix}_QC_final.RelatednessOutliers.txt  --make-bed --out Relatedness/${FilePrefix}_QC_final

mv Relatedness/${FilePrefix}_QC_final.b* ./
mv Relatedness/${FilePrefix}_QC_final.fam ./
