#!/bin/bash
set -e 

cd $OutDir
mkdir -p Relatedness
in_file="${FilePrefix}_QC_final"

if [[ ! -f "${FilePrefix}_QC_final.bim" ]]; then
  echo "Warning: There is no general QC output file (${in_file}.bim/bed/fam). "
  echo "         Using original input (${FilePrefix}.bim/bed/fam) instead."
  in_file="${FilePrefix}"
fi

## check for relatedness with other samples with KING
king -b ${in_file}.bed --kinship --prefix Relatedness/${in_file}_king

## check for relatedness with other samples with plink
plink --bfile $in_file --genome --out Relatedness/${in_file}_ibd

echo "[INFO] Genrating the relatedness plots..."

Rscript ${ScriptDir}/PlotRelatedness.R ./Relatedness $Kinship ${in_file}

if [ -f "Relatedness/${FilePrefix}_QC_final.RelatednessOutliers.txt" ]; then
  plink --bfile $in_file --remove Relatedness/${in_file}.RelatednessOutliers.txt  --make-bed --out Relatedness/${in_file}
fi
