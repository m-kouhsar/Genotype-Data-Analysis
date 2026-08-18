#!/bin/bash

set -eu


if [ "$CombinedInputs" == yes ]
then
  if [ "$InputFormat" == zip ]
  then
    echo "[INFO] Unzipping input data..."
    for i in {1..22}
    do
	    unzip -P $Password  ${InDir}/chr_${i}.zip -d $InDir 
    done
	
    echo "[INFO] Merging input data..."
    if [ -f ${OutDir}/binary.list.txt ]
    then
      rm ${OutDir}/binary.list.txt
    fi
    for i in {1..22}
    do 
      plink --vcf  ${InDir}/chr${i}.dose.vcf.gz --double-id --make-bed --out ${OutDir}/${FilePrefix}_chr${i}
      echo ${OutDir}/${FilePrefix}_chr${i} >> ${OutDir}/binary.list.txt 
    done
    plink --merge-list ${OutDir}/binary.list.txt --make-bed --out ${OutDir}/${FilePrefix}
    
  fi
  
  
  if [ "$InputFormat" == vcf ]
  then
	echo "[INFO] Merging input data..."
    if [ -f ${OutDir}/binary.list.txt ]
    then
      rm ${OutDir}/binary.list.txt
    fi
  
    for i in {1..22}
    do 
      plink --vcf  ${InDir}/chr${i}.dose.vcf.gz --double-id --make-bed --out ${OutDir}/${FilePrefix}_chr${i}
      echo ${OutDir}/${FilePrefix}_chr${i} >> ${OutDir}/binary.list.txt 
    done
    plink --merge-list ${OutDir}/binary.list.txt --make-bed --out ${OutDir}/${FilePrefix}
    
  fi
  
  if [ "$InputFormat" == ped-map ]
  then
	echo "[INFO] Merging input data..."
    if [ -f ${OutDir}/binary.list.txt ]
    then
      rm ${OutDir}/binary.list.txt
    fi
  
    for i in ${InDir}/*.map
    do 
	    Name=${i%".map"}
      plink --file ${InDir}/${Name} --make-bed --out  ${OutDir}/${FilePrefix}_${Name}
      echo ${OutDir}/${FilePrefix}_${Name} >> ${OutDir}/binary.list.txt 
    done
    plink --merge-list ${OutDir}/binary.list.txt --make-bed --out ${OutDir}/${FilePrefix}
    
  fi
  
  if [ "$InputFormat" == binary ]
  then
	echo "[INFO] Merging input data..."
    if [ -f ${OutDir}/binary.list.txt ]
    then
      rm ${OutDir}/binary.list.txt
    fi
  
    for i in ${InDir}/*.bed
    do 
	    Name=${i%".bed"}
      echo ${OutDir}/${Name} >> ${OutDir}/binary.list.txt 
    done
    plink --merge-list ${OutDir}/binary.list.txt --make-bed --out ${OutDir}/${FilePrefix}
    
  fi
fi


if [ "$CombinedInputs" == no ]
then
  
  if [ "$InputFormat" == vcf ]
  then
      plink --vcf  ${InputPrefix} --double-id --make-bed --out ${OutDir}/${FilePrefix}
  fi
  
  if [ "$InputFormat" == ped-map ]
  then
      plink --file ${InputPrefix} --make-bed --out  ${OutDir}/${FilePrefix}

  fi
  
  if [ "$InputFormat" == binary ]
  then
	    cp ${InputPrefix}.fam ${OutDir}/${FilePrefix}.fam
      cp ${InputPrefix}.bim ${OutDir}/${FilePrefix}.bim
      cp ${InputPrefix}.bed ${OutDir}/${FilePrefix}.bed
  fi
fi
