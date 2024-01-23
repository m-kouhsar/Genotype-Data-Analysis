
if [ "$CombinedInputs" == yes ]
then
  if [ "$InputFormat" == zip ]
  then
    echo "Unzipping input data..."
    for i in {1..22}
    do
	    unzip -P $Password  ${InDir}/chr_${i}.zip -d $InDir 
    done
	
    echo "Merging input data..."
    if [ -f ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt ]
    then
      rm ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt
    fi
    for i in {1..22}
    do 
      plink --vcf  ${InDir}/chr${i}.dose.vcf.gz --double-id --make-bed --out ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}_chr${i}
      echo ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}_chr${i} >> ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt 
    done
    plink --merge-list ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt --make-bed --out ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}
    
    while read p
    do
      rm ${p}.*
    done < ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt
  fi
  
  
  if [ "$InputFormat" == vcf ]
  then
	echo "Merging input data..."
    if [ -f ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt ]
    then
      rm ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt
    fi
  
    for i in {1..22}
    do 
      plink --vcf  ${InDir}/chr${i}.dose.vcf.gz --double-id --make-bed --out ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}_chr${i}
      echo ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}_chr${i} >> ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt 
    done
    plink --merge-list ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt --make-bed --out ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}
    
    while read p
    do
      rm ${p}.*
    done < ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt
    
  fi
  
  if [ "$InputFormat" == ped-map ]
  then
	echo "Merging input data..."
    if [ -f ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt ]
    then
      rm ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt
    fi
  
    for i in ${InDir}/*.map
    do 
	    Name=${i%".map"}
      plink --file ${InDir}/${Name} --make-bed --out  ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}_${Name}
      echo ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}_${Name} >> ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt 
    done
    plink --merge-list ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt --make-bed --out ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}
    
    while read p
    do
      rm ${p}.*
    done < ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt
    
  fi
  
  if [ "$InputFormat" == binary ]
  then
	echo "Merging input data..."
    if [ -f ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt ]
    then
      rm ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt
    fi
  
    for i in ${InDir}/*.bed
    do 
	    Name=${i%".bed"}
      echo ${OutDir}/QCoutput_${FilePrefix}/${Name} >> ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt 
    done
    plink --merge-list ${OutDir}/QCoutput_${FilePrefix}/binary.list.txt --make-bed --out ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}
    
  fi
fi


if [ "$CombinedInputs" == no ]
then
  if [ "$InputFormat" == zip ]
  then
    echo "Unzipping input data..."
    for i in ${InDir}/*.zip
    do
		  Name=${i#"${InDir}/}"}
		  Name=${i%",zip"}
	    unzip -P $Password  ${InDir}/chr_${i}.zip -d ${OutDir}/QCoutput_${FilePrefix}
		  mv ${OutDir}/QCoutput_${FilePrefix}/Name ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}
    done
  fi
  
  if [ "$InputFormat" == vcf ]
  then
    for i in ${InDir}/*.vcf.gz
    do 
      plink --vcf  ${InDir}/${i} --double-id --make-bed --out ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}
    done
  fi
  
  if [ "$InputFormat" == ped-map ]
  then
    for i in ${InDir}/*.map
    do 
	    Name=${i%".map"}
      plink --file ${InDir}/${Name} --make-bed --out  ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}
    done
  fi
  
  if [ "$InputFormat" == binary ]
  then
    for i in ${InDir}/*.bed
    do 
      Name=${i%".bed"}
	    cp ${Name}.fam ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}.fam
      cp ${Name}.bim ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}.bim
      cp ${Name}.bed ${OutDir}/QCoutput_${FilePrefix}/${FilePrefix}.bed
    done
  fi
fi
