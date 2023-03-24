#!/bin/sh
#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=1:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --job-name=unzip

Pass=${1}   ##${1} is password of the zip files
InDir=${2}
OutDir=${3}

set -e 

## print start date and time
echo Job started on:
date -u

mkdir -p ${OutDir}
for i in ${InDir}/*.zip
do
	unzip -P ${Pass}  $i -d ${OutDir} 
done

## print end date and time
echo Job finished:
date -u


