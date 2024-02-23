#!/bin/bash
#SBATCH -N 1 -c 10
#SBATCH --mem=60G 
#SBATCH -t 1-0 ## this requests one day of walltime 
#SBATCH --tmp=300G
#SBATCH --job-name=0-PREPROCESS
#SBATCH --output=job_output/0-PREPROCESS.%A_%a.out
#SBATCH --error=job_output/0-PREPROCESS.%A_%a.error
#SBATCH --array=1-870

cd temp
echo "Starting task $SLURM_ARRAY_TASK_ID"
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" 0-DOWNLOAD.sh)

echo $DIR

eval $DIR


