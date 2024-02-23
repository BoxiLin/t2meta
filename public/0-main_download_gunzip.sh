#!/bin/bash
#SBATCH --job-name=bgz-to-gz
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=5G
#SBATCH --output=job_output/0-bgz.%A.out
#SBATCH --error=job_output/0-bgz.%A.error
#SBATCH --array=1-870

# Move into the directory containing the .bgz files
cd temp


# Get the file name for this job
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" bgz_file_list.txt)


# Rename the file to .gz
mv "$FILE" "${FILE%.bgz}.gz"

# Gunzip the renamed file
gunzip "${FILE%.bgz}.gz"
