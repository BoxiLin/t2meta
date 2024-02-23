#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH --job-name=3-SUMMARY
#SBATCH --output=job_output/%x_%A_%a.out
#SBATCH --error=job_output/%x_%A_%a.err
#SBATCH --array=1-290


module load R/4.2.2
module load cairo
module load zlib

ARG=$SLURM_ARRAY_TASK_ID

echo $ARG
Rscript --no-save 1-1-QC.R $ARG
Rscript --no-save 1-2-TEST.R $ARG
Rscript --no-save 1-3-SUMMARY.R $ARG
