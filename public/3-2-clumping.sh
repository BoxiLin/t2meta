#!/bin/bash
#SBATCH --job-name=plink_clump_new
#SBATCH --array=1001-1016
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=job_output/plink_clump_new%A_%a.out
#SBATCH --error=job_output/plink_clump_new%A_%a.err

# Load necessary modules (if any)
module load plink/1.90b3x

# # Run once:
# find test_specific_summary -type f -name "*.txt" > files.txt

# Get the file name from the list
FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" files.txt)

# Remove the directory and extension from the file name to use it as output
BASENAME=$(basename $FILE .txt)
plink --bfile EUR_phase3_nodup --clump $FILE --clump-p1 5e-8 --clump-p2 5e-8 --clump-r2 0.01 --clump-kb 10000 --clump-field P --clump-snp-field variant --out locus_out_100kb/$BASENAME

