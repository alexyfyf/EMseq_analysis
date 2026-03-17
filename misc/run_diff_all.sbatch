#!/bin/bash
#SBATCH --job-name=diff_meth
#SBATCH --output=/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/diff_meth_%j.out
#SBATCH --error=/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/diff_meth_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16
#SBATCH --partition=regular

module purge
module load R/4.4.2

# echo "Running DSS..."
# Rscript /home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/run_diff_dss.R

echo "Running BSmooth..."
Rscript /home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/run_diff_bsmooth.R

# echo "Running MethylKit..."
# Rscript /home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/run_diff_methylkit.R

echo "All complete!"
