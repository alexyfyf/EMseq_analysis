#!/bin/bash
#SBATCH --job-name=plot_meth
#SBATCH --output=/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/plot_meth_%j.out
#SBATCH --error=/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/plot_meth_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=regular

module purge
module load R/4.4.2

echo "Running Visualization Script..."
Rscript /home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/collate_and_plot_meth.R
echo "All complete!"
