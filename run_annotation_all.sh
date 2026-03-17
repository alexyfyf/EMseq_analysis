#!/bin/bash
#SBATCH --job-name=anno_meth
#SBATCH --output=/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/logs/anno_%j.out
#SBATCH --error=/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/logs/anno_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --partition=regular

# Script to run annotation for each comparison folder
BASE_DIR="/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq"
DATA_DIR="${BASE_DIR}/data"
OUT_ROOT="${BASE_DIR}/output"

module purge
module load R/4.4.2

# Define cell lines
CELL_LINES=("CB" "KASU" "MDSL" "MOLM" "Pt")

# Define treatment pairs
PAIRS=("A:V" "A:G" "G:V")

for cell in "${CELL_LINES[@]}"; do
    for pair in "${PAIRS[@]}"; do
        trt1="${pair%%:*}"
        trt2="${pair##*:}"
        comp="${cell}_${trt1}_vs_${trt2}"
        out_dir="${OUT_ROOT}/${comp}"
        
        if [ -d "$out_dir" ]; then
            echo "Running annotation for ${comp}..."
            Rscript "${DATA_DIR}/annotate_results.R" "${out_dir}"
        else
            echo "Warning: Output directory ${out_dir} not found. Skipping."
        fi
    done
done

echo "All annotations completed."
