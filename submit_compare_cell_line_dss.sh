#!/bin/bash
#SBATCH --job-name=compare_dss
#SBATCH --output=/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/logs/compare_dss_%j.out
#SBATCH --error=/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/logs/compare_dss_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=regular

# Load R module
module purge
module load R/4.4.2

BASE_DIR="/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq"
SCRIPT="${BASE_DIR}/data/compare_cell_line_dss.R"
OUT_ROOT="${BASE_DIR}/output/compare_dss"

mkdir -p "${OUT_ROOT}"

# Run for each cell line
CELL_LINES=("CB" "KASU" "MDSL" "MOLM" "Pt")

for cell in "${CELL_LINES[@]}"; do
    echo "Processing ${cell}..."
    Rscript "${SCRIPT}" "${cell}" "${OUT_ROOT}"
done

echo "All comparisons complete."
