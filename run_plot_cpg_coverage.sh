#!/bin/bash
#SBATCH --job-name=plot_cpg_coverage
#SBATCH --output=logs/plot_cpg_coverage/%x_%a.out
#SBATCH --error=logs/plot_cpg_coverage/%x_%a.err
#SBATCH --array=1-15
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00

# Set up environment
module load R/4.4.2

# Define directories
PROJECT_DIR="/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq"
DATA_DIR="${PROJECT_DIR}/data"
RESULTS_DIR="${PROJECT_DIR}/results"
OUTPUT_DIR="${RESULTS_DIR}/cpg_coverage_plots"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${DATA_DIR}/logs/plot_cpg_coverage"

# Get the list of all .cov.gz files
FILE_LIST="${OUTPUT_DIR}/file_list.txt"
find "${RESULTS_DIR}" -name "*.cov.gz" ! -name "._*" | sort > "${FILE_LIST}"

# Get the file for this array task
INPUT_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FILE_LIST}")
SAMPLE_NAME=$(basename "${INPUT_FILE}" | sed 's/\.bismark\.cov\.gz//' | sed 's/\.deduplicated//')

echo "Processing sample: ${SAMPLE_NAME}"
echo "Input file: ${INPUT_FILE}"
echo "Output directory: ${OUTPUT_DIR}"

# Run the R script
Rscript "${DATA_DIR}/plot_cpg_region_coverage.R" \
    --input "${INPUT_FILE}" \
    --out-prefix "${OUTPUT_DIR}/${SAMPLE_NAME}"

echo "Finished processing ${SAMPLE_NAME}"
