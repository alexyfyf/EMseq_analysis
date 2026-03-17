#!/bin/bash
#SBATCH --job-name=meth_array
#SBATCH --output=/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/logs/array_%A_%a.out
#SBATCH --error=/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq/data/logs/array_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=512G
#SBATCH --cpus-per-task=16
#SBATCH --partition=regular
#SBATCH --array=0-14

# Master script to run differential methylation comparisons in PARALLEL via Slurm Array.
# Each task (0-14) handles one specific cell line and treatment pairing.

# --- Parameters ---
BASE_DIR="/home/users/allstaff/yan.a/davidson_longread/yan.a/emseq"
DATA_DIR="${BASE_DIR}/data"
RESULTS_DIR="${BASE_DIR}/results"
OUT_ROOT="${BASE_DIR}/output"

# Dry run flag (Set to "false" to actually run the jobs)
DRY_RUN="false"

# Define the comparison matrix
CELL_LINES=("CB:CB" "KASU:Kasu" "MDSL:MDSL" "MOLM:Molm" "Pt:Pt")
PAIRS=("A:V" "A:G" "G:V")

# Flatten matrix into a list for easy indexing
ALL_CONFIGS=()
for cell_entry in "${CELL_LINES[@]}"; do
    for pair_entry in "${PAIRS[@]}"; do
        ALL_CONFIGS+=("${cell_entry}:${pair_entry}")
    done
done

# Select the configuration for THIS specific array task
current_config="${ALL_CONFIGS[$SLURM_ARRAY_TASK_ID]}"

# Parse the components
# Format: CellPrefix:ResSuffix:Trt1:Trt2
cell_prefix=$(echo $current_config | cut -d: -f1)
res_suffix=$(echo $current_config | cut -d: -f2)
trt1=$(echo $current_config | cut -d: -f3)
trt2=$(echo $current_config | cut -d: -f4)

# Define file paths
name1="${cell_prefix}_${trt1}_1"
name2="${cell_prefix}_${trt2}_1"

file1="${RESULTS_DIR}/results_${res_suffix}_${trt1}/bismark/methylation_calls/methylation_coverage/${name1}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"
file2="${RESULTS_DIR}/results_${res_suffix}_${trt2}/bismark/methylation_calls/methylation_coverage/${name2}_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"

comparison_name="${cell_prefix}_${trt1}_vs_${trt2}"
out_dir="${OUT_ROOT}/${comparison_name}"

# Setup logs directory
mkdir -p "${DATA_DIR}/logs"

echo "=========================================================="
echo "Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Comparison:    ${comparison_name}"
echo "Group 1:       ${name1} (${file1})"
echo "Group 2:       ${name2} (${file2})"
echo "Output:        ${out_dir}"
echo "=========================================================="

# Load R module
module purge
module load R/4.4.2

if [ "$DRY_RUN" = "true" ]; then
    echo "[DRY RUN] Would create directory: ${out_dir}"
    echo "[DRY RUN] Would run MethylKit: Rscript ${DATA_DIR}/run_diff_methylkit.R ${file1} ${file2} ${name1} ${name2} ${out_dir}"
    echo "[DRY RUN] Would run BSmooth:   Rscript ${DATA_DIR}/run_diff_bsmooth.R ${file1} ${file2} ${name1} ${name2} ${out_dir}"
    echo "[DRY RUN] Would run DSS:       Rscript ${DATA_DIR}/run_diff_dss.R ${file1} ${file2} ${name1} ${name2} ${out_dir}"
    echo "[DRY RUN] Would run Metilene:  ${DATA_DIR}/run_diff_metilene.sh ${file1} ${file2} ${name1} ${name2} ${out_dir}"
    echo "[DRY RUN] Would run Plotting:  Rscript ${DATA_DIR}/collate_and_plot_meth.R ${out_dir}"
else
    mkdir -p "${out_dir}"
    
    echo "Running MethylKit..."
    Rscript "${DATA_DIR}/run_diff_methylkit.R" "${file1}" "${file2}" "${name1}" "${name2}" "${out_dir}"
    
    echo "Running BSmooth..."
    Rscript "${DATA_DIR}/run_diff_bsmooth.R" "${file1}" "${file2}" "${name1}" "${name2}" "${out_dir}"
    
    echo "Running DSS..."
    Rscript "${DATA_DIR}/run_diff_dss.R" "${file1}" "${file2}" "${name1}" "${name2}" "${out_dir}"
    
    echo "Running Metilene..."
    "${DATA_DIR}/run_diff_metilene.sh" "${file1}" "${file2}" "${name1}" "${name2}" "${out_dir}"
    
    echo "Running visualization (collate_and_plot_meth.R)..."
    Rscript "${DATA_DIR}/collate_and_plot_meth.R" "${out_dir}"
fi

echo "Task ${SLURM_ARRAY_TASK_ID} finished (DRY_RUN=${DRY_RUN})."
