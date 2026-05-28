#!/usr/bin/env bash
# ==============================================================================
# Script: generate_deeptools_heatmap.sh
# Description: Generates deepTools heatmaps of DNA methylation around CpG islands
#              using bigWig files in the 'bw' folder.
# Author: Antigravity
# Date: 2026-05-24
# ==============================================================================

# Job name
#SBATCH --job-name=deeptools
# Output and error files (jobid will be appended)
#SBATCH --output=logs/deeptools_%A_%a.out
#SBATCH --error=logs/deeptools_%A_%a.err
# Time and resources - tune as needed
# SBATCH --array=1-15
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G

set -euo pipefail

# ==============================================================================
# 1. CONFIGURATION PARAMETERS (Modify as needed)
# ==============================================================================

# Directory containing BigWig (.bw) files
BW_DIR="../bw"

# UCSC CpG Island annotation file (gzipped)
CPG_GZ="/home/users/allstaff/yan.a/davidson_longread/yan.a/XZ/emseq/dna_meth_collab/data_chr1/raw/cpgIslandExt_hg38.txt.gz"

# Output files
OUT_DIR="."
BED_OUT="${OUT_DIR}/hg38_cpg_islands.bed"
MATRIX_OUT="${OUT_DIR}/cpg_islands_matrix.gz"
HEATMAP_OUT="${OUT_DIR}/cpg_islands_heatmap.pdf"  # Supported formats: .pdf, .png, .svg, etc.

# Genomic scope for analysis:
# "chr1"      -> Analyze chr1 only (fast, matching the chr1 project subsets)
# "standard"  -> Analyze chr1-22, chrX, chrY, chrM
# "genome"    -> Analyze all contigs/chromosomes present in the annotation
SCOPE="genome"

# deepTools computeMatrix mode:
# "reference-point" -> Plot around a specific point (e.g., center, start, end)
# "scale-regions"    -> Scale the body of CpG islands to a fixed length
MATRIX_MODE="reference-point"

# Flanking regions upstream and downstream of CpG islands (in base pairs)
UPSTREAM=2000
DOWNSTREAM=2000

# Mode-specific parameters
# For "reference-point" mode: center, TSS, TES (center is standard for CpG islands)
REF_POINT="center"
# For "scale-regions" mode: length to scale the CpG island body to
REGION_BODY_LEN=1000

# deepTools plotHeatmap visual parameters
COLOR_MAP="RdYlBu_r"    # Reversed RdYlBu maps low=Blue, high=Red. Other options: coolwarm, bwr, viridis
Z_MIN=0                 # Minimum value for colorbar (DNA methylation percentage: 0)
Z_MAX=100               # Maximum value for colorbar (DNA methylation percentage: 100)
KMEANS=1                # Set > 1 to cluster regions (e.g. 3) to see distinct methylation patterns
WHAT_TO_SHOW="plot, heatmap and colorbar"  # Options: "plot, heatmap and colorbar", "heatmap and colorbar", "plot and heatmap"
BIN_SIZE=50             # Bin size for computeMatrix to handle sparse WGBS data (recommend 50)
MISSING_COLOR="white"   # Color for missing data (NaNs) in plotHeatmap
INTERPOLATION="nearest" # Set interpolation to nearest to prevent sparse data from blurring
THREADS=16               # Number of CPU threads to use for computeMatrix

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

# ==============================================================================
# 3. ENVIRONMENT SETUP
# ==============================================================================

log "Setting up environment..."

# Check and load deepTools module if available on cluster
if ! command -v computeMatrix &> /dev/null; then
    log "deepTools not found in PATH. Attempting to load system module..."
    if command -v module &> /dev/null; then
        module load deeptools/3.5.5 || module load deeptools || true
    fi
fi

# Double check if deepTools is loaded successfully
if ! command -v computeMatrix &> /dev/null; then
    log "Error: deepTools (computeMatrix) is not available. Please load the appropriate environment or module first."
    exit 1
fi
log "Using deepTools: $(computeMatrix --version)"

# Check and load bedtools module if available
if ! command -v bedtools &> /dev/null; then
    log "bedtools not found in PATH. Attempting to load system module..."
    if command -v module &> /dev/null; then
        module load bedtools || true
    fi
fi

# ==============================================================================
# 4. PREPARE CpG ISLAND BED FILE
# ==============================================================================

if [[ ! -f "$CPG_GZ" ]]; then
    log "Error: CpG Island annotation file not found at: $CPG_GZ"
    exit 1
fi

log "Preparing CpG Island BED file (Scope: ${SCOPE})..."

# Standard human chromosomes regex (chr1-22, chrX, chrY, chrM)
STANDARD_CHROM_REGEX='^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$'

case "$SCOPE" in
    chr1)
        log "Filtering for chr1 only..."
        gzip -cd "$CPG_GZ" \
            | awk 'BEGIN{OFS="\t"} $2 == "chr1" {print $2, $3, $4, "CGI_" ++n "_" $8 "CpGs", $7, "."}' \
            > "$BED_OUT"
        ;;
    standard)
        log "Filtering for standard human chromosomes..."
        gzip -cd "$CPG_GZ" \
            | awk -v reg="$STANDARD_CHROM_REGEX" 'BEGIN{OFS="\t"} $2 ~ reg {print $2, $3, $4, "CGI_" ++n "_" $8 "CpGs", $7, "."}' \
            > "$BED_OUT"
        ;;
    genome)
        log "Using all contigs/chromosomes..."
        gzip -cd "$CPG_GZ" \
            | awk 'BEGIN{OFS="\t"} {print $2, $3, $4, "CGI_" ++n "_" $8 "CpGs", $7, "."}' \
            > "$BED_OUT"
        ;;
    *)
        log "Error: Invalid SCOPE '$SCOPE'. Must be 'chr1', 'standard', or 'genome'."
        exit 1
        ;;
esac

log "Sorting CpG Island BED file..."
if command -v bedtools &> /dev/null; then
    bedtools sort -i "$BED_OUT" > "${BED_OUT}.tmp" && mv "${BED_OUT}.tmp" "$BED_OUT"
else
    log "bedtools not found. Falling back to system sort..."
    sort -k1,1 -k2,2n "$BED_OUT" > "${BED_OUT}.tmp" && mv "${BED_OUT}.tmp" "$BED_OUT"
fi

log "Saved $(wc -l < "$BED_OUT") CpG Island regions to $BED_OUT"

# ==============================================================================
# 5. SCAN AND CHECK BIGWIG FILES
# ==============================================================================

log "Scanning for BigWig files in $BW_DIR..."

if [[ ! -d "$BW_DIR" ]]; then
    log "Error: BigWig directory not found: $BW_DIR"
    exit 1
fi

BW_FILES=()
BW_LABELS=()

# Using nullglob to avoid literal expansion when no files exist
shopt -s nullglob
for f in "$BW_DIR"/*.bw; do
    # Only include non-empty files (deeptools fails on empty files)
    if [[ -s "$f" ]]; then
        BW_FILES+=("$f")
        # Extract a clean label: removes Bismark-specific naming details and extension
        label=$(basename "$f" | sed -E 's/_1_val_1_bismark_bt2_pe\.deduplicated\.bedGraph\.bw//; s/\.bw//')
        BW_LABELS+=("$label")
    else
        log "Warning: Skipping empty/missing BigWig file: $(basename "$f")"
    fi
done
shopt -u nullglob

if [[ ${#BW_FILES[@]} -eq 0 ]]; then
    log "Error: No non-empty BigWig (.bw) files found in $BW_DIR!"
    exit 1
fi

log "Found ${#BW_FILES[@]} valid BigWig file(s) for plotting:"
for i in "${!BW_FILES[@]}"; do
    log "  [$((i+1))] ${BW_FILES[$i]} -> Label: ${BW_LABELS[$i]}"
done

# ==============================================================================
# 6. RUN DEEPTOOLS COMPUTEMATRIX
# ==============================================================================

log "Running computeMatrix (Mode: ${MATRIX_MODE})..."

if [[ "$MATRIX_MODE" == "reference-point" ]]; then
    computeMatrix reference-point \
        --referencePoint "$REF_POINT" \
        -b "$UPSTREAM" -a "$DOWNSTREAM" \
        -bs "$BIN_SIZE" \
        -R "$BED_OUT" \
        -S "${BW_FILES[@]}" \
        --samplesLabel "${BW_LABELS[@]}" \
        -p "$THREADS" \
        -o "$MATRIX_OUT"
elif [[ "$MATRIX_MODE" == "scale-regions" ]]; then
    computeMatrix scale-regions \
        -b "$UPSTREAM" -a "$DOWNSTREAM" \
        --regionBodyLength "$REGION_BODY_LEN" \
        -bs "$BIN_SIZE" \
        -R "$BED_OUT" \
        -S "${BW_FILES[@]}" \
        --samplesLabel "${BW_LABELS[@]}" \
        -p "$THREADS" \
        -o "$MATRIX_OUT"
else
    log "Error: Invalid MATRIX_MODE '$MATRIX_MODE'. Must be 'reference-point' or 'scale-regions'."
    exit 1
fi

log "Successfully wrote matrix data to $MATRIX_OUT"

# ==============================================================================
# 7. RUN DEEPTOOLS PLOTHEATMAP
# ==============================================================================

log "Running plotHeatmap..."

PLOT_TITLE="DNA Methylation around CpG Islands (Scope: ${SCOPE})"

PLOT_CMD=(
    plotHeatmap
    -m "$MATRIX_OUT"
    -o "$HEATMAP_OUT"
    --colorMap "$COLOR_MAP"
    --zMin "$Z_MIN" --zMax "$Z_MAX"
    --whatToShow "$WHAT_TO_SHOW"
    --missingDataColor "$MISSING_COLOR"
    --interpolationMethod "$INTERPOLATION"
    --plotTitle "$PLOT_TITLE"
)

# Optional clustering
if [[ "$KMEANS" -gt 1 ]]; then
    log "Clustering CpG islands into $KMEANS groups using K-means..."
    PLOT_CMD+=(--kmeans "$KMEANS")
fi

# Execute plotHeatmap
"${PLOT_CMD[@]}"

log "=============================================================================="
log "Success! Heatmap successfully generated and saved to: $HEATMAP_OUT"
log "=============================================================================="
