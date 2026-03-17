#!/bin/bash

# Usage: ./run_diff_metilene.sh FILE_1 FILE_2 NAME_1 NAME_2 OUT_DIR
if [ "$#" -lt 5 ]; then
    echo "Usage: $0 FILE_1 FILE_2 NAME_1 NAME_2 OUT_DIR"
    exit 1
fi

SAMPLE_1=$1
SAMPLE_2=$2
NAME_1=$3
NAME_2=$4
OUT_DIR=$5

mkdir -p "${OUT_DIR}"

module purge
module load micromamba
eval "$(micromamba shell hook --shell bash)"
micromamba activate metilene

# Constants
MIN_COV=10
MIN_CPGS_PER_DMR=5
MIN_METH_DIFF=0.1
INPUT_METILENE="${OUT_DIR}/input_metilene.txt"

echo "Formatting input files for ${NAME_1} and ${NAME_2}..."

# Use OUT_DIR for temporary files to prevent collisions
TMP_A="${OUT_DIR}/tmp_A.bedgraph"
TMP_V="${OUT_DIR}/tmp_V.bedgraph"
TMP_UNION="${OUT_DIR}/tmp_union.bedgraph"

# Extract Bismark coverage fraction (meth / total) if min_cov is satisfied.
zcat "${SAMPLE_1}" | awk -v mincov="${MIN_COV}" '{if (($5+$6) >= mincov) print $1"\t"$2-1"\t"$2"\t"$4/100}' | sort -k1,1 -k2,2n > "${TMP_A}"
zcat "${SAMPLE_2}" | awk -v mincov="${MIN_COV}" '{if (($5+$6) >= mincov) print $1"\t"$2-1"\t"$2"\t"$4/100}' | sort -k1,1 -k2,2n > "${TMP_V}"

# Use bedtools unionbedg
bedtools unionbedg -i "${TMP_A}" "${TMP_V}" -names g1_A g2_V -filler NA > "${TMP_UNION}"

# Prep metilene input
# Metilene header needs to represent the groups. We use placeholders g1_A and g2_V 
# but metilene will just see them as two groups.
echo -e "chr\tpos\tg1_A\tg2_V" > "${INPUT_METILENE}"
awk '{if ($4 != "NA" && $5 != "NA") print $1"\t"$3"\t"$4"\t"$5}' "${TMP_UNION}" >> "${INPUT_METILENE}"

echo "Running metilene for DMRs..."
metilene -m "${MIN_CPGS_PER_DMR}" -d "${MIN_METH_DIFF}" -a g1_A -b g2_V "${INPUT_METILENE}" | sort -V -k1,1 -k2,2n > "${OUT_DIR}/metilene_dmr_results_10p.bed"

echo "Running metilene for DMCs..."
metilene -f 3 -d "${MIN_METH_DIFF}" -a g1_A -b g2_V "${INPUT_METILENE}" | sort -V -k1,1 -k2,2n > "${OUT_DIR}/metilene_dmc_results_10p.bed"

# Cleanup temporary files
rm "${TMP_A}" "${TMP_V}" "${TMP_UNION}"

echo "Metilene analysis for ${NAME_1} vs ${NAME_2} complete."