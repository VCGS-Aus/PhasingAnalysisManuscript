#!/usr/bin/env bash

# Usage extract_phase_sets_to_calculate_phasing_accuracy.sh manifest.txt out_dir
# manifest.txt: sample_id replicate_vcf truth_vcf

set -euo pipefail

INPUT_LIST=$1
OUTPUT_DIR=$2

while IFS=$'\t' read -r SAMPLE_ID SAMPLE_VCF TRUTH_VCF; do
    echo "Processing ${SAMPLE_ID}..."

    OUT_VCF="${OUTPUT_DIR}/common_variants_${SAMPLE_ID}_and_truth.vcf.gz"
    OUT_TSV="${OUTPUT_DIR}/ps_common_variants_${SAMPLE_ID}_and_truth.tsv"
    TRUTH_SAMPLE_ID=$(basename "$TRUTH_VCF" .vcf.gz)
    TRUTH_TSV="${OUTPUT_DIR}/ps_${TRUTH_SAMPLE_ID}_truth.tsv"

    # Step 1: intersection
    /misc/vcgs/seq/cpipe-wgs-staging/tools/bcftools/1.23/bcftools isec -n=2 -w1 \
        "${SAMPLE_VCF}" \
        "${TRUTH_VCF}" \
        -Oz -o "${OUT_VCF}"

    # Index (recommended for bcftools query)
    /misc/vcgs/seq/cpipe-wgs-staging/tools/bcftools/1.23/bcftools index -t "${OUT_VCF}"

    # Step 2: query + filter
    /misc/vcgs/seq/cpipe-wgs-staging/tools/bcftools/1.23/bcftools query -f '%CHROM\t%POS\t%FORMAT\t[%GT\t%PS]\n' "${OUT_VCF}" \
        | awk -F'\t' '$6 != "."' > "${OUT_TSV}"

    # Step 3: truth ps
    if [ -s $TRUTH_TSV ]; then
	echo "ps_${SAMPLE_ID}_truth.tsv exists and is not empty"
    else
	/misc/vcgs/seq/cpipe-wgs-staging/tools/bcftools/1.23/bcftools query -f '%CHROM\t%POS\t%FORMAT\t[%GT\t%PS]\n' "${TRUTH_VCF}" | awk -F'\t' '$6 != "."' > "${TRUTH_TSV}"
    fi

    echo "Finished ${SAMPLE_ID}"
done < "${INPUT_LIST}"
