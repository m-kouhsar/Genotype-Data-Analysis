#!/bin/bash

# PlinkBinaryLiftOver_hg19_to_hg38.sh
#
# Usage:
# ./PlinkBinaryLiftOver_hg19_to_hg38.sh \
#     input_prefix \
#     hg19ToHg38.over.chain.gz \
#     output_prefix
#
# Example:
# ./PlinkBinaryLiftOver_hg19_to_hg38.sh \
#     cohort \
#     hg19ToHg38.over.chain.gz \
#     cohort_hg38


set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "Usage:"
    echo "PlinkBinaryLiftOver input_prefix chain_file ref_genome_fasta_file output_prefix"
    echo "Example:"
    echo "PlinkBinaryLiftover ./Cohort1_hg19 hg19ToHg38.over.chain hg38.fasta ./cohort1_hg38"
    exit 1
fi

INPUT_PREFIX="$1"
CHAIN_FILE="$2"
REF_FASTA="$3"
OUTPUT_PREFIX="$4"

OUT_DIR=$(dirname "$OUTPUT_PREFIX")
mkdir -p "$OUT_DIR"


BASE=$(basename "$INPUT_PREFIX")


echo "========================================"
echo "Step 1: Remove duplicated SNP IDs"
echo "========================================"

plink2 \
    --bfile "$INPUT_PREFIX" \
    --rm-dup force-first \
    --make-bed \
    --out "${OUT_DIR}/${BASE}_nodup"



echo "========================================"
echo "Step 2: Create BED file for liftOver"
echo "========================================"


awk 'BEGIN{OFS="\t"}
{
    chr=$1
    if(chr !~ /^chr/)
        chr="chr"chr

    print chr,$4-1,$4,$2
}' \
"${OUT_DIR}/${BASE}_nodup.bim" \
> "${OUT_DIR}/${BASE}_nodup_bedlike.bed"



echo "========================================"
echo "Step 3: Run UCSC liftOver"
echo "========================================"


liftOver \
"${OUT_DIR}/${BASE}_nodup_bedlike.bed" \
"$CHAIN_FILE" \
"${OUT_DIR}/${BASE}.lifted.bed" \
"${OUT_DIR}/${BASE}.unmapped.bed"



echo "========================================"
echo "Step 4: Create PLINK update map"
echo "========================================"


# PLINK2 update format:
# SNP_ID   POSITION

awk 'BEGIN{OFS="\t"}
{
    chr=$1
    sub(/^chr/,"",chr)

    print $4,$3
}' \
"${OUT_DIR}/${BASE}.lifted.bed" \
> "${OUT_DIR}/new_pos.txt"

echo "========================================"
echo "Step 5: Create PLINK update chr"
echo "========================================"


# PLINK2 update format:
# SNP_ID   CHR

awk 'BEGIN{OFS="\t"}
{
    chr=$1
    sub(/^chr/,"",chr)

    print $4, chr
}' \
"${OUT_DIR}/${BASE}.lifted.bed" \
> "${OUT_DIR}/new_chr.txt"


echo "========================================"
echo "Step 6: Extract lifted SNPs"
echo "========================================"


awk '{print $4}' \
"${OUT_DIR}/${BASE}.lifted.bed" \
> "${OUT_DIR}/lifted_snps.txt"



echo "========================================"
echo "Step 7: Update PLINK coordinates and Sort"
echo "========================================"


plink2 \
    --bfile "${OUT_DIR}/${BASE}_nodup" \
    --extract "${OUT_DIR}/lifted_snps.txt" \
    --update-map "${OUT_DIR}/new_pos.txt" \
    --update-chr "${OUT_DIR}/new_chr.txt" \
    --sort-vars \
    --make-pgen \
    --out "${OUTPUT_PREFIX}_lifted_temp"

echo "============================================================"
echo "Step 8: Check reference alleles using reference fasta file"
echo "============================================================"

plink2 --pfile "${OUTPUT_PREFIX}_lifted_temp" \
       --fa $REF_FASTA \
       --ref-from-fa \
       --make-bed \
       --out "${OUTPUT_PREFIX}_lifted_temp_refCheck"

echo "========================================"
echo "Step 9: Standardise SNP IDs and Create BED"
echo "========================================"


plink2 \
    --bfile "${OUTPUT_PREFIX}_lifted_temp_refCheck" \
    --set-all-var-ids @:#:\$r:\$a \
    --make-bed \
    --out "${OUTPUT_PREFIX}_lifted_refCheck_ChrPos"


echo ""
echo "========================================"
echo "LiftOver completed"
echo "========================================"
echo "Final PLINK files:"
echo "${OUTPUT_PREFIX}_lifted_ChrPos.bed"
echo "${OUTPUT_PREFIX}_lifted_ChrPos.bim"
echo "${OUTPUT_PREFIX}_lifted_ChrPos.fam"
echo ""

echo "Lifted SNPs:"
wc -l "${OUT_DIR}/${BASE}.lifted.bed"

echo "Unmapped SNPs:"
wc -l "${OUT_DIR}/${BASE}.unmapped.bed"
