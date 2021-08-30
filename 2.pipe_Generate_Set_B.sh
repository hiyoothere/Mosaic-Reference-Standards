#!/bin/bash


## Input Parameters
SAMP_ID=$1
BASE_BAM=$2
SRC_BAM=$3
PC_CAND_VCF=$4
REF=$5
SEED=$6
OUT_DIR=$7


## Global Variables
BASE_DEP=1135
TGT_DEP=1100
PICARD="/opt/Yonsei/Picard/2.23.1/picard.jar"
GATK="/data/project/MRS/Resource/gatk-4.1.4.1/gatk"
OUT_SUB_ARR=("a.pc_cand" "b.pc_rem_base" "c.pc_ext_src" "d.merged" "e.iv_ext_reg")
pc_file=$(ls "${OUT_DIR}/${OUT_SUB_ARR[0]}" | egrep "M${INPUT:1:1}.*vcf")
PC_FILE="${OUT_DIR}/${OUT_SUB_ARR[0]}/${pc_file}"


echo "# Step 1: Downsample base (MRC5/S0) bam file"
out_file="${OUT_DIR}/${OUT_SUB_ARR[1]}/${SAMP_ID}.ds_base.${DEP}x.sd${SEED}.bam"
ds_ratio=$(echo "${DEP} / ${REF_DEP}" | bc -l)

if [[ -f "${out_file}" ]]; then
	echo "${out_file} exists"
else
	echo "** Generating ${out_file}"
	echo "****************"
	java -jar "${PICARD}" DownsampleSam \
	I="${BASE_BAM}" \
	O="${out_file}" \
	S="ConstantMemory" \
	P="${ds_ratio}" \
	CREATE_INDEX="true" \
	R="${SEED}"
fi


echo "# Step 2: Remove reads overlapping with positive control candidate loci from downsampled base bam file."
out_head="${OUT_DIR}/${OUT_SUB_ARR[1]}/${SAMP_ID}.ds_base.${DEP}x.sd${SEED}.tp_rem"
bedtools --version
samtools --verison

bedtools intersect -v -a "${out_dir}/${SAMP_ID}.ds_base.${DEP}x.sd${SEED}.bam" -b "${PC_FILE}" > "${out_head}.bam"
samtools sort -o "${out_head}.bam" "${out_head}.s.bam"
samtools index "${out_head}.s.bam"


echo "# Step 3: Extract reads overlapping with positive control candidate loci from sample mixture bam file."
out_head="${OUT_DIR}/${OUT_SUB_ARR[2]}/${SAMP_ID}.tp_ext"
bedtools --version
samtools --verison

bedtools intersect -a "${SRC_BAM}" -b "${PC_FILE}" > "${out_head}.bam"
samtools view -h "${out_head}.bam" | less | sed -e "s/^[^@]/RS/" | samtools view -b > "${out_head}.qt_renamed.bam"
samtools sort -o "${out_head}.qt_renamed.bam" "${out_head}.qt_renamed.s.bam"
samtools index "${out_head}.qt_renamed.s.bam"


echo "# Step 4: Merge the positive control removed base bam with positive control overlapping sample mixture bam."
out_head="${OUT_DIR}/${OUT_SUB_ARR[3]}/${SAMP_ID}.MRS"
samtools --verison

samtools merge -c "${out_head}.bam" \
"${OUT_DIR}/${OUT_SUB_ARR[1]}/${SAMP_ID}.ds_base.${DEP}x.sd${SEED}.tp_rem.s.bam" \
"${OUT_DIR}/${OUT_SUB_ARR[2]}/${SAMP_ID}.tp_ext.qt_renamed.s.bam"

java -jar "${PICARD}" SortSam \
SORT_ORDER=coordinate \
INPUT="${out_head}.bam" \
OUTPUT="${out_head}.s.bam" \
VALIDATION_STRINGENCY=LENIENT \
CREATE_INDEX=true

java -jar "${PICARD}" AddOrReplaceReadGroups \
INPUT="${out_head}.s.bam" \
OUTPUT="${out_head}.s.RGadd.bam" \
SORT_ORDER=coordinate \
RGLB="WES" \
RGPL="il" \
RGPU="MRS" \
RGSM=$INPUT \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=LENIENT

java -jar "${PICARD}" BuildBamIndex \
I="${out_head}.s.RGadd.bam"

"${GATK}" LeftAlignIndels \
-R "$REF" \
-I "${out_head}.s.RGadd.bam" \
-O "${out_head}.s.RGadd.LA.bam"


echo "# Step 5: Generate bed file of the region covered by the reads extracted from sample mixture bam overlapping with postive control loci."
out_head="${OUT_DIR}/${OUT_SUB_ARR[4]}/${SAMP_ID}.TP.ext"
bedtools --version
bedtools intersect -abam "${SRC_BAM}" -b "${PC_CAND_VCF}" > "${out_head}.bam"
bedtools bamtobed -i "${out_head}.bam" > "${out_head}.bed"
bedtools merge -i "${out_head}.bed" > "${out_head}.m.bed"