#!/bin/bash


## Input Parameters
SK_CAND_FILE=$1
DV_CAND_FILE=$2
OUT_DIR=$3

# Global Variables
SCRIPT_DIR="$PWD"
CNV_FILE="$HOME/CNVkit/6.CNVkit.0.3.s.m.bed"
MAN_IND_PASS_FILE="$HOME/ind_man_check/prev.ind.manpass.total.vcf"
MAN_IND_FAIL_FILE="$HOME/ind_man_check/prev.ind.manfail.total.vcf"

echo "## Step 5-A-1"
bedtools --version
sk_uniq_file="${OUT_DIR}/sk.uniq.${SK_CAND_FILE##*/}"
dv_uniq_file="${OUT_DIR}/dv.uniq.${DV_CAND_FILE##*/}"
conc_file="${OUT_DIR}/sk_dv.conc.${SK_CAND_FILE##*/}"
python3 "${SCRIPT_DIR}/1.B.5.pipe_Merge_Combined_Genotyping_For_All_Tools.py" \
"${SK_CAND_FILE}" "${DV_CAND_FILE}" "${OUT_DIR}" "sk_dv.conc"
bedtools intersect -v -header -a "${SK_CAND_FILE}" -b "${DV_CAND_FILE}" > "${sk_uniq_file}"
bedtools intersect -v -header -a "${DV_CAND_FILE}" -b "${SK_CAND_FILE}" > "${dv_uniq_file}"

echo "## Step 5-A-2"
bedtools --version
bedtools intersect -v -header -a "${conc_file}" -b "${CNV_FILE}" > "${conc_file::-4}.cnv_rem.vcf"
bedtools intersect -v -header -a "${sk_uniq_file}" -b "${CNV_FILE}" > "${sk_uniq_file::-4}.cnv_rem.vcf"
bedtools intersect -v -header -a "${dv_uniq_file}" -b "${CNV_FILE}" > "${dv_uniq_file::-4}.cnv_rem.vcf"

echo "## Step 5-A-3"
bedtools --version
sk_uniq_file="${sk_uniq_file::-4}.cnv_rem.vcf"
dv_uniq_file="${dv_uniq_file::-4}.cnv_rem.vcf"
python3 "${SCRIPT_DIR}/1.B.5-a-1.pipe_Extract_Indel_Loci.py" "${sk_uniq_file}"
python3 "${SCRIPT_DIR}/1.B.5-a-1.pipe_Extract_Indel_Loci.py" "${dv_uniq_file}"
sk_ind_file="${sk_uniq_file::-4}.ind.vcf"
dv_ind_file="${dv_uniq_file::-4}.ind.vcf"

bedtools intersect -header -a "${sk_ind_file}" -b "${MAN_IND_PASS_FILE}" > "${sk_ind_file::-4}.man_pass.vcf"
bedtools intersect -header -v -a "${sk_ind_file}" -b "${MAN_IND_PASS_FILE}" > "${sk_ind_file::-4}.man_fail.vcf"
bedtools intersect -header -v -a "${sk_ind_file::-4}.man_fail.vcf" -b "${MAN_IND_FAIL_FILE}" > "${sk_ind_file::-4}.man_fail.need_check.vcf"

bedtools intersect -header -a "${dv_ind_file}" -b "${MAN_IND_PASS_FILE}" > "${dv_ind_file::-4}.man_pass.vcf"
bedtools intersect -header -v -a "${dv_ind_file}" -b "${MAN_IND_PASS_FILE}" > "${dv_ind_file::-4}.man_fail.vcf"
bedtools intersect -header -v -a "${dv_ind_file::-4}.man_fail.vcf" -b "${MAN_IND_FAIL_FILE}" > "${dv_ind_file::-4}.man_fail.need_check.vcf"

echo "## Step 5-A-4"
out_file="${conc_file::-4}.man_uniq_ind.cnv_rem.vcf"
python3 "${SCRIPT_DIR}/1.B.5-a-2.pipe_Merge_Two_VCF.py" \
"${dv_ind_file::-4}.man_pass.vcf" "${sk_ind_file::-4}.man_pass.vcf" "${out_file::-4}.temp1.vcf"
python3 "${SCRIPT_DIR}/1.B.5-a-2.pipe_Merge_Two_VCF.py" \
"${out_file::-4}.temp1.vcf" "${conc_file::-4}.cnv_rem.vcf" "${out_file}"

echo "## Step 5-A-5"
grep -v "^chr[XxYy]" "${out_file}" > "${out_file::-4}.no_sex.vcf"