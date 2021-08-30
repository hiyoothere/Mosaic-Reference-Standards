#!/bin/bash


## Input Parameters
SK_PASS_FILE=$1
DV_PASS_FILE=$2
OUT_DIR=$3
MIX_SAMP_EXT_REG=$4

# Global Variables
SCRIPT_DIR="$PWD"

echo "## Step 5-C-1"
python3 "${SCRIPT_DIR}/1.B.5.pipe_Merge_Combined_Genotyping_For_All_Tools.py" \
"${SK_PASS_FILE}" "${DV_PASS_FILE}" "${OUT_DIR}" "sk_dv.conc"

echo "## Step 5-C-2"
bedtools --version
in_filename="${SK_PASS_FILE##*/}"
conc_file="${OUT_DIR}/sk_dv.conc.${in_filename}"
sk_fail_file="${SK_PASS_FILE::-9}.fail.vcf"
dv_fail_file="${DV_PASS_FILE::-9}.fail.vcf"
temp_file="${conc_file::-4}.temp.vcf"
out_file="${conc_file::-4}.fail_rem.vcf"
bedtools intersect -v -header -a "${conc_file}" -b "${sk_fail_file}" > "${temp_file}"
bedtools intersect -v -header -a "${temp_file}" -b "${dv_fail_file}" > "${out_file}"

echo "## Step 5-C-3"
bedtools --version
in_file="${out_file}"
out_file="${in_file::-4}.mx_ext_reg_rem.vcf"
bedtools intersect -v -header -a "${in_file}" -b "${MIX_SAMP_EXT_REG}" > "${out_file}"