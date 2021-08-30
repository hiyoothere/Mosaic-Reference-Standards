#!/bin/bash


## Input Parameters
IN_FILE=$1
OUT_DIR=$2

# Global Variables
SCRIPT_DIR="$PWD"

echo "## Step 1-1"
python3 "${SCRIPT_DIR}/1.B.1-a.pipe_Expand_Passing_Nonvar_Loci.py" "${IN_FILE}" "${OUT_DIR}"

echo "## Step 1-2"
bedtools --version
in_head="${OUT_DIR}/${IN_FILE##*/}"
if [[ "${IN_FILE##*.}" == "gz" ]]; then
    vcf_file="${IN_FILE::-14}.variants.vcf.gz"
    in_head="${in_head::-7}"
else
    vcf_file="${IN_FILE::-11}.variants.vcf"
    in_head="${in_head::-4}"
fi
in_file="${in_head}.pass.vcf"
pass_file="${in_file::-4}.wo_var.vcf"
fail_file="${in_file::-4}.w_var.vcf"
bedtools intersect -v -header -a "${in_file}" -b "${vcf_file}" > "${pass_file}"
bedtools intersect -header -a "${in_file}" -b "${vcf_file}" > "${fail_file}"

echo "## Step 1-3"
convert2bed --version
bedtools --version
in_file="${pass_file}"
out1_file="${in_file::-4}.temp.bed"
out2_file="${in_file::-4}.bed"
cat ${in_file} | convert2bed --input=vcf > "${out1_file}"
bedtools merge -i "${out1_file}" > "${out2_file}"
