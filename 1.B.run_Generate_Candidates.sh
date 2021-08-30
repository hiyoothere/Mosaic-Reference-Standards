#!/bin/bash


PROJ_DIR="$HOME"
SCRIPT_DIR="$PWD"
DATA_DIR="${PROJ_DIR}/Sect1.B/1.germline_gt_res"
OUT_DIR="${PROJ_DIR}/Sect1.B/2.ctrl_cand_preprocessing"
OUT_SUB_ARR=("a.nonvar_preprocessing" "b.var_preprocessing" "c.comb_cell_line_gt_res" "d.comb_tool_gt_res")
TOOL_ARR=("sk" "dv")
CTRL_ARR=("pc" "nc_wt" "nc_gl")
CHR_LIST=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")


echo "# Step 1 & 2"
for tool_id in "${TOOL_ARR[@]}"; do
    # Step 1: Extract non-variant loci with passing filters and not overlapping with variant loci
    mkdir -p "${OUT_DIR}/${OUT_SUB_ARR[0]}/${tool_id}"
    in_arr=($(ls ${DATA_DIR}/${tool_id}/*genome*))
    for in_file in "${in_arr[@]}"; do
        bash "${SCRIPT_DIR}/1.B.1.pipe_Extract_Nonvar_Loci.sh" "${in_file}"  \
        "${OUT_DIR}/${OUT_SUB_ARR[0]}/${tool_id}"
    done
    # Step 2: Extract variant loci with passing filters and diploid biallelic genotype
    mkdir -p "${OUT_DIR}/${OUT_SUB_ARR[1]}/${tool_id}"
    in_arr=($(ls ${DATA_DIR}/${tool_id}/*variants*))
    for in_file in "${in_arr[@]}"; do
        python3 "${SCRIPT_DIR}/1.B.2.pipe_Extract_Var_Loci.py" "${in_file}"  \
        "${OUT_DIR}/${OUT_SUB_ARR[1]}/${tool_id}"
    done
done


# Step 3: Prepare extracted variant/non-variant information for combined genotyping 
echo "# Step 3"
subsub_arr=("log" "pickle" "by_chr")
for tool_id in "${TOOL_ARR[@]}"; do
    for subsub_dir in "${subsub_arr[@]}"; do
        mkdir -p "${OUT_DIR}/${OUT_SUB_ARR[2]}/${tool_id}/${subsub_dir}"
    done
    in_arr=($(ls ${OUT_DIR}/${OUT_SUB_ARR[1]}/${tool_id}/*vcf))
    for in_file in "${in_arr[@]}"; do
        qsub -j "y" -o "${OUT_DIR}/${OUT_SUB_ARR[2]}/${tool_id}/${subsub_arr[0]}" \
        -N "pickle.${tool_id}.${in_file##*/}" \
        "${SCRIPT_DIR}/../wrapper/PyWrapper2.sh" \
        "${SCRIPT_DIR}/1.B.3.pipe_Pickle_VCF_Or_BED.py" \
        "${in_file}" "${OUT_DIR}/${OUT_SUB_ARR[2]}/${tool_id}/${subsub_arr[1]}"
        sleep 1
    done
    in_arr=($(ls ${OUT_DIR}/${OUT_SUB_ARR[0]}/${tool_id}/*pass.wo_var.bed))
    for in_file in "${in_arr[@]}"; do
        qsub -j "y" -o "${OUT_DIR}/${OUT_SUB_ARR[2]}/${tool_id}/${subsub_arr[0]}" \
        -N "pickle.${tool_id}.${in_file##*/}" \
        "${SCRIPT_DIR}/../wrapper/PyWrapper2.sh" \
        "${SCRIPT_DIR}/1.B.3.pipe_Pickle_VCF_Or_BED.py" \
        "${in_file}" "${OUT_DIR}/${OUT_SUB_ARR[2]}/${tool_id}/${subsub_arr[1]}"
        sleep 1
    done
done


# Step 4: Identify mutually exclusive variant or unanimous non-variant loci
echo "# Step 4"
cmb_ctrl_type=("nc_wt" "pc")
for tool_id in "${TOOL_ARR[@]}"; do
    for cur_chr in "${CHR_LIST[@]}"; do
        qsub -j "y" -o "${OUT_DIR}/${OUT_SUB_ARR[2]}/${tool_id}/${subsub_arr[0]}" \
        -N "cmb_cl.${tool_id}.chr${cur_chr}" -hold_jid "pickle.${tool_id}*" \
        "${SCRIPT_DIR}/../wrapper/PyWrapper4.sh" \
        "${SCRIPT_DIR}/1.B.4-a.pipe_Combined_Genotyping_By_Chr.py" \
        "${OUT_DIR}/${OUT_SUB_ARR[2]}/${tool_id}/${subsub_arr[1]}" "${cur_chr}" \
        "${OUT_DIR}/${OUT_SUB_ARR[2]}/${tool_id}/${subsub_arr[2]}"
        sleep 1
    done
    for ctrl_type in "${cmb_ctrl_type[@]}"; do
        qsub -j "y" -o "${OUT_DIR}/${OUT_SUB_ARR[2]}/${tool_id}/${subsub_arr[0]}" \
        -N "cmb_cl.${tool_id}.merge.${ctrl_type}" -hold_jid "cmb_cl.${tool_id}.chr*" \
        "${SCRIPT_DIR}/../wrapper/PyWrapper3.sh" \
        "${SCRIPT_DIR}/1.B.4-b.pipe_Merge_Combined_Genotyping_For_All_Chr.py" \
        "${OUT_DIR}/${OUT_SUB_ARR[2]}/${tool_id}/${subsub_arr[2]}" "${ctrl_type}" \
        "${OUT_DIR}/${OUT_SUB_ARR[2]}/${tool_id}/${ctrl_type}.vcf" \
        sleep 1
    done
done


# Please proceed to step 5 after step 4 has successfully finished
sleep 10000


# Step 5: Generate control candidates
echo "# Step 5"
for subsub_dir in "${CTRL_ARR[@]}"; do
    mkdir -p "${OUT_DIR}/${OUT_SUB_ARR[3]}/${subsub_dir}"
done

# Step 5-A: Generate positive control cadidates
echo "# Step 5-A: PC"
sk_pc_file="${OUT_DIR}/${OUT_SUB_ARR[2]}/sk/pc.vcf"
dv_pc_file="${OUT_DIR}/${OUT_SUB_ARR[2]}/dv/pc.vcf"
bash "${SCRIPT_DIR}/1.B.5-a.pipe_Generate_PC_Cand.sh" \
"${sk_pc_file}" "${dv_pc_file}" "${OUT_DIR}/${OUT_SUB_ARR[3]}/${CTRL_ARR[0]}"

# Step 5-B: Generate wildtype negative control candidates
echo "# Step 5-B: NC WT"
sk_nc_wt_file="${OUT_DIR}/${OUT_SUB_ARR[2]}/sk/nc_wt.vcf"
dv_nc_wt_file="${OUT_DIR}/${OUT_SUB_ARR[2]}/dv/nc_wt.vcf"
python3 "${SCRIPT_DIR}/1.B.5.Merge_Combined_Genotyping_For_All_Tools.py" \
"${sk_nc_wt_file}" "${dv_nc_wt_file}" \
"${OUT_DIR}/${OUT_SUB_ARR[3]}/${CTRL_ARR[1]}" "sk_dv.conc"

# Step 5-C: Generate germline negative control candidates
echo "# Step 5-C: NC GL"
sk_s0_pass_file="${OUT_DIR}/${OUT_SUB_ARR[1]}/sk/S0.variants.pass.vcf"
dv_s0_pass_file="${OUT_DIR}/${OUT_SUB_ARR[1]}/dv/S0.variants.pass.vcf"
bash "${SCRIPT_DIR}/1.B.5-c.pipe_Generate_NC_GL_Cand.sh" \
"${sk_s0_pass_file}" "${dv_s0_pass_file}" "${OUT_DIR}/${OUT_SUB_ARR[3]}/${CTRL_ARR[2]}"
