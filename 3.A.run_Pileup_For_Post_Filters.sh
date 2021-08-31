#!/bin/bash

PROJ_DIR="$HOME"
SCRIPT_DIR="$PWD"
LOG_DIR="$PWD/log"
DATA_DIR="${PROJ_DIR}/samples"
OUT_DIR="${PROJ_DIR}/Sect3.A"
OUT_SUB_ARR=("a.pileup" "b.parsed_raw" "c.parsed_by_tag" "d.parsed_all_samp" "p.py_pickle")
REF="/data/resource/reference/human/NCBI/GRCh38_GATK/BWAIndex/genome.fa"
CHR_LIST=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y")


## Create Input
create_one_single() {
    # Example input: "1-3_6"
    mos_var_str=$1
    mos=$(echo $mos_var_str | cut -d "-" -f 1)
    var=$(echo $mos_var_str | cut -d "-" -f 2)
    read -a var_range < <(echo $var | tr '_' ' ')
    ret=""
    for (( i = ${var_range[0]}; i <= ${var_range[1]}; i++ )); do
        ret+="M${mos}-${i};"
    done
    echo $ret >&2
    echo $ret
}

create_single() {
    # Example input: "1 2 3"
    mos_str=$1
    read -a mos_arr <<< $mos_str
    ret=""
    for mos in "${mos_arr[@]}"; do
        case $mos in
            1 )
                ret+=$(create_one_single "1-1_9");;
            2 )
                ret+=$(create_one_single "2-1_12");;
            3 )
                ret+=$(create_one_single "3-1_18");;
        esac
    done
    echo $ret >&2
    echo $ret
}

read -a input_arr < <(create_single "1 2 3" | tr ';' ' ')


## Run
mkdir -p "$LOG_DIR"
mkdir -p "${OUT_DIR}/${OUT_SUB_ARR[0]}"
for i in ${!input_arr[*]}; do
    job_id="PU.${input_arr[$i]}"
    in_bam="${DATA_DIR}/${input_arr[$i]}.bam"
    out_pup="${OUT_DIR}/${OUT_SUB_ARR[0]}/${input_arr[$i]}.pup"

	qsub -j "y" -o $LOG_DIR -N $job_id "${SCRIPT_DIR}/3.A.1.pipe_Pileup.sh" \
	"${in_bam}" "${out_pup}" "${REF}"
	sleep 1

    qsub -j "y" -o $LOG_DIR -hold_jid "${job_id}" -N "${job_id}.sc" \
    "${SCRIPT_DIR}/wrapper/PyWrapper1.sh"\
    "${SCRIPT_DIR}/3.A.2.pipe_Split_Pileup_By_Chr.py" "${out_pup}"
    sleep 1
done


out_dir="${OUT_DIR}/${OUT_SUB_ARR[1]}"
mkdir -p out_dir="${OUT_DIR}/${OUT_SUB_ARR[2]}"
for cur_input in ${input_arr[*]}; do
    job_id="PU.${input_arr[$i]}"
    mkdir -p "${out_dir}/by_chr_${cur_input}"
    
    for cur_chr in ${CHR_LIST[*]}; do
        cur_pup="${OUT_DIR}/${OUT_SUB_ARR[0]}/${input_arr[$i]}.chr${cur_chr}.pup"

        qsub -j "y" -o $LOG_DIR -hold_jid "${job_id}.sc" -N "${job_id}.parse_${cur_chr}" \
        "${SCRIPT_DIR}/wrapper/PyWrapper4.sh" \
        "${SCRIPT_DIR}/3.A.3.pipe_Parse_Pileup_By_Chr.py" \
        "${OUT_DIR}/${OUT_SUB_ARR[4]}" "${cur_pup}" "${out_dir}/by_chr_${cur_input}" "false"
        sleep 1
    done

    qsub -j "y" -o $LOG_DIR -hold_jid "${job_id}.parse_${cur_chr}*" -N "${job_id}.merge_chr" \
    "${SCRIPT_DIR}/wrapper/PyWrapper3.sh" \
    "${SCRIPT_DIR}/3.A.4.pipe_Merge_Pileup_For_All_Chr.py" \
    "${out_dir}" "${cur_input}" "ND"
    sleep 1
    
    qsub -j "y" -o $LOG_DIR -hold_jid "${job_id}.merge_chr" -N "${job_id}.split_tag" \
    "${SCRIPT_DIR}/wrapper/PyWrapper5.sh" \
    "${SCRIPT_DIR}/3.A.5.pipe_Split_Pileup_By_Tag.py" \
    "${out_dir}" "${cur_input}" "ND" "${OUT_DIR}/${OUT_SUB_ARR[2]}" "false"
    sleep 1
done


input_arr=("nc_low" "nc_snv" "nc_ind" "pc_snv" "pc_ind")
for cur_input in ${input_arr[*]}; do
    job_id="0PU.merge_${input_arr[$i]}"
    qsub -j "y" -o $LOG_DIR -hold_jid "^PU*" -N "${job_id}" \
    "${SCRIPT_DIR}/wrapper/PyWrapper4.sh" \
    "${SCRIPT_DIR}/3.A.6.pipe_Merge_Pileup_For_All_Samp.py" \
    "${OUT_DIR}/${OUT_SUB_ARR[2]}" "${cur_input}" "${OUT_DIR}/${OUT_SUB_ARR[3]}" "all_samp_mer"
    sleep 1
done