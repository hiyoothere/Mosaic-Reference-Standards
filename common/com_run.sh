#!/bin/bash

## Key Variables
# ROOT_DIR="/data/project/MRS"
# OUT_DIR="${ROOT_DIR}/4.Analysis_Mosaic"

## Function Collection
create_one_single() {
    # 1-3_6
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
    # 1 2 3
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

create_one_pair() {
    # 1_2
    pair_str=$1
    case_mos=$(echo $pair_str | cut -d "_" -f 1)
    ctrl_mos=$(echo $pair_str | cut -d "_" -f 2)
    read -a case_arr < <(create_single $case_mos | tr ';' ' ')
    read -a ctrl_arr < <(create_single $ctrl_mos | tr ';' ' ')
    ret=""
    for case_e in "${case_arr[@]}"; do
        for ctrl_e in "${ctrl_arr[@]}"; do
            ret+="${case_e}_${ctrl_e};"
        done
    done
    echo $ret >&2
    echo $ret
}

create_pair() {
    # 1_2 3_2 2_3
    pair_str=$1
    read -a pair_arr <<< $pair_str
    ret=""
    for pair in "${pair_arr[@]}"; do
        ret+=$(create_one_pair $pair)
    done
    echo $ret >&2
    echo $ret
}

create_range_pair() {
    # 3-1_9 1-1_9
    case_str=$1
    ctrl_str=$2
    read -a case_arr < <(create_one_single ${case_str} | tr ';' ' ')
    read -a ctrl_arr < <(create_one_single ${ctrl_str} | tr ';' ' ')
    ret=""
    for case_e in "${case_arr[@]}"; do
        for ctrl_e in "${ctrl_arr[@]}"; do
            ret+="${case_e}_${ctrl_e};"
        done
    done
    echo $ret >&2
    echo $ret
}

grab_out_dir() {
    base_dir=$1
    ds=$2
    if [[ $ds != "ND" ]]; then
        echo "${base_dir}/2.DS"
    else
        echo "${base_dir}/0.DF"
    fi
}

grab_analysis_dir() {
    base_dir=$1
    ds=$2
    if [[ $ds != "ND" ]]; then
        echo "${base_dir}/3.DS_A"
    else
        echo "${base_dir}/1.DF_A"
    fi
}

create_other() {
    in_dir=$1
    in_file_arr=($(ls ${in_dir}/*))
    ret=""
    for in_file in "${in_file_arr[@]}"; do
        ret+="${in_file};"
    done
    echo $ret >&2
    echo $ret
}
