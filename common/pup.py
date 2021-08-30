"""Collection of functions to aid processing pileup files."""


# Import Modules
from typing import Dict, List
import os.path


# Key Variables
REF_PAT_SET = set(".,")
ALT_PAT_SET = set("atgcATGC")
IND_PAT_SET = set("+-")
BQ_SKIP_PAT_SET = set("*<>")
BQ_MIN_CNST: int = ord("!")
MIN_TRI_DEP: float = 0.01


def process_base(base_str: str):
    ref_cnt: int = 0
    alt_cnt_dic: Dict[str, int] = {"A": 0, "T": 0, "G": 0, "C": 0}
    i: int = 0
    while i < len(base_str):
        cur_base: str = base_str[i]
        if cur_base in REF_PAT_SET:
            ref_cnt += 1
        elif cur_base in ALT_PAT_SET:
            alt_cnt_dic[cur_base.upper()] += 1
        elif cur_base == "^":
            i += 1
        elif cur_base in IND_PAT_SET:
            ind_len_str: str = base_str[i + 1]
            skip: int = i + 2
            while base_str[skip].isdigit():
                ind_len_str = "{}{}".format(ind_len_str, base_str[skip])
                skip += 1
            i += int(ind_len_str) + len(ind_len_str)
        i += 1
    maj_alt_cnt: int = 0
    for nt, nt_cnt in alt_cnt_dic.items():
        if nt_cnt > maj_alt_cnt:
            maj_alt_cnt = nt_cnt
    maj_alt_nts: List[str] = list()
    for nt, nt_cnt in alt_cnt_dic.items():
        if nt_cnt == maj_alt_cnt:
            maj_alt_nts.append(nt)
    af: float = maj_alt_cnt / (maj_alt_cnt + ref_cnt) if maj_alt_cnt > 0 else 0.0
    af = round(af, 8)
    return {
        "af": af, "ref_cnt": ref_cnt, "alt_al": tuple(maj_alt_nts),
        "alt_cnt": maj_alt_cnt, "debug": alt_cnt_dic}


def process_indel_base(base_str: str):
    ref_cnt: int = 0
    alt_cnt_dic: Dict[str, int] = dict()
    i: int = 0
    while i < len(base_str):
        cur_base: str = base_str[i]
        if cur_base in REF_PAT_SET:
            ref_cnt += 1
        elif cur_base == "^":
            i += 1
        elif cur_base in IND_PAT_SET:
            ind_len_str: str = base_str[i + 1]
            ind_alt: str = ""
            skip: int = i + 2
            while base_str[skip].isdigit():
                ind_len_str = "{}{}".format(ind_len_str, base_str[skip])
                skip += 1
            ind_len = int(ind_len_str)
            for j in range(int(ind_len)):
                ind_alt = "{}{}".format(ind_alt, base_str[skip + j])
            i += ind_len + len(ind_len_str)
            ind_alt = ind_alt.upper()
            if ind_alt not in alt_cnt_dic:
                alt_cnt_dic[ind_alt] = 0
            alt_cnt_dic[ind_alt] += 1
        i += 1
    maj_alt_cnt: int = 0
    for ind, ind_cnt in alt_cnt_dic.items():
        if ind_cnt > maj_alt_cnt:
            maj_alt_cnt = ind_cnt
    maj_alt_inds: List[str] = list()
    for ind, ind_cnt in alt_cnt_dic.items():
        if ind_cnt == maj_alt_cnt:
            maj_alt_inds.append(ind)
    af: float = maj_alt_cnt / ref_cnt
    af = round(af, 8)
    pup_info = {
        "af": af, "ref_cnt": ref_cnt - maj_alt_cnt,
        "alt_al": tuple(maj_alt_inds), "alt_cnt": maj_alt_cnt,
        "debug": alt_cnt_dic}
    return pup_info


def process_base_qual(base_str, qual_str):
    base_idx = 0
    qual_idx = 0
    alt_cnt_dic = {"A": 0, "T": 0, "C": 0, "G": 0}
    alt_bq_dic = {"A": 0.0, "T": 0.0, "C": 0.0, "G": 0.0}
    ref_cnt = 0
    ref_bq = 0.0
    while base_idx < len(base_str) and qual_idx < len(qual_str):
        cur_nt = base_str[base_idx]
        cur_bq = ord(qual_str[qual_idx]) - BQ_MIN_CNST
        if cur_nt == "." or cur_nt == ",":
            ref_cnt += 1
            ref_bq += cur_bq
            qual_idx += 1
        elif cur_nt in ALT_PAT_SET:
            alt_cnt_dic[cur_nt.upper()] += 1
            alt_bq_dic[cur_nt.upper()] += cur_bq
            qual_idx += 1
        elif cur_nt == "^":
            base_idx += 1
        elif cur_nt == "*" or cur_nt == "<" or cur_nt == ">":
            qual_idx += 1
        elif cur_nt == "+" or cur_nt == "-":
            skip = base_idx + 1
            indel_len_str = ""
            while(base_str[skip].isdigit()):
                indel_len_str += base_str[skip]
                skip += 1
            indel_len = int(indel_len_str)
            base_idx = skip + indel_len - 1
        base_idx += 1
    maj_alt_nt = "."
    maj_alt_cnt = 0
    for nt in alt_cnt_dic:
        if alt_cnt_dic[nt] > maj_alt_cnt:
            maj_alt_nt = nt
            maj_alt_cnt = alt_cnt_dic[nt]
    af = maj_alt_cnt / (maj_alt_cnt + ref_cnt) if maj_alt_cnt > 0 else 0.0
    if maj_alt_cnt > 0 and af >= MIN_TRI_DEP:
        for nt in alt_cnt_dic:
            if nt != maj_alt_nt and alt_cnt_dic[nt] >= maj_alt_cnt:
                print("WARNING: Potential triallelic position: {}".format(alt_cnt_dic))
                break
    ref_avg_bq = ref_bq / ref_cnt if ref_cnt > 1 else 0.0
    maj_alt_avg_bq = 0.0 if maj_alt_nt == "." else alt_bq_dic[maj_alt_nt] / maj_alt_cnt
    return {"af": af, "ref_cnt": ref_cnt, "ref_bq": ref_avg_bq,
        "alt_nt": maj_alt_nt, "alt_cnt": maj_alt_cnt, "alt_bq": maj_alt_avg_bq}


def skip_to_chr(in_chr: str, in_file: str, in_fp):
    cur_chr: str = in_chr if "chr" in in_chr else "chr{}".format(in_chr)
    skip: int = 0
    idx_file = in_file if ".pui" in in_file else "{}.pui".format(in_file)
    if os.path.isfile(idx_file):
        with open(idx_file, "r") as pui_f:
            for pui_l in pui_f:
                pui_chr, pui_skip = pui_l.split()
                if pui_chr == cur_chr:
                    skip = int(pui_skip)
                    break
    if skip > 1:
        in_fp.long_next(skip - 1)
    return


def parse_snv(base_str: str):
    tot_alt_cnt: int = 0
    alt_cnt_dic: Dict[str, int] = dict()
    i: int = 0
    while i < len(base_str):
        cur_base: str = base_str[i]
        if cur_base in ALT_PAT_SET:
            tot_alt_cnt += 1
            if cur_base.upper() not in alt_cnt_dic:
                alt_cnt_dic[cur_base.upper()] = 0
            alt_cnt_dic[cur_base.upper()] += 1
        elif cur_base in IND_PAT_SET:
            tot_alt_cnt += 1
            ind_len_str: str = base_str[i + 1]
            skip: int = i + 2
            while base_str[skip].isdigit():
                ind_len_str = "{}{}".format(ind_len_str, base_str[skip])
                skip += 1
            i += int(ind_len_str) + len(ind_len_str)
        elif cur_base == "^":
            i += 1
        i += 1
    maj_alt: str = "*"
    maj_alt_cnt: int = 0
    for alt, alt_cnt in alt_cnt_dic.items():
        cnt_diff: int = alt_cnt - maj_alt_cnt
        if cnt_diff > 0:
            maj_alt = alt
            maj_alt_cnt = alt_cnt
        elif alt_cnt > 0 and cnt_diff == 0:
            maj_alt = "{},{}".format(maj_alt, alt)
    return {
        "tot_alt_cnt": tot_alt_cnt, "alt_dict": alt_cnt_dic,
        "maj_alt": maj_alt, "maj_alt_cnt": maj_alt_cnt}


def parse_ind(base_str: str):
    tot_alt_cnt: int = 0
    alt_cnt_dic: Dict[str, int] = dict()
    i: int = 0
    while i < len(base_str):
        cur_base: str = base_str[i]
        if cur_base in IND_PAT_SET:
            tot_alt_cnt += 1
            ind_len_str: str = base_str[i + 1]
            skip: int = i + 2
            while base_str[skip].isdigit():
                ind_len_str = "{}{}".format(ind_len_str, base_str[skip])
                skip += 1
            ind_len: int = int(ind_len_str)
            ind_alt: str = "{}{}".format(cur_base, base_str[skip:skip+ind_len].upper())
            i += ind_len + len(ind_len_str)
            if ind_alt not in alt_cnt_dic:
                alt_cnt_dic[ind_alt] = 0
            alt_cnt_dic[ind_alt] += 1
        elif cur_base in ALT_PAT_SET:
            tot_alt_cnt += 1
        elif cur_base == "^":
            i += 1
        i += 1
    maj_alt: str = "*"
    maj_alt_cnt: int = 0
    for alt, alt_cnt in alt_cnt_dic.items():
        cnt_diff: int = alt_cnt - maj_alt_cnt
        if cnt_diff > 0:
            maj_alt = alt
            maj_alt_cnt = alt_cnt
        elif alt_cnt > 0 and cnt_diff == 0:
            maj_alt = "{},{}".format(maj_alt, alt)
    return {
        "tot_alt_cnt": tot_alt_cnt, "alt_dict": alt_cnt_dic,
        "maj_alt": maj_alt, "maj_alt_cnt": maj_alt_cnt}


def __count_helper(cur_base, cnt_dict):
    if cur_base not in cnt_dict:
        cnt_dict[cur_base] = 0
    cnt_dict[cur_base] += 1


def parse_base_and_bq(base_str: str, bq_str: str):
    cnt_dict: Dict[str, int] = dict()
    bq_dict = dict()
    base_idx: int = 0
    bq_idx: int = 0
    while base_idx < len(base_str):
        cur_base: str = base_str[base_idx]
        if cur_base in REF_PAT_SET:
            bq_idx += 1
        elif cur_base in ALT_PAT_SET:
            cur_bq: int = ord(bq_str[bq_idx]) - BQ_MIN_CNST
            __count_helper(cur_base.upper(), cnt_dict)
            if cur_base not in bq_dict:
                bq_dict[cur_base] = list()
            bq_dict[cur_base].append(cur_bq)
            bq_idx += 1
        elif cur_base in IND_PAT_SET:
            ind_len_str: str = base_str[base_idx + 1]
            skip: int = base_idx + 2
            while base_str[skip].isdigit():
                ind_len_str = "{}{}".format(ind_len_str, base_str[skip])
                skip += 1
            ind_len: int = int(ind_len_str)
            ind_base: str = "{}{}".format(cur_base, base_str[skip:skip + ind_len])
            __count_helper(ind_base.upper(), cnt_dict)
            base_idx += len(ind_len_str) + ind_len
        elif cur_base == "^":
            base_idx += 1
        elif cur_base in BQ_SKIP_PAT_SET:
            bq_idx += 1
        base_idx += 1
    tot_alt_cnt: int = 0
    maj_alt: str = "*"
    maj_alt_cnt: int = 0
    for alt, alt_cnt in cnt_dict.items():
        tot_alt_cnt += alt_cnt
        cnt_diff: int = alt_cnt - maj_alt_cnt
        if cnt_diff > 0:
            maj_alt = alt
            maj_alt_cnt = alt_cnt
        elif alt_cnt > 0 and cnt_diff == 0:
            maj_alt = "{},{}".format(maj_alt, alt)
    return {
        "alt_cnt_dict": cnt_dict, "alt_bq_dict": bq_dict,
        "tot_alt_cnt": tot_alt_cnt, "maj_alt": maj_alt, "maj_alt_cnt": maj_alt_cnt}
