"""Parse by chromosome pileup file to extract variant/non-variant loci information."""


# Import Modules
import sys
sys.path.append(sys.path[0] + "/../common")
import os
import glob
import bioparse as bp
import pup


# Input Parameters
REF_DIR = sys.argv[1]
IN_FILE = sys.argv[2]
OUT_DIR = sys.argv[3]
NO_REFHOM = sys.argv[4]


# Global Variables
IN_NAME = IN_FILE.split("/")[-1].split(".")[:-1]
if "D" in IN_NAME[0]:
    IN_SAMP = ".".join(IN_NAME[:-1])
    IN_MOS = IN_NAME[-2].split("-")[0]
else:
    IN_SAMP = IN_NAME[0]
    IN_MOS = IN_SAMP.split("-")[0]
IN_CHR = IN_NAME[-1]
OUT_FILE = f"{OUT_DIR}/by_chr_{IN_SAMP}/{IN_SAMP}.{IN_CHR}.parse.txt"
ARG_TRUE_SET = {"true", "True", "TRUE"}
HQ_BQ_THRESHOLD = 30
no_refhom = NO_REFHOM in ARG_TRUE_SET


def __grab_ref():
    ref_list = list()
    ref_file_list = glob.glob(f"{REF_DIR}/*.p")
    for ref_file in ref_file_list:
        if ("pc" in ref_file and IN_MOS in ref_file) or "pc" not in ref_file:
            if no_refhom and "refhom" in ref_file:
                continue
            print(ref_file)
            ref_list.append(bp.load_pickle(ref_file, False))
    return tuple(ref_list)


def __grab_bq(cur_alt, alt_bq_dict):
    hq_bq_cnt: int = 0
    avg_bq: float = 0.0
    if cur_alt in alt_bq_dict:
        bq_sum: int = 0
        for bq in alt_bq_dict[cur_alt]:
            bq_sum += bq
            if bq >= HQ_BQ_THRESHOLD:
                hq_bq_cnt += 1
        avg_bq = round((bq_sum / len(alt_bq_dict[cur_alt])), 8)
    return hq_bq_cnt, avg_bq


if __name__ == '__main__':
    if not os.path.isfile(IN_FILE):
        print("ERROR: Invliad input file")
        sys.exit()
    if not os.path.isdir(f"{OUT_DIR}/by_chr_{IN_SAMP}"):
        os.mkdir(f"{OUT_DIR}/by_chr_{IN_SAMP}")
    print(f"** Input: {IN_FILE}")
    print(f"** Output: {OUT_FILE}")
    in_fp = bp.FileParser(IN_FILE)
    print("** Loading reference...")
    ref_tup = __grab_ref()
    print(in_fp.get_line().strip())
    out_f = open(OUT_FILE, "w")
    for in_items in in_fp.iter():
        in_fp.print_prog()
        in_chr = in_items[0]
        if in_chr != IN_CHR:
            break
        in_pos = int(in_items[1])
        ref_info = None
        for ref_tree_dict in ref_tup:
            if in_chr in ref_tree_dict:
                q_res = ref_tree_dict[in_chr].search(in_pos)
                if q_res:
                    ref_info = q_res.split(";")
                    break
        if not ref_info:
            continue
        dep: int = int(in_items[3])
        pup_info = pup.parse_base_and_bq(in_items[-2], in_items[-1])
        tot_alt_cnt: int = pup_info["tot_alt_cnt"]
        uniq_alt_cnt: int = len(pup_info["alt_cnt_dict"])
        if "refhom" in ref_info[0]:
            alt = pup_info["maj_alt"]
            alt_cnt = pup_info["maj_alt_cnt"]
        elif "snv" in ref_info[0]:
            alt = ref_info[2]
            if alt in pup_info["alt_cnt_dict"]:
                alt_cnt = pup_info["alt_cnt_dict"][alt]
            else:
                alt_cnt = 0
        else:
            if ref_info[1] > ref_info[2]:
                alt = f"-{ref_info[1][1:]}"
            else:
                alt = f"+{ref_info[2][1:]}"
            if alt in pup_info["alt_cnt_dict"]:
                alt_cnt = pup_info["alt_cnt_dict"][alt]
            else:
                alt_cnt = 0
        alt_hq_bq_cnt, alt_avg_bq = __grab_bq(alt, pup_info["alt_bq_dict"])
        maj_alt_hq_bq_cnt, maj_alt_avg_bq = __grab_bq(pup_info["maj_alt"], pup_info["alt_bq_dict"])
        if dep == 0:
            print("WARNING: Depth is zero")
            print(in_fp.get_line().strip())
            continue
        af: float = round((alt_cnt / dep), 8)
        maj_alt_af: float = round((pup_info["maj_alt_cnt"] / dep), 8)
        out_list = in_items[:2]
        out_list.append(str(dep))
        out_list.append(str(dep - tot_alt_cnt))
        out_list.append(str(alt))
        out_list.append(str(alt_cnt))
        out_list.append(str(af))
        out_list.append(str(alt_hq_bq_cnt))
        out_list.append(str(alt_avg_bq))
        out_list.append(str(pup_info["maj_alt"]))
        out_list.append(str(pup_info["maj_alt_cnt"]))
        out_list.append(str(maj_alt_af))
        out_list.append(str(maj_alt_hq_bq_cnt))
        out_list.append(str(maj_alt_avg_bq))
        out_list.append(str(uniq_alt_cnt))
        alt_dict_str = ""
        for cur_alt, cur_alt_cnt in pup_info["alt_cnt_dict"].items():
            cur_hq_bq_cnt, cur_avg_bq = __grab_bq(cur_alt, pup_info["alt_bq_dict"])
            alt_dict_str = f"{alt_dict_str};{cur_alt}_{cur_alt_cnt}_{cur_hq_bq_cnt}_{cur_avg_bq}"
        out_list.append(alt_dict_str[1:] if len(alt_dict_str) > 0 else "*")
        out_list.append(ref_info[0])
        out_str = "\t".join(out_list)
        out_f.write(f"{out_str}\n")
    out_f.close()
    print("** All operations complete")
