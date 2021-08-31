"""Identify and extract loci with low reproducibility across samples."""


# Import Modules
import sys
sys.path.append(sys.path[0] + "/../common")
import bioparse as bp


# Input Parameters
CTRL_FILE = sys.argv[1]
IN_FILE = sys.argv[2]
THRESHOLD = sys.argv[3]
OUT_DIR = sys.argv[4]
IGNORE_NULL = sys.argv[5] in {"True", "true", "T", "t"}


# Global Variables
OUT_FILE = f"{OUT_DIR}/{IN_FILE.split('/')[-1][:-4]}.low_rep.vcf"
HOLDER_STR = "\t".join(tuple("......"))
threshold = float(THRESHOLD)
PC_SIZE_DICT = {"V1": 39, "V2": 9, "V3": 30, "V4": 12, "V5": 18}


def __test_rep(in_items, ctrl_info):
    check_idx = 1
    tot_cnt = 39
    ctrl_type = ctrl_info.split("_")[0].upper()
    if ctrl_type[0] == "V":
        check_idx = 3
        tot_cnt = PC_SIZE_DICT[ctrl_type]
    val_cnt = 0
    raw_list = in_items[-1].split(";")
    for raw_item in raw_list:
        if int(raw_item.split("_")[check_idx]) > 0:
            val_cnt += 1
    cur_thres = tot_cnt * threshold
    if val_cnt < cur_thres:
        return (ctrl_type, cur_thres, val_cnt)
    return None


def __print_and_write(out_f, chrpos, tag_count_dict, rep_test=None):
    if rep_test:
        ctrl_info = rep_test[0]
        if ctrl_info not in tag_count_dict:
            tag_count_dict[ctrl_info] = 0
        tag_count_dict[ctrl_info] += 1
        print(f"{chrpos} | CurThres: {rep_test[1]:<6.3f} | Count: {rep_test[2]}")
    else:
        tag_count_dict["na"] += 1
        print(f"{chrpos} | CurThres: ??? | Count: 0")
    out_f.write(f"{chrpos[0]}\t{chrpos[1]}\t{HOLDER_STR}\n")


if __name__ == '__main__':
    print(f"# Reference: {CTRL_FILE}")
    print(f"# Input: {IN_FILE}")
    print(f"# Reproducibility threshold: {threshold}")
    ref_fp = bp.FileParser(CTRL_FILE)
    in_fp = bp.FileParser(IN_FILE)
    out_f = open(OUT_FILE, "w")
    bp.write_vcf_hd(out_f)
    tag_count_dict = {"na": 0}
    for ref_items in ref_fp.iter():
        if in_fp.term:
            if not IGNORE_NULL:
                __print_and_write(out_f, ref_items, tag_count_dict)
            continue
        cmp_res = bp.cmp_lines(ref_items, in_fp.get_line())
        while cmp_res > 0:
            in_fp.next()
            if in_fp.term:
                break
            cmp_res = bp.cmp_lines(ref_items, in_fp.get_line())
        if cmp_res == 0:
            rep_test_res = __test_rep(in_fp.get_items(), ref_items[2])
            if rep_test_res:
                __print_and_write(out_f, ref_items, tag_count_dict, rep_test_res)
        else:
            if not IGNORE_NULL:
                __print_and_write(out_f, ref_items, tag_count_dict)
    in_fp.close()
    out_f.close()
    print("# Statistics")
    tot_count = 0
    for tag in sorted(tag_count_dict.keys()):
        print(f"Tag: {tag} | Count: {tag_count_dict[tag]}")
        tot_count += tag_count_dict[tag]
    print(f"Total: {tot_count}")
