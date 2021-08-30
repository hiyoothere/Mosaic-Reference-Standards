"""Identify per chromosome mutually exclusive variant or unanimous non-variant loci from pickled files."""


# Import Modules
import sys
import glob
import time
sys.path.append(sys.path[0] + "/../common")
import bioparse as bp


# Input Parameters
IN_DIR = sys.argv[1]
IN_CHR = sys.argv[2]
OUT_DIR = sys.argv[3]


# Global Variables
IN_CHR = IN_CHR if "chr" in IN_CHR else f"chr{IN_CHR}"
OUT_ID = f"cmb_cl.{IN_CHR}"
OUT_PC_FILE = f"{OUT_DIR}/pc.{OUT_ID}.vcf"
OUT_NC_WT_FILE = f"{OUT_DIR}/nc_wt.{OUT_ID}.vcf"
PC_HOLDER_STR = "\t".join(tuple("..."))
NC_WT_HOLDER_STR = "\t".join(tuple("....."))
GT_TO_PL_DICT = {"0/1": "het", "1/0": "het", "1/1": "hom"}


def __get_info_for_chr(in_dir, in_chr):
    p_list = glob.glob(f"{in_dir}/*p")
    p_dict = {"var_fail": dict(), "var_pass": dict(), "nonvar": dict()}
    for p_file in p_list:
        p_filename = p_file.split("/")[-1]
        p_type = None
        if "genome.pass" in p_filename:
            p_type = "nonvar"
        elif "variants.pass" in p_filename:
            p_type = "var_pass"
        elif "variants.fail" in p_filename:
            p_type = "var_fail"
        if p_type:
            print(p_filename)
            p_id = p_filename.split(".")[0]
            dict_tree = bp.load_pickle(p_file, False)
            if in_chr in dict_tree:
                p_dict[p_type][p_id] = dict_tree[in_chr]
    return p_dict


def __grab_info_per_locus(info_dict, cur_pos):
    ret_dict = dict()
    for info_type in info_dict:
        for cl_id in info_dict[info_type]:
            q_res = info_dict[info_type][cl_id].search(cur_pos)
            if not q_res:
                continue
            if info_type not in ret_dict:
                ret_dict[info_type] = dict()
            ret_dict[info_type][cl_id] = q_res
    return ret_dict


def __process_per_locus_info(locus_info):
    if "var_fail" in locus_info or "nonvar" not in locus_info:
        return None
    if len(locus_info["nonvar"]) == cl_cnt:
        return ("nonvar", "nonvar")
    elif len(locus_info["nonvar"]) == cl_cnt - 1:
        if "var_pass" not in locus_info:
            return None
        if len(locus_info["var_pass"]) != 1:
            return None
        var_cl_id = next(iter(locus_info["var_pass"]))
        if var_cl_id in locus_info["nonvar"]:
            return None
        return (var_cl_id, locus_info["var_pass"][var_cl_id])
    return None


def __validate_var_info(var_cl_id, var_info, info_dict, cl_cnt, cur_pos):
    if var_cl_id == "S0" or var_cl_id == "s0":
        return None
    cur_ref, cur_alt, cur_gt = var_info.split("-")
    cur_pl = GT_TO_PL_DICT[cur_gt]
    if len(cur_ref) == 1 and len(cur_alt) == 1:
        return f"{var_cl_id}_{cur_pl}_snv\t{cur_ref}\t{cur_alt}"
    elif len(cur_ref) > len(cur_alt):
        for shift in range(len(cur_alt), len(cur_ref)):
            chk_pos = cur_pos + shift
            chk_locus_info = __grab_info_per_locus(info_dict, chk_pos)
            # Check deletion loci for any variant calls, passed or failed
            if "var_pass" in chk_locus_info or "var_fail" in chk_locus_info:
                return None
        return f"{var_cl_id}_{cur_pl}_ind\t{cur_ref}\t{cur_alt}"
    elif len(cur_ref) < len(cur_alt):
        return f"{var_cl_id}_{cur_pl}_ind\t{cur_ref}\t{cur_alt}"
    return None


def __process_info(info_dict, inv_max, cl_cnt, cur_chr, pc_file, nc1_file):
    init_time = time.time()
    pc_f = open(pc_file, "w")
    nc1_f = open(nc1_file, "w")
    bp.write_vcf_hd(pc_f)
    bp.write_vcf_hd(nc1_f)
    for cur_pos in range(1, inv_max + 1):
        bp.print_prog(cur_pos - 1, init_time)
        check_locus = False
        out_l_head = f"{cur_chr}\t{cur_pos}"
        locus_info = __grab_info_per_locus(info_dict, cur_pos)
        proc_res = __process_per_locus_info(locus_info)
        if not proc_res:
            check_locus = True
        elif "nonvar" in proc_res:
            nc1_f.write(f"{out_l_head}\tnonvar\t{NC_WT_HOLDER_STR}\n")
        else:
            var_info = __validate_var_info(proc_res[0], proc_res[1], info_dict, cl_cnt, cur_pos)
            if var_info:
                pc_f.write(f"{out_l_head}\t{var_info}\t{PC_HOLDER_STR}\n")
            else:
                check_locus = True
        if check_locus:
            locus_cl_sum = 0
            for info_type in locus_info:
                locus_cl_sum += len(locus_info[info_type])
            if locus_cl_sum > cl_cnt:
                print(f"ERROR: Cell line count exceeds given > {locus_info}")
    pc_f.close()
    nc1_f.close()


if __name__ == '__main__':
    print(f"# Input directory: {IN_DIR}")
    info_dict = __get_info_for_chr(IN_DIR, IN_CHR)
    print("# Processing combined genotype for all cell lines...")
    inv_max = 0
    cl_cnt = 0
    for info_type in info_dict:
        cur_cl_cnt = len(info_dict[info_type])
        if cur_cl_cnt > cl_cnt:
            cl_cnt = cur_cl_cnt
        for cl_id in info_dict[info_type]:
            cur_max = info_dict[info_type][cl_id].root.max
            if cur_max > inv_max:
                inv_max = cur_max
    print(f"Total cell line count: {cl_cnt}, Chromosome max locus: {inv_max}")
    __process_info(info_dict, inv_max, cl_cnt, IN_CHR, OUT_PC_FILE, OUT_NC_WT_FILE)
    print("# All operations complete")
