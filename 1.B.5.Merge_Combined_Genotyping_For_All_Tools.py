"""Merge combined genotyping results for all germline genotyping tools."""


# Import Modules
import sys
sys.path.append(sys.path[0] + "/../common")
import bioparse as bp


# Input Parameters
SK_FILE = sys.argv[1]
DV_FILE = sys.argv[2]
OUT_DIR = sys.argv[3]
OUT_ID = sys.argv[4]


# Global Variables
IN_FILENAME = SK_FILE.split("/")[-1]
IN_FILENAME = IN_FILENAME[:-7] if "gz" == IN_FILENAME[-2:] else IN_FILENAME[:-4]
OUT_FILE = f"{OUT_DIR}/{OUT_ID}.{IN_FILENAME}.vcf"
HET_GT_SET = {"1/0", "0/1", "1|0", "0|1"}
HOM_GT_SET = {"1/1", "1|1", "1"}
HOLDER_STR = "\t".join(tuple("..."))


def __check_if_proc_vcf(in_items):
    return in_items[2] != "." and in_items[-1] == "."


def __grab_gt(format_str):
    raw_gt = format_str.split(":")[0]
    if raw_gt in HET_GT_SET:
        return "het"
    elif raw_gt in HOM_GT_SET:
        return "hom"
    else:
        return "non-dip"


def __check_bases(ref_items, in_items):
    return ref_items[3] == in_items[3] and ref_items[4] == in_items[4]


def __compare_vcf_info(ref_items, in_items):
    if not __check_bases(ref_items, in_items):
        return None
    elif not __check_if_proc_vcf(ref_items):
        ref_gt = __grab_gt(ref_items[-1])
        in_gt = __grab_gt(in_items[-1])
        if ref_gt == "non-dip" or ref_gt != in_gt:
            return None
        out_list = ref_items[:2]
        var_id = ref_gt
        var_id += "_ind" if len(ref_items[3]) != len(ref_items[4]) else "_snv"
        out_list.append(var_id)
        out_list += ref_items[3:5]
        return "{}\t{}\n".format("\t".join(out_list), HOLDER_STR)
    elif sk_items[2] == in_items[2]:
        return "{}\n".format("\t".join(ref_items))
    return None


if __name__ == '__main__':
    print(f"# SK Input: {SK_FILE}")
    print(f"# DV Input: {DV_FILE}")
    sk_fp = bp.FileParser(SK_FILE)
    dv_fp = bp.FileParser(DV_FILE)
    out_f = open(OUT_FILE, "w")
    bp.write_vcf_hd(out_f)
    for sk_items in sk_fp.iter():
        if dv_fp.term:
            break
        cmp_res = bp.cmp_lines(sk_items, dv_fp.get_items())
        while cmp_res > 0:
            dv_fp.next()
            if dv_fp.term:
                break
            cmp_res = bp.cmp_lines(sk_items, dv_fp.get_items())
        if dv_fp.term:
            break
        elif cmp_res == 0:
            vcf_cmp_res = __compare_vcf_info(sk_items, dv_fp.get_items())
            if vcf_cmp_res:
                out_f.write(f"{vcf_cmp_res}")
    sk_fp.close()
    dv_fp.close()
    out_f.close()
    print("# All operations complete")
