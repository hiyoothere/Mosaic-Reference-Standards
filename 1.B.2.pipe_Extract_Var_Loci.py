"""Extract variant loci with passing filters and diploid biallelic genotype from VCF file."""


# Import Modules
import sys
sys.path.append(sys.path[0] + "/../common")
import bioparse as bp


# Input Parameters
IN_FILE = sys.argv[1]
OUT_DIR = sys.argv[2]

# Global Variables
IN_FILENAME = IN_FILE.split("/")[-1]
IN_FILENAME = IN_FILENAME[:-7] if "gz" == IN_FILENAME[-2:] else IN_FILENAME[:-4]
PASS_FILE = f"{OUT_DIR}/{IN_FILENAME}.pass.vcf"
FAIL_FILE = f"{OUT_DIR}/{IN_FILENAME}.fail.vcf"
VALID_GT_SET = {"0/1", "1/0", "1/1"}


def __convert_gt(gt_str: str):
    if "." in gt_str:
        return None
    if "|" in gt_str:
        gt_str = gt_str.replace("|", "/")
    if "/" not in gt_str:
        gt_str = f"{gt_str}/{gt_str}"
    return gt_str


def __validate_line(in_items) -> bool:
    if in_items[6] != "PASS":
        return False
    cur_gt = __convert_gt(in_items[-1].split(":")[0])
    if not cur_gt or cur_gt not in VALID_GT_SET:
        return False
    alt_base_list = in_items[4].split(",")
    if "<*>" in alt_base_list:
        alt_base_list.remove("<*>")
    if len(alt_base_list) > 1:
        return False
    return True


if __name__ == '__main__':
    print(f"# Input: {IN_FILE}")
    in_fp = bp.FileParser(IN_FILE)
    pass_f = open(PASS_FILE, "w")
    fail_f = open(FAIL_FILE, "w")
    bp.write_vcf_hd(pass_f)
    bp.write_vcf_hd(fail_f)
    prev_chrpos: str = "chr1.0"
    for in_items in in_fp.iter():
        in_fp.print_prog()
        in_l = in_fp.get_line().rstrip()
        cur_chrpos: str = f"{in_items[0]}.{in_items[1]}")
        if cur_chrpos == prev_chrpos:
            print(f"WARNING: Multiallelic site: {in_l}")
            fail_f.write(f"{in_l}\n")
        elif __validate_line(in_items):
            pass_f.write(f"{in_l}\n")
        else:
            print(f"WARNING: Filter fail or non-diploid site: {in_l}")
            fail_f.write(f"{in_l}\n")
        prev_chrpos = cur_chrpos
    in_fp.close()
    pass_f.close()
    fail_f.close()
    print("# All operations complete")
