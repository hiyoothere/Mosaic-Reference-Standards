"""Expand GVCF block loci into single position loci VCF split based on filter pass status."""


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
VALID_FIL_SET = {"PASS", "."}
VALID_GT_SET = {"0/0", "0|0", "0"}
HOLDER_STR = "\t".join(tuple("......"))


def __validate_line(in_items) -> bool:
    if in_items[6] not in VALID_FIL_SET:
        return False
    cur_gt = in_items[9].split(":")[0]
    if cur_gt not in VALID_GT_SET:
        return False
    return True


def __grab_block_end(in_str: str, pos: int):
    return int(in_str.split(";")[0].split("=")[1]) if "END" in in_str else pos


if __name__ == '__main__':
    print(f"# Input: {IN_FILE}")
    in_fp = bp.FileParser(IN_FILE)
    pass_f = open(PASS_FILE, "w")
    fail_f = open(FAIL_FILE, "w")
    bp.write_vcf_hd(pass_f)
    bp.write_vcf_hd(fail_f)
    for in_items in in_fp.iter():
        in_fp.print_prog()
        start_pos = int(in_items[1])
        end_pos = __grab_block_end(in_items[7], start_pos)
        if __validate_line(in_items):
            for cur_pos in range(start_pos, end_pos + 1):
                pass_f.write(f"{in_items[0]}\t{cur_pos}\t{HOLDER_STR}\n")
        else:
            for cur_pos in range(start_pos, end_pos + 1):
                fail_f.write(f"{in_items[0]}\t{cur_pos}\t{HOLDER_STR}\n")
    in_fp.close()
    pass_f.close()
    fail_f.close()
    print("# All operations complete")
