# Import Modules
import sys
import os
sys.path.append(sys.path[0] + "/../common")
import bioparse as bp

# Input Parameters
REF_FILE = sys.argv[1]
IN_FILE = sys.argv[2]
OUT_FILE = sys.argv[3]


def __check_content(CHECK_FILE):
    check = False
    with open(CHECK_FILE, "r") as chck_f:
        for chck_l in chck_f:
            if chck_l[0] != "#":
                check = True
                break
    return check


def __both(out_f):
    ref_fp = bp.FileParser(REF_FILE)
    in_fp = bp.FileParser(IN_FILE)
    bp.write_vcf_hd(out_f)
    for ref_items in ref_fp.iter():
        if in_fp.term:
            out_f.write(f"{ref_fp.get_line().strip()}\n")
        else:
            cmp_res = bp.cmp_lines(ref_items, in_fp.get_line())
            while cmp_res > 0:
                out_f.write(f"{in_fp.get_line().strip()}\n")
                in_fp.next()
                if in_fp.term:
                    break
                cmp_res = bp.cmp_lines(ref_items, in_fp.get_line())
            if in_fp.term:
                out_f.write(f"{ref_fp.get_line().strip()}\n")
            elif cmp_res == 0:
                print(f"WARNING: Overlap for {ref_fp.get_line().strip()}")
            else:
                out_f.write(f"{ref_fp.get_line().strip()}\n")
    ref_fp.close()
    for _ in in_fp.iter():
        out_f.write(f"{in_fp.get_line().strip()}\n")
    in_fp.close()


if __name__ == '__main__':
    print(f"# Reference: {REF_FILE}")
    print(f"# Input: {IN_FILE}")
    print(f"# Output: {OUT_FILE}")
    can_run = False
    if os.path.isfile(REF_FILE) and __check_content(REF_FILE):
        if os.path.isfile(IN_FILE) and __check_content(IN_FILE):
            can_run = True
            out_f = open(OUT_FILE, "w")
            __both(out_f)
            out_f.close()
    if not can_run:
        print("ERROR: Invalid input files")
        sys.exit(1)
    print("# All operations complete")
