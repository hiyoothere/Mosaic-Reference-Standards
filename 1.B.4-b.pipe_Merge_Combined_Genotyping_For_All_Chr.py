"""Merge combined genotyping results for all chromosomes."""

# Import Modules
import sys
import glob
sys.path.append(sys.path[0] + "/../common")
import bioparse as bp


# Input Parameters
IN_DIR = sys.argv[1]
IN_ID = sys.argv[2]
OUT_FILE = sys.argv[3]


def __grab_all_input(in_dir, in_id):
    in_wc = in_id if "*" in in_id else f"{in_id}*"
    return glob.glob(f"{in_dir}/{in_wc}")


if __name__ == '__main__':
    in_file_list = __grab_all_input(IN_DIR, IN_ID)
    print(f"# Inputs: {len(in_file_list)}")
    for in_file in in_file_list:
        print(in_file)
    out_f = open(OUT_FILE, "w")
    bp.write_vcf_hd(out_f)
    for chr_num in bp.iter_chrom():
        cur_file = None
        for in_file in in_file_list:
            if f".chr{chr_num}." in in_file:
                cur_file = in_file
                break
        if not cur_file:
            continue
        cur_fp = bp.FileParser(cur_file)
        for cur_item in cur_fp.iter():
            out_f.write(f"{cur_fp.get_line().strip()}\n")
        cur_fp.close()
    out_f.close()
