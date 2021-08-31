"""Merge pileup parsed per loci information for all chromosomes by sample."""


# Import Modules
import sys
sys.path.append(sys.path[0] + "/../common")
import glob
import bioparse as bp


# Input Parameters
IN_DIR = sys.argv[1]
IN_ID = sys.argv[2]
IN_DS = sys.argv[3]


# Global Variables
IN_ID = IN_ID if "N" in IN_DS else f"{IN_DS}.{IN_ID}"


if __name__ == '__main__':
    in_file_list = glob.glob(f"{IN_DIR}/by_chr_{IN_ID}/*{IN_ID}*.txt")
    print(f"** Input: {len(in_file_list)}")
    for in_file in in_file_list:
        print(in_file)
    print("** Merging for all chromosomes")
    out_file = f"{IN_DIR}/{IN_ID}.parse.txt"
    out_f = open(out_file, "w")
    out_f.write("#CHR\tPOS\tDEP\tREF_CNT\tALT\tALT_CNT\tALT_AF\tALT_HQ_BQ_CNT\tALT_AVG_BQ\tMAJ_ALT\tMAJ_ALT_CNT\tMAJ_ALT_AF\tMAJ_HQ_BQ_CNT\tMAJ_ALT_AVG_BQ\tUNIQ_ALT_CNT\tALT_INFO\tTAG\n")
    for chr_num in bp.iter_chrom():
        cur_file = None
        for by_chr_file in in_file_list:
            if f".chr{chr_num}." in by_chr_file:
                cur_file = by_chr_file
                break
        if not cur_file:
            print(f"Skipping chr{chr_num}")
        else:
            cur_fp = bp.FileParser(cur_file)
            for cur_item in cur_fp.iter():
                out_f.write(f"{cur_fp.get_line().strip()}\n")
    out_f.close()
    print("** All operations complete")
