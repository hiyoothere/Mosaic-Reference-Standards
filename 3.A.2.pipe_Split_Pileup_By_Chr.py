"""Split samtools pileup files by chromosomes for parallel processing."""


# Import Modules
import sys
sys.path.append(sys.path[0] + "/../common")
import time
import os


# Input Parameters
IN_PUP = sys.argv[1]


# Global Variables
IN_PATH = IN_PUP.split("/")
IN_DIR = "/".join(IN_PATH[:-1])
IN_ID = IN_PATH[-1][:-4]
OUT_DIR = f"{IN_DIR}/{IN_ID}"
TERM_CHR_SET = {"chrM", "chrm"}


if __name__ == '__main__':
    if not os.path.isdir(OUT_DIR):
        os.mkdir(OUT_DIR)
    time.sleep(0.1)
    out_f = None
    with open(IN_PUP, "r") as in_f:
        prev_chr = "chr0"
        init_time = time.time()
        for in_l in in_f:
            cur_chr = in_l.split()[0]
            if prev_chr != cur_chr:
                if cur_chr in TERM_CHR_SET:
                    break
                if out_f:
                    out_f.close()
                out_file = f"{OUT_DIR}/{IN_ID}.{cur_chr}.pup"
                run_time = (time.time()-init_time)/60
                print(f"** Writing output: {out_file} | Run time: {run_time:.2f}m")
                out_f = open(out_file, "w")
                prev_chr = cur_chr
            out_f.write(f"{in_l.strip()}\n")
        if out_f:
            out_f.close()
    run_time = (time.time()-init_time)/60
    print(f"** Total run time: {run_time:.2f}m")
    print("** All operations done")
