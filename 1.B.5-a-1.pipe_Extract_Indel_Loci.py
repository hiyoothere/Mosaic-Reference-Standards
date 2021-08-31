# Import Modules
import sys
sys.path.append(sys.path[0] + "/../common")
import bioparse as bp


# Input Parameters
IN_FILE = sys.argv[1]

# Global Variables
IN_FILENAME = IN_FILE[:-7] if "gz" == IN_FILE[-2:] else IN_FILE[:-4]
OUT_FILE = f"{IN_FILENAME}.ind.vcf"


if __name__ == '__main__':
    print(f"# Input: {IN_FILE}")
    in_fp = bp.FileParser(IN_FILE)
    out_f = open(OUT_FILE, "w")
    bp.write_vcf_hd(out_f)
    for in_items in in_fp.iter():
        if "ind" in in_items[2]:
            out_f.write(in_fp.get_line().strip() + "\n")
    in_fp.close()
    out_f.close()
    print("# All operations complete")
