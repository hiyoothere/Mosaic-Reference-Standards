"""Identify and extract loci with low average depth of coverage across samples."""


# Import Modules
import sys
import statistics
sys.path.append(sys.path[0] + "/../common")
import bioparse as bp


# Input Parameters
IN_FILE = sys.argv[1]
THRESHOLD = sys.argv[2]
OUT_DIR = sys.argv[3]


# Global Variables
OUT_FILE = f"{OUT_DIR}/{IN_FILE.split('/')[-1][:-4]}.low_avg_dp.vcf"
threshold = int(THRESHOLD)
HOLDER_STR = "\t".join(tuple("......"))


if __name__ == '__main__':
    print(f"# Input: {IN_FILE}")
    print(f"# Depth threshold: {threshold}")
    in_fp = bp.FileParser(IN_FILE)
    out_f = open(OUT_FILE, "w")
    bp.write_vcf_hd(out_f)
    tag_count_dict = dict()
    for in_items in in_fp.iter():
        in_fp.print_prog()
        raw_list = in_items[-1].split(";")
        target_list = list()
        for raw_item in raw_list:
            target_list.append(int(raw_item.split("_")[1]))
        cur_tag = in_items[-2]
        cur_avg_dp = statistics.mean(target_list)
        if cur_avg_dp < threshold:
            if in_items[-2] not in tag_count_dict:
                tag_count_dict[in_items[-2]] = 0
            tag_count_dict[in_items[-2]] += 1
            print(f"Avg: {cur_avg_dp:<6.3f} | Raw: {in_fp.get_line().strip()}")
            out_f.write(f"{in_items[0]}\t{in_items[1]}\t{HOLDER_STR}\n")
    in_fp.close()
    out_f.close()
    print("# Statistics")
    tot_count = 0
    for tag in sorted(tag_count_dict.keys()):
        print(f"Tag: {tag} | Count: {tag_count_dict[tag]}")
        tot_count += tag_count_dict[tag]
    print(f"Total: {tot_count}")
