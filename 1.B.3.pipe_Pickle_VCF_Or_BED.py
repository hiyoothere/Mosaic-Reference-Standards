"""Convert VCF or BED file into by chromosome dictionary + red-black tree data structure pickle file."""


# Import Modules
import sys
sys.path.append(sys.path[0] + "/../common")
import tree
import bioparse as bp


# Input Parameters
IN_FILE = sys.argv[1]
OUT_DIR = sys.argv[2]


# Global Variables
HET_GT_SET = {"1/0", "0/1", "1|0", "0|1"}
HOM_GT_SET = {"1/1", "1|1", "1"}
OUT_FILE = f"{OUT_DIR}/{IN_FILE.split('/')[-1]}.p"


def construct_from_bed(in_file: str):
    dict_tree = dict()
    in_fp = bp.FileParser(in_file)
    for in_items in in_fp.iter():
        in_fp.print_prog()
        cur_chr: str = in_items[0]
        cur_start: int = int(in_items[1]) + 1
        cur_end: int = int(in_items[2])
        if cur_chr not in dict_tree:
            dict_tree[cur_chr] = tree.RBITree()
        q_res = dict_tree[cur_chr].search(cur_start)
        if q_res:
            print(f"WARNING: Overlap found: {cur_chr}:{cur_start} - {q_res}")
        dict_tree[cur_chr].insert(cur_start, "0/0", cur_end)
    in_fp.close()
    return dict_tree


def process_vcf_line(in_items):
    cur_ref: str = in_items[3]
    cur_alt: str = in_items[4]
    if "," in cur_alt:
        return ("non-dip",)
    cur_gt: str = in_items[-1].split(":")[0]
    if cur_gt in HET_GT_SET:
        return (cur_ref, cur_alt, "0/1")
    elif cur_gt in HOM_GT_SET:
        return (cur_ref, cur_alt, "1/1")
    else:
        return ("inval_gt",)


def construct_from_vcf(in_file: str):
    dict_tree = dict()
    in_fp = bp.FileParser(in_file)
    for in_items in in_fp.iter():
        in_fp.print_prog()
        cur_chr: str = in_items[0]
        cur_pos: int = int(in_items[1])
        if cur_chr not in dict_tree:
            dict_tree[cur_chr] = tree.RBITree()
        q_res = dict_tree[cur_chr].search(cur_pos)
        if q_res:
            print(f"WARNING: Overlap found: {cur_chr}:{cur_pos} - {q_res}")
        proc_res = process_vcf_line(in_items)
        cur_info: str = "-".join(proc_res)
        dict_tree[cur_chr].insert(cur_pos, cur_info)
    in_fp.close()
    return dict_tree


if __name__ == '__main__':
    print(f"# Input: {IN_FILE}")
    dict_tree = None
    if IN_FILE[-3:] == "vcf":
        dict_tree = construct_from_vcf(IN_FILE)
    elif IN_FILE[-3:] == "bed":
        dict_tree = construct_from_bed(IN_FILE)
    else:
        print(f"ERROR: Invalid file format > {IN_FILE}")
        sys.exit(1)
    dict_tree_size = 0
    for cur_chr in dict_tree:
        dict_tree_size += dict_tree[cur_chr].size()
    print(f"# Output: {OUT_FILE} of dict-tree size {dict_tree_size}...")
    bp.dump_pickle(dict_tree, OUT_FILE, False)
    print("# All operations complete")
