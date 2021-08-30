"""Merge by control tag pileup parsed information for all samples."""


# Import Modules
import sys
sys.path.append(sys.path[0] + "/../common")
import glob
import re
import bioparse as bp
import tree


# Input Parameters
IN_DIR = sys.argv[1]
IN_ID = sys.argv[2]
OUT_DIR = sys.argv[3]
OUT_ID = sys.argv[4]


# Global Variables
SAMP_PAT = re.compile("M[1-3]-[0-9]+")
ALT_SAMP_SET = {"S0", "V1", "V2", "V3", "V4", "V5"}
FLOAT_PREC = 6


def __grab_all_input(in_dir, in_id):
    all_in_path_list = glob.glob(f"{in_dir}/*{in_id}.txt")
    in_file_dict = dict()
    for in_path in all_in_path_list:
        file_name = in_path.split("/")[-1]
        samp_id = SAMP_PAT.findall(file_name)
        if len(samp_id) > 0:
            in_file_dict[in_path] = samp_id[0]
        else:
            for alt_samp_id in ALT_SAMP_SET:
                if alt_samp_id in file_name:
                    in_file_dict[in_path] = alt_samp_id
                    break
    return in_file_dict


def __construct_tree_dict(in_file, in_id, merge_tree_dict):
    in_fp = bp.FileParser(in_file)
    for in_items in in_fp.iter():
        in_fp.print_prog()
        cur_chr = in_items[0][3:]
        cur_pos = int(in_items[1])
        cur_tag = in_items[-1]
        if "V" in cur_tag:
            type_tag = "pc_ind" if "ind" in cur_tag else "pc_snv"
        elif "refhom" in cur_tag:
            type_tag = "nc_low"
        else:
            type_tag = "nc_ind" if "ind" in cur_tag else "nc_snv"
        af = round(float(in_items[5]), FLOAT_PREC)
        cur_info = f"{in_id}_{'_'.join(in_items[2:5])}_{af}_{'_'.join(in_items[6:8])}"
        if type_tag not in merge_tree_dict:
            merge_tree_dict[type_tag] = dict()
        if cur_chr not in merge_tree_dict[type_tag]:
            merge_tree_dict[type_tag][cur_chr] = tree.RBTree()
        prev_info = merge_tree_dict[type_tag][cur_chr].search(cur_pos)
        if prev_info:
            if prev_info.split(";")[0] != cur_tag:
                print("WARNING: Inconsistent tag for locus")
            cur_info = f"{prev_info};{cur_info}"
        else:
            cur_info = f"{cur_tag};{cur_info}"
        merge_tree_dict[type_tag][cur_chr].insert(cur_pos, cur_info)
    return


def __create_out_dict(type_tag_list, in_id, out_dir, out_id):
    out_dict = dict()
    for type_tag in type_tag_list:
        out_dict[type_tag] = f"{out_dir}/{out_id}.{type_tag}.txt"
    return out_dict


def __write_out_file(out_file, tag_merge_tree_dict):
    out_f = open(out_file, "w")
    out_f.write("#CHR\tPOS\tVAR\tMERGED\n")
    for chr_key in bp.iter_chrom():
        if chr_key not in tag_merge_tree_dict:
            print(f"WARNING: Skipping chr{chr_key}")
            continue
        for pos_val, info_str in tag_merge_tree_dict[chr_key].inorder():
            info_list = info_str.split(";")
            out_info = f"{info_list[0]}\t{';'.join(info_list[1:])}"
            expexted_count = 40
            if "nc" not in info_list[0]:
                if "V2" in info_list[0]:
                    expexted_count = 10
                elif "V3" in info_list[0]:
                    expexted_count = 31
                elif "V4" in info_list[0]:
                    expexted_count = 13
                elif "V5" in info_list[0]:
                    expexted_count = 19
            if len(info_list) != expexted_count:
                print(f"WARNING: Inconsistent for {chr_key}:{pos_val} > Expected: {expexted_count} / Parsed: {len(info_list)}")
                print(sorted(info_list))
            out_f.write(f"chr{chr_key}\t{pos_val}\t{out_info}\n")
    out_f.close()


if __name__ == '__main__':
    in_file_dict = __grab_all_input(IN_DIR, IN_ID)
    print(f"** All input files: {len(in_file_dict)}")
    for in_file in in_file_dict:
        print(in_file)
    print("** Constructing merged tree-dict...")
    merge_tree_dict = dict()
    for in_file, in_id in in_file_dict.items():
        print(f"** Processing {in_file} as {in_id}")
        __construct_tree_dict(in_file, in_id, merge_tree_dict)
    out_dict = __create_out_dict(merge_tree_dict.keys(), IN_ID, OUT_DIR, OUT_ID)
    for type_tag, out_file in out_dict.items():
        print(f"** Writing {out_file}...")
        __write_out_file(out_file, merge_tree_dict[type_tag])
    print("** All operations complete")
