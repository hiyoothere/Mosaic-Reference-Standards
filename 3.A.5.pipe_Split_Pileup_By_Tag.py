"""Split all chromosome merged parsed pileup information by control set tag."""


# Import Modules
import sys
sys.path.append(sys.path[0] + "/../common")
import glob
import bioparse as bp


# Input Parameters
IN_DIR = sys.argv[1]
IN_ID = sys.argv[2]
IN_DS = sys.argv[3]
OUT_DIR = sys.argv[4]
NO_REFHOM = sys.argv[5]


# Global Variables
IN_MOS = IN_ID.split("-")[0]
IN_ID = IN_ID if "N" in IN_DS else f"{IN_DS}.{IN_ID}"
NAIVE = False
ARG_TRUE_SET = {"true", "True", "TRUE"}
no_refhom = NO_REFHOM in ARG_TRUE_SET
if no_refhom:
    OUT_ID_TUP = ("nc_snv", "nc_ind", "pc_snv", "pc_ind")
else:
    OUT_ID_TUP = ("nc_low", "nc_snv", "nc_ind", "pc_snv", "pc_ind")
SAMP_TO_VAR = {"M1": ("V1", "V2"), "M2": ("V1", "V3", "V4"), "M3": ("V1", "V3", "V5")}


def __create_out_dict():
    out_dict = dict()
    for out_id in OUT_ID_TUP:
        out_file = f"{OUT_DIR}/{IN_ID}.parse.{out_id}.txt"
        print(out_file)
        out_dict[out_id] = open(out_file, "w")
    return out_dict


if __name__ == '__main__':
    in_file_list = glob.glob(f"{IN_DIR}/*{IN_ID}.parse.txt")
    if len(in_file_list) != 1:
        print("ERROR: Invalid input file")
        sys.exit()
    else:
        print(f"** Input: {in_file_list[0]}")
    out_dict = __create_out_dict()
    in_fp = bp.FileParser(in_file_list[0])
    in_var_tup = SAMP_TO_VAR[IN_MOS] if IN_MOS in SAMP_TO_VAR else tuple()
    for in_items in in_fp.iter():
        in_fp.print_prog()
        cur_id = in_items[-1]
        if cur_id[0] == "[":
            cur_id = cur_id[2:-2].split(";")[0]
        if NAIVE:
            out_l = "{}\t{}".format("\t".join(in_items[:-1]), cur_id)
        else:
            out_l = "{}\t{}\t{}".format("\t".join(in_items[:3]), "\t".join(in_items[4:9]), cur_id)
        if not no_refhom and "refhom" in cur_id:
            out_dict["nc_low"].write(f"{out_l}\n")
        elif cur_id.split("_")[0] in in_var_tup:
            if "ind" in cur_id:
                out_dict["pc_ind"].write(f"{out_l}\n")
            else:
                out_dict["pc_snv"].write(f"{out_l}\n")
        elif "V" not in cur_id:
            if "ind" in cur_id:
                out_dict["nc_ind"].write(f"{out_l}\n")
            else:
                out_dict["nc_snv"].write(f"{out_l}\n")
    for out_f in out_dict.values():
        out_f.close()
    print("** All operations complete")
