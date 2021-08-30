import sys
import bioparse as bp


# Input Parameters
REF_FILE = sys.argv[1]
IN_FILE = sys.argv[2]

# Global Variables
REF_TAG_DICT = dict()
COM_TAG_DICT = dict()


def __add_to_tag_dict(cur_items, tag_dict):
    cur_id = cur_items[2]
    if cur_id not in tag_dict:
        tag_dict[cur_id] = 0
    tag_dict[cur_id] += 1


if __name__ == '__main__':
    print("# Reference: {}".format(REF_FILE))
    print("# Input: {}".format(IN_FILE))
    ref_fp = bp.FileParser(REF_FILE)
    in_fp = bp.FileParser(IN_FILE)
    for ref_items in ref_fp.iter():
        if in_fp.term:
            __add_to_tag_dict(ref_items, REF_TAG_DICT)
            continue
        cmp_res = bp.cmp_lines(ref_items, in_fp.get_line())
        while cmp_res > 0:
            in_fp.next()
            if in_fp.term:
                break
            cmp_res = bp.cmp_lines(ref_items, in_fp.get_line())
        if cmp_res == 0:
            __add_to_tag_dict(ref_items, COM_TAG_DICT)
        else:
            __add_to_tag_dict(ref_items, REF_TAG_DICT)
    ref_fp.close()
    in_fp.close()
    print("# Statistics")
    tot_count = 0
    for tag in sorted(REF_TAG_DICT.keys()):
        print("RefOnly | Tag: {} | Count: {}".format(tag, REF_TAG_DICT[tag]))
        tot_count += REF_TAG_DICT[tag]
    print("RefOnly | Total: {}".format(tot_count))
    tot_count = 0
    for tag in sorted(COM_TAG_DICT.keys()):
        print("Common | Tag: {} | Count: {}".format(tag, COM_TAG_DICT[tag]))
        tot_count += COM_TAG_DICT[tag]
    print("Common | Total: {}".format(tot_count))
