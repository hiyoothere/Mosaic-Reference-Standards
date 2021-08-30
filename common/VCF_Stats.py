import sys
import bioparse as bp


# Input Parameters
VCF_FILE = sys.argv[1]


# Global Variables
TAG_DICT = dict()


if __name__ == '__main__':
    print("# Input: {}".format(VCF_FILE))
    vcf_f = bp.FileParser(VCF_FILE)
    for vcf_items in vcf_f.iter():
        cur_tag = vcf_items[2]
        if cur_tag not in TAG_DICT:
            TAG_DICT[cur_tag] = 0
        TAG_DICT[cur_tag] += 1
    vcf_f.close()
    for cur_tag in sorted(TAG_DICT.keys()):
        print("{}: {}".format(cur_tag, TAG_DICT[cur_tag]))