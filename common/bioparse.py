"""Collection of classes and functions to aid parsing files."""

# Import Modules
import os.path
from datetime import date
import time
import gzip
import pickle
from typing import Dict, List

# Key Variables
PROG: int = 1000000
NON_AUTO_CHR_DIC: Dict[str, int] = {"x": 23, "y": 24, "m": 25}
ORD_CHR_LIST: List[str] = [*(map(str, range(1, 23))), "X", "Y", "M"]


class FileParser(object):
    """Mitigate file parsing."""

    def __init__(self, filepath: str, is_skip_hd: bool = True) -> None:
        """Open file and skip header by default."""
        self.obj = None
        self.line: str = None
        self.term: bool = True
        self.l_cnt: int = 0
        self.init_time = time.time()
        if os.path.getsize(filepath) > 0:
            self.term = False
            if filepath[-2:] == "gz":
                self.obj = gzip.open(filepath, "rt")
            else:
                self.obj = open(filepath, "r")
            self.line = next(self.obj)
            if is_skip_hd:
                while self.line and self.line[0] == "#":
                    self.next()
        return

    def get_line(self) -> str:
        return None if self.term else self.line

    def get_items(self) -> List[str]:
        return None if self.term else self.line.split()

    def next(self) -> None:
        if not self.term:
            try:
                self.l_cnt += 1
                self.line = next(self.obj)
            except StopIteration:
                self.line = None
                self.term = True
                self.close()
        return

    def close(self) -> None:
        if self.obj:
            self.obj.close()
        return

    def iter(self) -> List[str]:
        while not self.term:
            yield self.get_items()
            self.next()
        return

    def copy_hd(self, f) -> None:
        while self.line[0] == "#":
            f.write("{}\n".format(self.line.strip()))
            self.next()
        return

    def print_prog(self, custom_hd: str = None) -> None:
        print_prog(self.l_cnt, self.init_time, custom_hd)

    def long_next(self, skip: int) -> None:
        if not self.term:
            try:
                if skip > 1:
                    for i in range(skip - 1):
                        next(self.obj)
                self.line = next(self.obj)
            except StopIteration:
                self.line = None
                self.term = True
        return


def load_pickle(p_filepath: str, print_time: bool = True):
    if print_time:
        init_time = time.time()
    with open(p_filepath, "rb") as p_f:
        data_inst = pickle.load(p_f)
    if print_time:
        print_runtime(init_time)
    return data_inst


def dump_pickle(dump_data, dump_filepath: str, print_time: bool = True):
    if print_time:
        init_time = time.time()
    with open(dump_filepath, "wb") as dump_f:
        pickle.dump(dump_data, dump_f)
    if print_time:
        print_runtime(init_time)
    return


def chrom_to_int(chrom_str: str) -> int:
    if chrom_str.isdigit():
        chrom_num: int = int(chrom_str)
        if 1 <= chrom_num or chrom_num <= 22:
            return int(chrom_str)
    elif chrom_str.lower() in NON_AUTO_CHR_DIC:
        return NON_AUTO_CHR_DIC[chrom_str.lower()]
    return -1


def iter_chrom() -> str:
    for chr_key in ORD_CHR_LIST:
        yield chr_key
    return


def skip_hd(f) -> None:
    line: str = next(f)
    while line[0] == "#":
        line = next(f)
    return line


def write_vcf_hd(f, extra_fields: List[str] = None) -> None:
    f.write("##fileformat=VCFv4.3\n")
    f.write("##fileDate={}\n".format(date.today().strftime("%Y%m%d")))
    hd = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if extra_fields:
        hd = "{}\t{}".format(hd, "\t".join(extra_fields))
    f.write("{}\n".format(hd))


def line_to_chrpos(line) -> str:
    e: List[str] = line.split() if isinstance(line, str) else line
    return "{}-{}".format(e[0], e[1])


def cmp_lines(ref_l, other_l) -> int:
    ref_e: List[str] = ref_l.split() if isinstance(ref_l, str) else ref_l
    other_e: List[str] = other_l.split() if isinstance(other_l, str) else other_l
    ref_chr_num: int = chrom_to_int(ref_e[0][3:])
    other_chr_num: int = chrom_to_int(other_e[0][3:])
    if ref_chr_num <= 0 or other_chr_num <= 0:
        return None
    chr_diff: int = ref_chr_num - other_chr_num
    if chr_diff == 0:
        pos_diff: int = int(ref_e[1]) - int(other_e[1])
        if pos_diff > 0:
            return 1
        elif pos_diff < 0:
            return -1
        else:
            return 0
    elif chr_diff > 0:
        return 1
    else:
        return -1


def print_prog(l_cnt: int, init_time, custom_hd: str = None) -> None:
    if l_cnt % PROG == 0:
        if custom_hd:
            print(custom_hd, end=" ")
        print("Progress: {:>6.1f}m line | ".format(l_cnt / PROG), end="")
        print_runtime(init_time)
    return


def print_runtime(init_time, custom_hd: str = None) -> None:
    run_time = (time.time() - init_time) / 60
    if custom_hd:
        print(custom_hd, end=" ")
    print("Runtime: {:>7.2f}m".format(run_time))
