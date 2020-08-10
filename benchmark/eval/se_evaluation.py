import sys
from argparse import ArgumentParser
from .fastq import Fastq

def se_cmp_results(input_fnames, ans_fnames, read_len):
    input_f = open(input_fnames[0], "r")
    ans_f = open(ans_fnames[0], "r")
    TP = 0  # reads contain adapter & trimmer trims to the right length
    TN = 0  # reads don't contain adapter & trimmer does not trim
    FP = 0  # reads contain adapter & trimmer trims too much adapter
    FN = 0  # reads contain adapter & trimmer trims too little adapter

    input_fq1 = Fastq()
    ans_fq1 = Fastq()
    is_eof = False

    num_overtrim = 0
    num_undertrim = 0

    while not input_fq1.get_obj(input_f) and not ans_fq1.get_obj(ans_f):
        diff = len(ans_fq1.seq) - len(input_fq1.seq)
        if diff > 0:
            # ans is longer than input -> overtrim
            num_overtrim += diff
        else:
            num_undertrim -= diff

        while input_fq1.name != ans_fq1.name:
            diff = len(ans_fq1.seq)
            if diff > 0:
                # ans is longer than input -> overtrim
                num_overtrim += diff

            if len(ans_fq1.seq) == 0:
                TP += 1
            else:
                # ans len > 0
                FP += 1

            is_eof = ans_fq1.get_obj(ans_f)
            if is_eof:
                break

        assert not is_eof, "Opps, something wrong happens!"

        # seq
        if len(input_fq1.seq) == len(ans_fq1.seq) == read_len:
            TN += 1
        elif len(input_fq1.seq) == len(ans_fq1.seq):
            TP += 1
        elif len(input_fq1.seq) > len(ans_fq1.seq):
            FN += 1
        elif len(input_fq1.seq) < len(ans_fq1.seq):
            FP += 1

    # cal_matrix(TP, TN, FP, FN)
    return TP, TN, FP, FN, num_overtrim, num_undertrim

if __name__ == "__main__":
    parser = ArgumentParser(description="A program to evaluate how adapter trimmer works.")
    parser.add_argument("input_fnames", type=str, nargs=1, help="Two input files used to evaluate")
    parser.add_argument("ans_fnames", type=str, nargs=1, help="Two answer files used to evaluate")
    parser.add_argument("read_length", type=int, nargs=1, help="The length of the read")
    parser.add_argument("-o", "--output", help="output filename's prefix", dest="o_Fname", default="default")

    args = parser.parse_args()
    se_cmp_results(args.input_fnames, args.ans_fnames, args.read_length[0])


