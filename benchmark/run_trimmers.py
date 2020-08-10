import os
import re

def run_SE(exe, in1, ad1, earrings_idx, thread_num=1, out_min_len=0, earrings_min_len=50, disable_default_ad=0, enable_polyX_fq=0, record_mm=0):
    command = "{} -a1 {} -i {} -t {} -m {} -mp {} -da {} -px {} -rc {}".format(in1, ad1, earrings_idx, thread_num, out_min_len, earrings_min_len, disable_default_ad, enable_polyX_fq, record_mm)
    os.system("python3 {} {}".format(exe, command))

def run_PE(exe, in1, in2, ad1, ad2, thread_num=1, out_min_len=0, match_rate=0.7, seq_cmp_rate=0.9, adapter_cmp_rate=0.8, disable_default_ad=0, enable_polyX_fq=0, record_mm=0):
    command = "{} -f2 {} -a1 {} -a2 {} -t {} -m {} -mr {} -sc {} -ac {} -da {} -px {} -rc {}".format(in1, in2, ad1, ad2, thread_num, out_min_len, match_rate, seq_cmp_rate, adapter_cmp_rate, disable_default_ad, enable_polyX_fq, record_mm)
    os.system("python3 {} {}".format(exe, command))

def run_aligner_SE(exe, idx, input_f, output_sam, thread=64):
    os.system("{} -x {} -U {} -S {} -p {}".format(exe, idx, input_f, output_sam, thread))

def run_aligner_PE(exe, idx, input_f1, input_f2, output_sam, thread=64):
    os.system("{} -x {} -1 {} -2 {} -S {} -p {}".format(exe, idx, input_f1, input_f2, output_sam, thread))

def atropos_retrieve_ad(f):
    kmer_largest = ""
    kmer_match_largest = 0

    line = f.readline()
    if line[:8] == "detector":
        f.readline()

    for _ in range(5):
        f.readline()
    line = f.readline()
    found = re.findall('[0-9]+', line)
    if len(found) == 0:
        return kmer_largest, kmer_match_largest
    records = int(found[0])
    for i in range(records):
        kmer = f.readline().split(' ')[-1][:-1]
        f.readline()
        line = f.readline()
        found = re.findall('Name', line)
        if len(found) != 0:
            found = []
            while len(found) == 0:
                line = f.readline()
                found = re.findall('Known', line)
            for _ in range(3):
                f.readline()
        line = f.readline()
        kmer_match = int(line.split(' ')[-1][:-1])
        if kmer_match > kmer_match_largest:
            kmer_largest = kmer
            kmer_match_largest = kmer_match

    return kmer_largest, kmer_match_largest