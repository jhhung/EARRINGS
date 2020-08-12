import os
import sys
import argparse
from path import *

PWD = os.path.dirname(os.path.realpath(__file__))

def build_all_required_index(num_threads=32):
    # hg38 chr1 for simulation benchmarking
    print("Prebuild index of hg38_chr1 for EARRINGS...", flush=True)
    command = "{} build -r {} -p {}".format(EARRINGS_EXE, HG38_CHR1_REF, HG38_CHR1_EARRINGS_IDX)
    os.system(command)
    print("Finish prebuilding index of hg38_chr1 for EARRINGS...", flush=True)
    
    # hg38 for SRR529095
    print("Prebuild index of hg38 for EARRINGS...", flush=True)
    command = "{} build -r {} -p {}".format(EARRINGS_EXE, HG38_REF, HG38_EARRINGS_IDX)
    os.system(command)
    print("Finish prebuilding index of hg38 for EARRINGS...", flush=True)

    print("Prebuild index of hg38 for hisat2...", flush=True)
    command = "{} {} {} -p {}".format(HISAT2_BUILD_EXE, HG38_REF, HG38_HISAT2_IDX, num_threads)
    os.system(command)
    print("Finish prebuilding index of hg38 for hisat2...", flush=True)

    # C.elegans for SRR014966/SRR5000681
    print("Prebuild index of C.elegans for EARRINGS...", flush=True)
    command = "{} build -r {} -p {}".format(EARRINGS_EXE, C_ELEGANS_REF, C_ELEGAN_EARRINGS_IDX)
    os.system(command)
    print("Finish prebuilding index of C.elegans for EARRINGS...", flush=True)

    print("Prebuild index of C.elegans for bowtie2...", flush=True)
    command = "{} {} {} --threads {}".format(BOWTIE2_BUILD_EXE, C_ELEGANS_REF, C_ELEGAN_BOWTIE2_IDX, num_threads)
    os.system(command)
    print("Finish prebuilding index of C.elegans for bowtie2...", flush=True)

    # drosophila for SRR330569
    print("Prebuild index of drosophila for hisat2...", flush=True)
    command = "{} {} {} -p {}".format(HISAT2_BUILD_EXE, D_SIMULANS_REF, D_SIMULANS_HISAT2_IDX, num_threads)
    os.system(command)
    print("Finish prebuilding index of drosophila for hisat2...", flush=True)

if __name__ == "__main__":
    build_all_required_index()