import os
import sys
import glob

from prog_out import *
from run_trimmers import *
from path import *
from eval.cal_mean_std import *
from eval.cut_length import *

PWD = os.path.dirname(os.path.realpath(__file__))
RUN_TRIMMERS = glob.glob(os.path.join(RUN_TRIMMER_DIR, "run_*"))
ALIGNER_THREAD = 64
THREAD_NUM = 32

def run_SE_align(prog_out, in1, ad1, earrings_idx, aligner_idx, aligner_exe, file_prefix, is_low_complexity, min_len=50):
    for trimmer in RUN_TRIMMERS:
        run_SE(trimmer, in1, ad1, earrings_idx, THREAD_NUM, OUTPUT_MIN, min_len, enable_polyX_fq=1)
    
    TR_path = os.path.join(PWD, TR[1])
    if os.path.isfile(TR_path):
        cut_length([TR_path], OUTPUT_MIN)

    for fs in prog_out:
        f = os.path.join(PWD, fs[1])
        if os.path.isfile(f):
            print(fs[0], flush=True)
            base_name = os.path.join(PWD, "{}_{}".format(fs[0], file_prefix))
            output_sam = "{}.sam".format(base_name)
            run_aligner_SE(aligner_exe, aligner_idx, f, output_sam)
            os.remove(output_sam) 

    print("RAW", flush=True)
    output_sam = os.path.join(PWD, "{}_raw.sam".format(file_prefix))
    os.system("{} -x {} -U {} -S {} -p {}".format(aligner_exe, aligner_idx, in1, output_sam, ALIGNER_THREAD))
    os.remove(output_sam)
            

def run_PE_align(prog_out, in1, in2, ad1, ad2, aligner_idx, aligner_exe, file_prefix, ref_fa):
    for trimmer in RUN_TRIMMERS:
        run_PE(trimmer, in1, in2, ad1, ad2, THREAD_NUM, OUTPUT_MIN)

    for i in range(2, 4):
        SP_path = os.path.join(PWD, SP[i])
        PE_path = os.path.join(PWD, PE[i])
        PE_UNZIP_path = os.path.join(PWD, PE_UNZIP[i])
        if os.path.isfile(SP_path):
            unzip(SP_path, os.path.join(PWD, SP_UNZIP[i]))
        if os.path.isfile(PE_path):
            rename_PE(PE_path, PE_UNZIP_path)

    TR_path1 = os.path.join(PWD, TR[2])
    TR_path2 = os.path.join(PWD, TR[3])
    if os.path.isfile(TR_path1) and os.path.isfile(TR_path2):
        cut_length([TR_path1, TR_path2], OUTPUT_MIN)

    for fs in prog_out:
        f1 = os.path.join(PWD, fs[2])
        f2 = os.path.join(PWD, fs[3])

        if os.path.isfile(f1) and os.path.isfile(f2):
            print(fs[0], flush=True)
            base_name = os.path.join(PWD, "{}_{}".format(fs[0], file_prefix))
            output_sam = "{}.sam".format(base_name)
            run_aligner_PE(aligner_exe, aligner_idx, f1, f2, output_sam)
            sorted_sam = cal_insert_picard(base_name, ref_fa, ads)
            os.remove(output_sam)
            os.remove(sorted_sam)
    
    print("RAW", flush=True)
    output_sam = os.path.join(PWD, "{}_raw.sam".format(file_prefix))
    os.system("{} -x {} -1 {} -2 {} -S {} -p {}".format(aligner_exe, aligner_idx, in1, in2, output_sam, ALIGNER_THREAD))    
    os.remove(output_sam)

# We disable auto-detect for atropos and fastp here.
# Since in this benchmarking, we focus on the ability of trimmers to trim real data.

### modify the following path to correct path ###
# ref_fa: path to the reference genome
# in1/in2: path to the input fastq file(s)
# er_idx: path to the prebuild index of EARRINGS
# idx: path to the prebuild index of hisat2/bowtie2

# SRR529095: RNA IP-Seq of Argonaute proteins from Homo sapiens
# aligner: hisat2
ad1 = "AAAAAAAAAAAAAAAAAAAA"  # low complexity adapters
# ref_fa = ...
# in1 = ...
# er_idx = ...
# idx = ...
OUTPUT_MIN = 18
run_SE_align(SE_PROG_OUT, in1, ad1, er_idx, idx, HISAT2_EXE, "SRR529095", 1, min_len=50)

# SRR014966: miRNA of C.elegans
# aligner: bowtie2
ad1 = "TCGTATGCCGTCTTCTGCT"
# ref_fa = ...
# in1 = ... 
# er_idx = ...
# idx = ...
OUTPUT_MIN = 18
run_SE_align(SE_PROG_OUT, in1, ad1, er_idx, idx, BOWTIE2_EXE, "SRR014866", 0, min_len=18)

# SRR330569: RNA-Seq of Gonads and Carcasses in D. simulans and D. pseudoobscura
# aligner: hisat2
ad1 = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCG"
ad2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGAT"
# ref_fa = ...
# in1 = ... 
# in2 = ...
# idx = ... 
OUTPUT_MIN = 18
PE_PROG_OUT_AKA = append_AKA("SRR330569.3_1.fastq", "SRR330569.3_2.fastq", PE_PROG_OUT, "fastq")
run_PE_align(PE_PROG_OUT_AKA, in1, in2, ad1, ad2, idx, HISAT2_EXE, "SRR330569", ref_fa, [ad1, ad2])

# SRR5000681: ATAC-Seq of Early Embryo replicate in C.elegans
# aligner: bowtie2
ad1 = "CTGTCTCTTATACACATCTCCGAGCCCACGAG"
ad2 = "CTGTCTCTTATACACATCTGACGCTGCCGACG"
# ref_fa = ...
# in1 = ...
# in2 = ... 
# idx = ... 
OUTPUT_MIN = 18
PE_PROG_OUT_AKA = append_AKA("SRR5000681.1_1.fastq", "SRR5000681.1_2.fastq", PE_PROG_OUT, "fastq")
run_PE_align(PE_PROG_OUT_AKA, in1, in2, ad1, ad2, idx, BOWTIE2_EXE, "SRR5000681", ref_fa, [ad1, ad2])
