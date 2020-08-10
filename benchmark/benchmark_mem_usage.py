import os
import sys
import glob
import numpy as np
import gzip
import shutil

from prog_out import *
from run_trimmers import *
from path import *
from eval.pe_evaluation import *
from eval.se_evaluation import *
from eval.cal_confusion import *
from eval.order import *

PWD = os.path.dirname(os.path.realpath(__file__))
RUN_TRIMMERS = glob.glob(os.path.join(RUN_TRIMMER_DIR, "run_*"))

AD1_SIM = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCACCTAATCTCGTATGCCGTCTTCTGCTTG"
AD2_SIM = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"
AD1 = AD1_SIM[:32]
AD2 = AD1_SIM[:32]

READS_NUM = 10000000
READ_LEN = 100
INSERT_MEAN = 100
INSERT_DEV = 50

IDX = HG38_CHR1_EARRINGS_IDX

MEM_DIR = os.path.join(PWD, 'benchmark_mem')
OUTPUT_DATA_DIR = os.path.join(MEM_DIR, 'data')
OUTPUT_DIR = [os.path.join(PWD, 'benchmark_mem', 'benchmark_mem_se'), 
              os.path.join(PWD, 'benchmark_mem', 'benchmark_mem_pe')]
if __name__ == "__main__":
    if not os.path.exists(MEM_DIR):
        os.makedirs(MEM_DIR)
    if not os.path.exists(OUTPUT_DATA_DIR):
        os.makedirs(OUTPUT_DATA_DIR)

    for seed in range(0, 3):
        for i in range(2):
            dir_name = OUTPUT_DIR[i] + "_" + str(seed)
            if not os.path.exists(dir_name):
                os.makedirs(dir_name)
        simu_fname =  "mean_{}_std_{}_seed_{}".format(str(INSERT_MEAN), str(INSERT_DEV), str(seed))
        simu_fpath = os.path.join(OUTPUT_DATA_DIR, simu_fname)
        simu_log = os.path.join(OUTPUT_DATA_DIR, "{}.log".format(simu_fname))
        simu_read_info = os.path.join(OUTPUT_DATA_DIR, "{}.read.info".format(simu_fname))
        GEN1 = "{}_1.fq".format(simu_fpath)
        GEN2 = "{}_2.fq".format(simu_fpath)
        print("mean: {}, dev: {}, seed: {}".format(str(INSERT_MEAN), str(INSERT_DEV), str(seed)), flush=True)
        command = "python3 {} {} -r {} -n {} -m {} -std {} -q 1 -t 1 -a 1 -o {} -s {} -O {}"\
                .format(RUN_SIM_EXE, str(seed), HG38_CHR1_REF, READS_NUM, INSERT_MEAN, INSERT_DEV, simu_log, simu_fname, OUTPUT_DATA_DIR)
        os.system(command)
        
        for thread in [1, 2, 4, 8, 16, 32]:
            for trimmer in RUN_TRIMMERS:
                run_SE(trimmer, GEN1, AD1, IDX, thread, 0, 50, 0, 0, 1)
                run_PE(trimmer, GEN1, GEN2, AD1, AD2, thread, 0, 0.7, 0.9, 0.8, 0, 0, 1)

        os.system("mv *se_mm* {}".format(OUTPUT_DIR[0] + "_" + str(seed)))
        os.system("mv *pe_mm* {}".format(OUTPUT_DIR[1] + "_" + str(seed)))
        
        os.remove(GEN1)
        os.remove(GEN2)
        os.remove(ANS1)
        os.remove(ANS2)
        os.remove(simu_read_info)
        