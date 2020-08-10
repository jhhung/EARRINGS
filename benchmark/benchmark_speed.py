import os
import glob
from prog_out import *
from run_trimmers import *
from timeit import default_timer as timer
import numpy as np
from path import *

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

OUTPUT_DIR = os.path.join(PWD, 'benchmark_speed')
OUTPUT_DATA_DIR = os.path.join(OUTPUT_DIR, 'data')

THREADS = [1, 2, 4, 8, 16, 32]
IDX = HG38_CHR1_EARRINGS_IDX
OUTPUT_DIR = os.path.join(PWD, 'benchmark_speed')
OUTPUT_DATA_DIR = os.path.join(OUTPUT_DIR, 'data')
fname_results = [os.path.join(OUTPUT_DIR, "SE_times.npy"), os.path.join(OUTPUT_DIR, "PE_times.npy")]

if __name__ == "__main__":
    SE_times = np.zeros((len(RUN_TRIMMERS), len(THREADS), 10))
    PE_times = np.zeros((len(RUN_TRIMMERS), len(THREADS), 10))
    # AR, CA, FP, PE, SK, AP, AKA
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    if not os.path.exists(OUTPUT_DATA_DIR):
        os.makedirs(OUTPUT_DATA_DIR)
    
    for seed in range(10):
        print("Seed: {}".format(seed), flush=True)
        simu_fname =  "mean_{}_std_{}_seed_{}".format(str(INSERT_MEAN), str(INSERT_DEV), str(seed))
        simu_fpath = os.path.join(OUTPUT_DATA_DIR, simu_fname)
        simu_log = os.path.join(OUTPUT_DATA_DIR, "{}.log".format(simu_fname))
        simu_read_info = os.path.join(OUTPUT_DATA_DIR, "{}.read.info".format(simu_fname))
        GEN1 = "{}_1.fq".format(simu_fpath)
        GEN2 = "{}_2.fq".format(simu_fpath)
        ANS1 = "{}_ans_1.fq".format(simu_fpath)
        ANS2 = "{}_ans_2.fq".format(simu_fpath)
        command = "python3 {} {} -r {} -n {} -m {} -std {} -ad1 {} -ad2 {} -q 1 -t 1 -a 1 -o {} -s {} -O {}"\
                .format(RUN_SIM_EXE, str(seed), HG38_CHR1_REF, READS_NUM, INSERT_MEAN, INSERT_DEV, AD1, AD2, simu_log, simu_fname, OUTPUT_DATA_DIR)
        os.system(command)
        
        for i, trimmer in enumerate(RUN_TRIMMERS):
            for j, thread in enumerate(THREADS):
                print(trimmer)
                SE_start = timer()
                run_SE(trimmer, GEN1, AD1, IDX, thread)
                SE_end = timer()
                SE_times[i][j][seed] = SE_end - SE_start
                print(SE_end - SE_start)
                
                PE_start = timer()
                run_PE(trimmer, GEN1, GEN2, AD1, AD2, thread)
                PE_end = timer()
                PE_times[i][j][seed] = PE_end - PE_start
                print(PE_end - PE_start)
        
        os.remove(GEN1)
        os.remove(GEN2)
        os.remove(ANS1)
        os.remove(ANS2)
        os.remove(simu_read_info)
    
    f = open(os.path.join(OUTPUT_DIR, "speed_order.txt"), "w")
    for trimmer in RUN_TRIMMERS:
        f.write(trimmer + '\n')
    f.close()

    np.save(fname_results[0], SE_times)
    np.save(fname_results[1], PE_times)


