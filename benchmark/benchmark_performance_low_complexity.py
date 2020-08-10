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
THREAD_NUM = 32

# we use low complexity to reproduce data similar to SRR529095
AD1_SIM = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
AD2_SIM = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

AD1 = AD1_SIM[:32]
AD2 = AD2_SIM[:32]

READS_NUM = 10000000
READ_LEN = 100
INSERT_MEAN = 100

SE_RECORD = os.path.join(PWD, "se_record.txt")
PE_RECORD = os.path.join(PWD, "pe_record.txt")

IDX = HG38_CHR1_EARRINGS_IDX

NUM_EXP_GRPS = 1  # total: 11 groups
NUM_EXP_REP = 3

START_STD = 50
STD_INCRE = 5
END_STD = 51

OUTPUT_DIR = os.path.join(PWD, 'benchmark_per_low_complexity')
OUTPUT_DATA_DIR = os.path.join(OUTPUT_DIR, 'data')
fname_names = [os.path.join(OUTPUT_DIR, 'se_per_order.txt'),
               os.path.join(OUTPUT_DIR, 'pe_per_order.txt')]
fname_results = [os.path.join(OUTPUT_DIR, 'se_per_results.npy'), 
                 os.path.join(OUTPUT_DIR, 'pe_per_results.npy')]
fname_matrices = [os.path.join(OUTPUT_DIR, 'se_per_matrices.npy'),
                  os.path.join(OUTPUT_DIR, 'pe_per_matrices.npy')]
fname_trim = [os.path.join(OUTPUT_DIR, 'se_per_trim.npy'), 
              os.path.join(OUTPUT_DIR, 'pe_per_trim.npy')]

if __name__ == "__main__":
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    if not os.path.exists(OUTPUT_DATA_DIR):
        os.makedirs(OUTPUT_DATA_DIR)

    se_trim_results = np.zeros((len(SE_PROG_SORT_OUT), NUM_EXP_GRPS, NUM_EXP_REP, 2))
    pe_trim_results = np.zeros((len(PE_PROG_SORT_OUT) + 1, NUM_EXP_GRPS, NUM_EXP_REP, 2))

    se_per_results = np.zeros((len(SE_PROG_SORT_OUT), NUM_EXP_GRPS, NUM_EXP_REP, 4))
    pe_per_results = np.zeros((len(PE_PROG_SORT_OUT) + 1, NUM_EXP_GRPS, NUM_EXP_REP, 4))

    se_per_matrices = np.zeros((len(SE_PROG_SORT_OUT), NUM_EXP_GRPS, NUM_EXP_REP, 5))
    pe_per_matrices = np.zeros((len(PE_PROG_SORT_OUT) + 1, NUM_EXP_GRPS, NUM_EXP_REP, 5))

    for std_idx, std in enumerate(range(START_STD, END_STD, STD_INCRE)):
        for seed_idx, seed in enumerate(range(NUM_EXP_REP)):
            simu_fname =  "mean_{}_std_{}_seed_{}".format(str(INSERT_MEAN), str(std), str(seed))
            simu_fpath = os.path.join(OUTPUT_DATA_DIR, simu_fname)
            simu_log = os.path.join(OUTPUT_DATA_DIR, "{}.log".format(simu_fname))
            simu_read_info = os.path.join(OUTPUT_DATA_DIR, "{}.read.info".format(simu_fname))
            GEN1 = "{}_1.fq".format(simu_fpath)
            GEN2 = "{}_2.fq".format(simu_fpath)
            ANS1 = "{}_ans_1.fq".format(simu_fpath)
            ANS2 = "{}_ans_2.fq".format(simu_fpath)
            PE_PROG_LIST = append_AKA(GEN1, GEN2, PE_PROG_SORT_OUT)
            PE_PROG_LIST = append_AKA(os.path.join(PWD, "{}_1.fq".format(simu_fname)), os.path.join(PWD, "{}_2.fq".format(simu_fname)), PE_PROG_SORT_OUT)
            print("dev: {}, seed: {}".format(str(std), str(seed)), flush=True)
            command = "python3 {} {} -r {} -n {} -m {} -std {} -ad1 {} -ad2 {} -q 1 -t 1 -a 1 -o {} -s {} -O {}"\
                    .format(RUN_SIM_EXE, seed, HG38_CHR1_REF, READS_NUM, INSERT_MEAN, std, AD1_SIM, AD2_SIM, simu_log, simu_fname, OUTPUT_DATA_DIR)
            os.system(command)
            # Note that here, atropos and fastp are not able to auto-detect adapters
            # Therefore, we disable auto-detect for atropos and fastp here.
            for trimmer in RUN_TRIMMERS:
                # trim without auto-detect
                if trimmer[-5:-3] == "FP":
                    run_SE(trimmer, GEN1, AD1, IDX, THREAD_NUM, disable_default_ad=1, enable_polyX_fq=1)
                    run_PE(trimmer, GEN1, GEN2, AD1, AD2, THREAD_NUM, disable_default_ad=1, enable_polyX_fq=1)
                else:
                    run_SE(trimmer, GEN1, AD1, IDX, THREAD_NUM)
                    run_PE(trimmer, GEN1, GEN2, AD1, AD2, THREAD_NUM)
                # trim with auto-detect
                # if trimmer[-5:-3] == "AP":
                #     p = os.path.join(PWD, AP[1])
                #     os.system("atropos detect -se {} | tee {}".format(GEN1, p))
                #     f = open(p, "r")
                #     kmer_largest, kmer_match_largest = atropos_retrieve_ad(f)
                #     run_SE(trimmer, GEN1, kmer_largest[:32], IDX, THREAD_NUM)

                #     p = os.path.join(PWD, AP[2])
                #     os.system("atropos detect -pe1 {} -pe2 {} | tee {}".format(GEN1, GEN2, p))
                #     f = open(p, "r")
                #     kmer_largest1, kmer_match_largest = atropos_retrieve_ad(f)
                #     kmer_largest2, kmer_match_largest = atropos_retrieve_ad(f)
                #     run_PE(trimmer, GEN1, GEN2, kmer_largest1[:32], kmer_largest2[:32], THREAD_NUM)
                # else:
                #     run_SE(trimmer, GEN1, AD1, IDX, THREAD_NUM, disable_default_ad=1)
                #     run_PE(trimmer, GEN1, GEN2, AD1, AD2, THREAD_NUM, disable_default_ad=1)

            # unzip SeqPurge
            for i in range(2, 4):
                SP_path = os.path.join(PWD, SP[i])
                if os.path.isfile(SP_path):
                    unzip(SP_path, os.path.join(PWD, SP_UNZIP[i]))
            # concat the dropped reads
            concat_TR_PE(PWD)
            # sort atropos, fastp SE & PE
            for i in range(1, 4):
                AP_path = os.path.join(PWD, AP[i])
                FP_path = os.path.join(PWD, FP[i])
                SP_path = os.path.join(PWD, SP_UNZIP[i])
                TR_path = os.path.join(PWD, TR[i])
                
                if os.path.isfile(AP_path): 
                    ordering(AP_path, os.path.join(PWD, AP_SORT[i]))
                if os.path.isfile(FP_path):
                    ordering(FP_path, os.path.join(PWD, FP_SORT[i]))
                if os.path.isfile(SP_path):
                    ordering(SP_path, os.path.join(PWD, SP_SORT[i]))
                if os.path.isfile(TR_path):
                    ordering(TR_path, os.path.join(PWD, TR_SORT[i]))


            print("==SE==", flush=True)
            for idx, f_name in enumerate(SE_PROG_SORT_OUT):
                f_path = os.path.join(PWD, f_name[1])
                if not os.path.isfile(f_path):
                    continue
                print(f_name[0], flush=True)
                tp, tn, fp, fn, num_overtrim, num_undertrim = se_cmp_results([f_path], [ANS1], READ_LEN)
                se_trim_results[idx][std_idx][seed_idx] = [num_overtrim, num_undertrim]
                se_per_results[idx][std_idx][seed_idx] = [tp, tn, fp, fn]
                ACC, SEN, SPC, PPV, MCC = cal_matrix(tp, tn, fp, fn)
                se_per_matrices[idx][std_idx][seed_idx] = [ACC, SEN, SPC, PPV, MCC]

            print("==PE==", flush=True)
            for idx, f_name in enumerate(PE_PROG_LIST):
                f1_path = os.path.join(PWD, f_name[2])
                f2_path = os.path.join(PWD, f_name[3])
                if not os.path.isfile(f1_path) or not os.path.isfile(f2_path):
                    continue
                    
                print(f_name[0], flush=True)
                tp, tn, fp, fn, num_overtrim, num_undertrim = pe_cmp_results([f1_path, f2_path], [ANS1, ANS2], READ_LEN)
                pe_trim_results[idx][std_idx][seed_idx] = [num_overtrim, num_undertrim]
                pe_per_results[idx][std_idx][seed_idx] = [tp, tn, fp, fn]
                ACC, SEN, SPC, PPV, MCC = cal_matrix(tp, tn, fp, fn)
                pe_per_matrices[idx][std_idx][seed_idx] = [ACC, SEN, SPC, PPV, MCC]
            
            os.remove(simu_read_info)
            os.remove(GEN1)
            os.remove(GEN2)
            os.remove(ANS1)
            os.remove(ANS2)
            
    se_order = open(fname_names[0], "w")
    pe_order = open(fname_names[1], "w")
    for f_name in SE_PROG_SORT_OUT:
        se_order.write(f_name[0] + ',')
    se_order.close()

    for f_name in PE_PROG_LIST:
        pe_order.write(f_name[0] + ',')
    pe_order.close()

    np.save(fname_trim[0], se_trim_results)
    np.save(fname_trim[1], pe_trim_results)
    np.save(fname_results[0], se_per_results)
    np.save(fname_matrices[0], se_per_matrices)
    np.save(fname_results[1], pe_per_results)
    np.save(fname_matrices[1], pe_per_matrices)

