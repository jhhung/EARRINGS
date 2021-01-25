import os
import sys
from path import *
PWD = os.path.dirname(os.path.realpath(__file__))

READS_NUM = 10000
INSERT_MEAN = 150

IDX_PREFIX = HG38_CHR1_EARRINGS_IDX
OUTPUT_DIRS = [os.path.join(PWD, "benchmark_no_adapter_outputs"),
               os.path.join(PWD, "benchmark_adapter_outputs"),
               os.path.join(PWD, "benchmark_adapter_lowcom_outputs")]

AD1_LOW = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
AD2_LOW = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

# Generate reads without/with normal adapters/ with low complexity adapter
commands = ["python3 {} -n {} -m {} -ad1 '' -ad2 '' -q 1 -t 1 -a 1".format(RUN_SIM_EXE, READS_NUM, INSERT_MEAN),
            "python3 {} -n {} -m {} -q 1 -t 1 -a 1".format(RUN_SIM_EXE, READS_NUM, INSERT_MEAN),
            "python3 {} -n {} -m {} -ad1 {} -ad2 {} -q 1 -t 1 -a 1".format(RUN_SIM_EXE, READS_NUM, INSERT_MEAN, AD1_LOW, AD2_LOW)]

if __name__ == "__main__":
    for i in range(0, 3):
        OUTPUT_DIR = OUTPUT_DIRS[i]
        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)

        OUTPUT_DATA_DIR = os.path.join(OUTPUT_DIR, 'data')
        if not os.path.exists(OUTPUT_DATA_DIR):
            os.makedirs(OUTPUT_DATA_DIR)
        for sd in range(18, 36):
            output_files = ["{}/{}_pe_adapter_atropos.txt".format(OUTPUT_DIR, str(sd)),
                            "{}/{}_se_adapter_atropos.txt".format(OUTPUT_DIR, str(sd)),
                            "{}/{}_pe_adapter_fastp.txt".format(OUTPUT_DIR, str(sd)),
                            "{}/{}_se_adapter_fastp.txt".format(OUTPUT_DIR, str(sd)),
                            "{}/{}_se_adapter_earrings_sen.txt".format(OUTPUT_DIR, str(sd)),
                            "{}/{}_pe_adapter_earrings_sen.txt".format(OUTPUT_DIR, str(sd)),
                            "{}/{}_se_adapter_earrings.txt".format(OUTPUT_DIR, str(sd)),
                            "{}/{}_pe_adapter_earrings.txt".format(OUTPUT_DIR, str(sd))]

            for f in output_files:
                if os.path.isfile(f):
                    os.remove(f)

            for seed in range(10):
                simu_fname =  "mean_{}_std_{}_seed_{}".format(str(INSERT_MEAN), str(sd), str(seed))
                simu_fpath = os.path.join(OUTPUT_DATA_DIR, simu_fname)
                simu_log = os.path.join(OUTPUT_DATA_DIR, "{}.log".format(simu_fname))
                simu_read_info = os.path.join(OUTPUT_DATA_DIR, "{}.read.info".format(simu_fname))
                GEN1 = "{}_1.fq".format(simu_fpath)
                GEN2 = "{}_2.fq".format(simu_fpath)
                ANS1 = "{}_ans_1.fq".format(simu_fpath)
                ANS2 = "{}_ans_2.fq".format(simu_fpath)

                run_all = []
            
                run_all.append("atropos detect -pe1 {} -pe2 {} | tee -a {}".format(GEN1, GEN2, output_files[0]))
                run_all.append("atropos detect -se {} | tee -a {}".format(GEN1, output_files[1]))
                # We cout the the detected adapters in fastp
                run_all.append("fastp -i {} -I {} --detect_adapter_for_pe | tee -a {}".format(GEN1, GEN2, output_files[2]))
                run_all.append("fastp -i {} | tee -a {}".format(GEN1, output_files[3]))
                
                # Sensitive mode: is able to successfully detect adapter even when the adapter contamination ratio is low.
                run_all.append("{} single -p {} --sensitive -1 {} | tee -a {}".format(EARRINGS_EXE, IDX_PREFIX, GEN1, output_files[4]))
                run_all.append("{} paired -1 {} -2 {} --sensitive -o ./earrings_pe_sen | tee -a {}".format(EARRINGS_EXE, GEN1, GEN2, output_files[5]))
                # Non-sensitive mode
                run_all.append("{} single -p {} -1 {} | tee -a {}".format(EARRINGS_EXE, IDX_PREFIX, GEN1, output_files[6]))
                run_all.append("{} paired -1 {} -2 {} -o ./earrings_pe | tee -a {}".format(EARRINGS_EXE, GEN1, GEN2, output_files[7]))
                
                print("dev: {}, seed: {}".format(str(sd), str(seed)), flush=True)
                command = commands[i] + " -r {} -std {} -s {} -O {} -o {} {}".format(HG38_CHR1_REF, str(sd), simu_fname, OUTPUT_DATA_DIR, simu_log, str(seed))
                print(command, flush=True)
                os.system(command)

                for run_file in run_all:
                    os.system(run_file)
                os.remove(GEN1)
                os.remove(GEN2)
                os.remove(ANS1)
                os.remove(ANS2)
                os.remove(simu_read_info)
