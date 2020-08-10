import os
import sys
from path import *

if __name__ == "__main__":
    hg38_chr1_fw = HG38_CHR1_EARRINGS_IDX + "table"
    hg38_chr1_bw = HG38_CHR1_EARRINGS_IDX + "rc_table"
    if not os.path.isfile(hg38_chr1_fw) or not os.path.isfile(hg38_chr1_bw):
        print("Prebuild index of hg38_chr1 of EARRINGS...", flush=True)
        command = "{} build -r {} -p {}".format(EARRINGS_EXE, HG38_CHR1_REF, HG38_CHR1_EARRINGS_IDX)
        os.system(command)
        print("Finish prebuilding index of hg38_chr1 of EARRINGS...", flush=True)
    else:
        print("Prebuild index of hg38_chr1 of EARRINGS already exists. Skipping ...", flush=True)

    hg38_fw = HG38_EARRINGS_IDX + "table"
    hg38_bw = HG38_EARRINGS_IDX + "rc_table"
    print("Prebuild index of hg38 of EARRINGS...", flush=True)
    if not os.path.isfile(hg38_fw) or not os.path.isfile(hg38_bw):
        command = "{} build -r {} -p {}".format(EARRINGS_EXE, HG38_REF, HG38_EARRINGS_IDX)
        os.system(command)
        print("Finish prebuilding index of hg38 of EARRINGS...", flush=True)
    else:
        print("Prebuild index of hg38 of EARRINGS already exists. Skipping ...", flush=True)
    print("Finish!", flush=True)