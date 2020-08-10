#!/usr/bin/env python3
import os
from parser import parse_args

NAME = "TM"
EXE = "java -jar /usr/local/bin/trimmomatic-0.39.jar"
BASE = "trimmomatic_"
ADAPTER_FNAME = "adapters.fa"

args = parse_args()
command = EXE

f = open(ADAPTER_FNAME, "w")
f.write(">adapter/1\n")
f.write(args.a1)
is_paired = False
# ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>
#             :<simple clip threshold>:<minAdapterLength>:<keepBothReads>
if args.f2 is not None:
    is_paired = True
    BASE += "pe.fastq"
    command += " PE {} {} -baseout {}".format(args.f1, args.f2, BASE)
    command += " ILLUMINACLIP:{}:2:30:7:1:keepBothReads".format(ADAPTER_FNAME)
    f.write("\n>adapter/2\n")
    f.write(args.a2)
else:
    BASE += "se.fastq"
    command += " SE {} {}".format(args.f1, BASE)
    command += " ILLUMINACLIP:{}:2:30:7".format(ADAPTER_FNAME)

f.close()

if args.e is not None:
    command += " -r {}".format(str(args.e))

command += " -threads {} MINLEN:{}".format(str(args.t), str(args.m))
# ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>
#             :<simple clip threshold>:<minAdapterLength>:<keepBothReads>
# command += " ILLUMINACLIP:{}:2:30:7:1:keepBothReads".format(ADAPTER_FNAME)

print(command)
if args.rc:
    output_fname = ""
    if is_paired:
        output_fname = "TM_pe_mm_{}.txt".format(str(args.t))
    else:
        output_fname = "TM_se_mm_{}.txt".format(str(args.t))
    os.system("/usr/bin/time -o {} -v {}".format(output_fname, command))
else:
    os.system(command)

