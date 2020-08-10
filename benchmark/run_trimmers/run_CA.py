#!/usr/bin/env python3
import os
from parser import parse_args

# EXEC = "~/.local/bin/cutadapt"
EXE = "cutadapt"
NAME = "CA"
args = parse_args()
command = EXE
files = args.f1
is_paired = False

if args.f2 is not None:
    is_paired = True
    command += " -o cutadapt_pe_1.fastq -p cutadapt_pe_2.fastq"
    files += " " + args.f2
else:
    command += " -o ./cutadapt_se.fastq"

if args.a1 is not None:
    command += " -a {}".format(args.a1)

if args.a2 is not None:
    command += " -A {}".format(args.a2)

if args.e is not None:
    command += " -e {}".format(str(args.e))

command += " -j {} -m {} {}".format(str(args.t), str(args.m), files)

print(command)
if args.rc:
    output_fname = ""
    if is_paired:
        output_fname = "CA_pe_mm_{}.txt".format(str(args.t))
    else:
        output_fname = "CA_se_mm_{}.txt".format(str(args.t))
    os.system("/usr/bin/time -o {} -v {}".format(output_fname, command))
else:
    os.system(command)

