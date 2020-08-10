#!/usr/bin/env python3
import os
from parser import parse_args

NAME = "PE"
command = None
args = parse_args()
is_paired = False
PWD = os.path.dirname(os.path.realpath(__file__))
EXE = "PEAT"

if args.f2 is not None:
    is_paired = True
    command = "{} paired -1 {} -2 {} --output_1 PEAT_pe1.fastq --output_2 PEAT_pe2.fastq -n {} -l {}".format(EXE, args.f1, args.f2, str(args.t), str(args.m))
else:
    command = "{} single -i {} -a {} -o PEAT_se.fastq -n {} -q SANGER".format(EXE, args.f1, args.a1, str(args.t))
    
print(command)
if args.rc:
    output_fname = ""
    if is_paired:
        output_fname = "PE_pe_mm_{}.txt".format(str(args.t))
    else:
        output_fname = "PE_se_mm_{}.txt".format(str(args.t))
    os.system("/usr/bin/time -o {} -v {}".format(output_fname, command))
else:
    os.system(command)

