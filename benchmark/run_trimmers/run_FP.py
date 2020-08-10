#!/usr/bin/env python3
import os
from parser import parse_args

args = parse_args()
NAME = "FP"
EXE = "fastp"
command = "{} -i {}".format(EXE, args.f1)
is_paired = False
if args.f2 is not None:
    is_paired = True
    command += " -I {} -o fastp_pe1.fastq -O fastp_pe2.fastq".format(args.f2)
else:
    command += " -o ./fastp_se.fastq "

if not args.da:
    if args.a1 is not None:
        command += " --adapter_sequence {}".format(args.a1)
    if args.a2 is not None:
        command += " --adapter_sequence_r2 {}".format(args.a2)
else:
    command += " --detect_adapter_for_pe"
if args.px:
    command += " -x "

command += " -w {} -l {} -Q".format(str(args.t), str(args.m))

print(command)
if args.rc:
    output_fname = ""
    if is_paired:
        output_fname = "FP_pe_mm_{}.txt".format(str(args.t))
    else:
        output_fname = "FP_se_mm_{}.txt".format(str(args.t))
    os.system("/usr/bin/time -o {} -v {}".format(output_fname, command))
else:
    os.system(command)

