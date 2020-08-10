#!/usr/bin/env python3
import os
from parser import parse_args

args = parse_args()
NAME = "AP"
EXE = "atropos"

command = EXE
is_paired = False

if args.da:
    command += " detect "
else:
    command += " trim "
    if args.a1 is not None:
        command += " -a {}".format(args.a1)

    if args.a2 is not None:
        command += " -A {}".format(args.a2)

if args.af is not None:
    command += "--{}".format(args.af)

if args.f2 is not None:
    is_paired = True
    command += " -pe1 {} -pe2 {} -o atropos_pe1.fastq -p atropos_pe2.fastq".format(args.f1, args.f2)
else:
    command += " -se {} -o atropos_se.fastq".format(args.f1)

args.t = 2 if args.t == 1 else args.t
command += " -m {} -T {}".format(str(args.m), str(args.t))

'''
if is_paired:
    command += " --aligner insert"
'''
print(command)
if args.rc:
    output_fname = ""
    if is_paired:
        output_fname = "AP_pe_mm_{}.txt".format(str(args.t))
    else:
        output_fname = "AP_se_mm_{}.txt".format(str(args.t))
    os.system("/usr/bin/time -o {} -v {}".format(output_fname, command))
else:
    os.system(command)
