#!/usr/bin/env python3
import os
from parser import parse_args

# SeqPurge (2019_11)
# only for PE
NAME = "SP"
EXE = "SeqPurge"

args = parse_args()
command = "{} -in1 {}".format(EXE, args.f1)

if args.f2 is not None:
    command += " -in2 {} -out1 seqpurge_pe1.fastq.gz -out2 seqpurge_pe2.fastq.gz -summary seqpurge_summary.txt".format(args.f2)

    if args.a1 is not None:
        command += " -a1 {}".format(args.a1)

    if args.a2 is not None:
        command += " -a2 {}".format(args.a2)

    command += " -threads {} -min_len {}".format(str(args.t), str(args.m))

    print(command)
    if args.rc:
        output_fname = "SP_pe_mm_{}.txt".format(str(args.t))
        os.system("/usr/bin/time -o {} -v {}".format(output_fname, command))
    else:
        os.system(command)

