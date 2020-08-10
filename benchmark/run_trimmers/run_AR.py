#!/usr/bin/env python3
import os
from parser import parse_args

NAME = "AR"
EXE = "AdapterRemoval"
BASE = "AdapterRemoval_"

args = parse_args()
command = "{} --file1 {}".format(EXE, args.f1)
is_paired = False
if args.f2 is not None:
    is_paired = True
    BASE += "pe"
    command += " --file2 {}".format(args.f2)
else:
    BASE += "se"

if args.a1 is not None:
    command += " --adapter1 {}".format(args.a1)

if args.a2 is not None:
    command += " --adapter2 {}".format(args.a2)

if args.e is not None:
    command += " --mm {}".format(str(args.e))

command += " --threads {} --minlength {} --basename {}".format(str(args.t), str(args.m), BASE)

print(command)
if args.rc:
    output_fname = ""
    if is_paired:
        output_fname = "AR_pe_mm_{}.txt".format(str(args.t))
    else:
        output_fname = "AR_se_mm_{}.txt".format(str(args.t))
    os.system("/usr/bin/time -o {} -v {}".format(output_fname, command))
else:
    os.system(command)

