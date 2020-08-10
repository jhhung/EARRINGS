#!/usr/bin/env python3
import os
from parser import parse_args

NAME = "SK"
EXE = "skewer"
BASE = "skewer_"
is_paired = False
args = parse_args()
command = EXE

if args.f2 is not None:
    is_paired = True
    BASE += "pe"
    command += " -m pe {} {}".format(args.f1, args.f2)
else:
    BASE += "se"
    command += " {}".format(args.f1)

if args.a1 is not None:
    command += " -x {}".format(args.a1)

if args.a2 is not None:
    command += " -y {}".format(args.a2)

if args.e is not None:
    command += " -r {}".format(str(args.e))

command += " -t {} -l {} -o {}".format(str(args.t), str(args.m), BASE)

print(command)
if args.rc:
    output_fname = ""
    if is_paired:
        output_fname = "SK_pe_mm_{}.txt".format(str(args.t))
    else:
        output_fname = "SK_se_mm_{}.txt".format(str(args.t))
    os.system("/usr/bin/time -o {} -v {}".format(output_fname, command))
else:
    os.system(command)

