#!/usr/bin/env python3
import os
from parser import parse_args

NAME = "EA"
command = ""
args = parse_args()
is_paired = False

PWD = os.path.dirname(os.path.realpath(__file__))
EXE = os.path.join(PWD, "..", "..", "build", "EARRINGS")

if args.f2 is not None:
    is_paired = True
    command += "{} paired -1 {} -2 {} -o earrings_pe -t {} -m {}".format(EXE, args.f1, args.f2, str(args.t), str(args.m))
    command += " -M {} -s {} -S {} -f {}".format(str(args.mr), str(args.sc), str(args.ac), str(args.pf))
    if args.se:
        command += " --sensitive "
else:
    if args.se:
        command += "{} single -p {} -d {} -f {} --sensitive -1 {} -t {} -m {} -o earrings_se".format(EXE, args.idx, str(args.mp), str(args.pf), args.f1, str(args.t), str(args.m))
    else:
        command += "{} single -p {} -d {} -f {} -1 {} -t {} -m {} -o earrings_se".format(EXE, args.idx, str(args.mp), str(args.pf), args.f1, str(args.t), str(args.m))

    if args.e is not None:
        command += " -r {}".format(str(args.e))


print(command)
if args.rc:
    output_fname = ""
    if is_paired:
        output_fname = "EA_pe_mm_{}.txt".format(str(args.t))
    else:
        output_fname = "EA_se_mm_{}.txt".format(str(args.t))
    os.system("/usr/bin/time -o {} -v {}".format(output_fname, command))
else:
    os.system(command)
