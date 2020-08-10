import os
from parser import parse_args

NAME = "AKA"
EXE = "AKATrim"
args = parse_args()

# disable other trimming/filtering method
command = "{} --no-tail-n-trim --max-n=-1 --no-quality-trim ".format(EXE)
command += "--length-range {}:500 ".format(args.m)

if args.f2 is not None:
    command += "--read1 {} --read2 {} --procs {} --adapter1 {} --adapter2 {}"\
                        .format(str(args.f1), str(args.f2), str(args.t), args.a1, args.a2)
    print(command)
    if args.rc:
        output_fname = "AKA_pe_mm_{}.txt".format(str(args.t))
        os.system("/usr/bin/time -o {} -v {}".format(output_fname, command))
    else:
        os.system(command)