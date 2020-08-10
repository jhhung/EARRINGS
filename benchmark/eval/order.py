import numpy as np
import sys
import operator
from .fastq import Fastq

def ordering(input_fname, output_fname):
    input_f = open(input_fname, "r")
    fastqs = []
    fq = Fastq()

    while not fq.get_obj(input_f):
        fastqs.append(Fastq(fq))
    fastqs.sort()

    output_f = open(output_fname, "w")

    for f in fastqs:
        output_f.write(f.name + "\n")
        output_f.write(f.seq + "\n")
        output_f.write(f.plus + "\n")
        output_f.write(f.qual + "\n")

    output_f.close()

if __name__ == "__main__":
    ordering(sys.argv[1], sys.argv[2])
