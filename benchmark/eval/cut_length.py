from .fastq import Fastq
import os
# this is created for TM, since the min output option of TR does not work properly.

tmp_fname = ["tmp1", "tmp2"]
def cut_length(input_fs, min_len):
    if len(input_fs) > 1:
        cut_length_2f(input_fs, min_len)
    else:
        cut_length_1f(input_fs[0], min_len)

    
def cut_length_2f(input_fs, min_len):
    ifs1 = open(input_fs[0], "r")
    ifs2 = open(input_fs[1], "r")
    fq1 = Fastq()
    fq2 = Fastq()
    ofs1 = open(tmp_fname[0], "w")
    ofs2 = open(tmp_fname[1], "w")
    while not fq1.get_obj(ifs1) and not fq2.get_obj(ifs2):
        if len(fq1.seq) < min_len or len(fq2.seq) < min_len:
            continue
        fq1.write(ofs1)
        fq2.write(ofs2)
    ofs1.close()
    ofs2.close()

    for i in range(2):
        os.remove(input_fs[i])
        os.rename(tmp_fname[i], input_fs[i])
    

def cut_length_1f(input_f, min_len):
    ifs = open(input_f, "r")
    fq = Fastq()
    ofs = open(tmp_fname[0], "w")
    while not fq.get_obj(ifs):
        if len(fq.seq) < min_len:
            continue
        fq.write(ofs)
    ofs.close()

    os.remove(input_f)
    os.rename(tmp_fname[0], input_f)