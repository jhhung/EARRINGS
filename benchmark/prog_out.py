import gzip
import shutil
import os

AD = ["AD", "AdapterRemoval_se.truncated", "AdapterRemoval_pe.pair1.truncated", "AdapterRemoval_pe.pair2.truncated"]
SK = ["SK", "skewer_se-trimmed.fastq", "skewer_pe-trimmed-pair1.fastq", "skewer_pe-trimmed-pair2.fastq"]
CA = ["CA", "cutadapt_se.fastq", "cutadapt_pe_1.fastq", "cutadapt_pe_2.fastq"]
FP = ["FP", "fastp_se.fastq", "fastp_pe1.fastq", "fastp_pe2.fastq"]
FP_SORT = ["FP", "fastp_se_sorted.fastq", "fastp_pe1_sorted.fastq", "fastp_pe2_sorted.fastq"]
AP = ["AP", "atropos_se.fastq", "atropos_pe1.fastq", "atropos_pe2.fastq"]
AP_SORT = ["AP", "atropos_se_sorted.fastq", "atropos_pe1_sorted.fastq", "atropos_pe2_sorted.fastq"]
AKATrim = ["AKA", "None", ".akatrim.", ".akatrim."]
EA = ["EA", "earrings-trimmed.fastq", "earrings_pe1.fastq", "earrings_pe2.fastq"]
PE = ["PE", "PEAT_se.fastq", "PEAT_pe1.fastq.gz", "PEAT_pe2.fastq.gz"]
PE_UNZIP = ["PE", "PEAT_se.fastq", "PEAT_pe1.fastq", "PEAT_pe2.fastq"]
# TR = ["TR", "trimmomatic_se.fastq", "trimmomatic_pe_1P.fastq", "trimmomatic_pe_2P.fastq", "trimmomatic_pe_1U.fastq", "trimmomatic_pe_2U.fastq"]
TR = ["TR", "trimmomatic_se.fastq", "trimmomatic_pe_1.fastq", "trimmomatic_pe_2.fastq"]
TR_SORT = ["TR", "trimmomatic_se.fastq", "trimmomatic_pe_1_sorted.fastq", "trimmomatic_pe_2_sorted.fastq"]
SP = ["SP", "None", "seqpurge_pe1.fastq.gz", "seqpurge_pe2.fastq.gz"]
SP_UNZIP = ["SP", "None", "seqpurge_pe1.fastq", "seqpurge_pe2.fastq"]
SP_SORT = ["SP", "None", "seqpurge_pe1_sorted.fastq", "seqpurge_pe2_sorted.fastq"]

SE_PROG_OUT = [AD, SK, CA, FP, AP, EA, PE, TR]
SE_PROG_SORT_OUT = [AD, SK, CA, FP_SORT, AP_SORT, EA, PE, TR_SORT]

# don't forget to append AKA & SP
PE_PROG_OUT = [AD, SK, CA, FP, AP, EA, PE_UNZIP, TR, SP_UNZIP] 
PE_PROG_SORT_OUT = [AD, SK, CA, FP_SORT, AP_SORT, EA, PE, TR_SORT, SP_SORT]

def concat_files(filenames, output_f):
    for f in filenames:
        if not os.path.isfile(f):
            return
            
    with open(output_f, 'w') as outfile:
        for fname in filenames:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

def concat_TR_PE(directory):
    filenames = [os.path.join(directory, 'trimmomatic_pe_1P.fastq'), os.path.join(directory, 'trimmomatic_pe_1U.fastq')]
    concat_files(filenames, os.path.join(directory, TR[2]))
    filenames = [os.path.join(directory, 'trimmomatic_pe_2P.fastq'), os.path.join(directory, 'trimmomatic_pe_2U.fastq')]
    concat_files(filenames, os.path.join(directory, TR[3]))
    
def append_AKA(in1, in2, prog_list, suffix="fq"):
    new_prog_list = list(prog_list)
    AKAtmp = [AKATrim[0], AKATrim[1]]
    AKAtmp.append(in1[:-len(suffix)-1] + AKATrim[2] + suffix)
    AKAtmp.append(in2[:-len(suffix)-1] + AKATrim[3] + suffix)
    new_prog_list.append(AKAtmp)

    return new_prog_list

def unzip(in_f, out_f):
    with gzip.open(in_f, 'rb') as f_in:
        with open(out_f, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

# rename the paired-end output files of PEAT
# since the output is *.fastq instead of *.fastq.gz
def rename_PE(f_before, f_after):
    os.system("mv {} {}".format(f_before, f_after))