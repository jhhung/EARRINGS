import os

### modify the following path to the path where you put your real data. eg. SRR529095, hg38.fa ###
REAL_DATA_DIR = os.path.join(os.sep, "real_data")

# dir that store reference
HG38_DIR = os.path.join(REAL_DATA_DIR, "ref", "hg38")
HG38_CHR1_DIR = os.path.join(REAL_DATA_DIR, "ref", "hg38_chr1")
C_ELEGANS_DIR = os.path.join(REAL_DATA_DIR, "ref", "c_elegans")
D_SIMULANS_DIR = os.path.join(REAL_DATA_DIR, "ref", "drosophila")

# reference path
HG38_REF = os.path.join(HG38_DIR, "hg38.fa")
HG38_CHR1_REF = os.path.join(HG38_CHR1_DIR, "chr1.fa")
C_ELEGANS_REF = os.path.join(C_ELEGANS_DIR, "c_elegans.fa")
D_SIMULANS_REF = os.path.join(D_SIMULANS_DIR, "d_simu.fa")

# earrings indexing path
HG38_CHR1_EARRINGS_IDX = os.path.join(HG38_CHR1_DIR, "earrings_hg38_chr1.")
HG38_EARRINGS_IDX = os.path.join(HG38_DIR, "earrings_hg38.")
C_ELEGAN_EARRINGS_IDX = os.path.join(C_ELEGANS_DIR, "earrings_c_elegans.")

# dir where aligner locates
BOWTIE2_DIR = os.path.join(REAL_DATA_DIR, "..", "bowtie2")
HISAT2_DIR = os.path.join(REAL_DATA_DIR, "..", "hisat2")

# aligner exe
BOWTIE2_EXE = os.path.join(BOWTIE2_DIR, "bowtie2")
BOWTIE2_BUILD_EXE = os.path.join(BOWTIE2_DIR, "bowtie2-build")
HISAT2_EXE = os.path.join(HISAT2_DIR, "hisat2")
HISAT2_BUILD_EXE = os.path.join(HISAT2_DIR, "hisat2-build")

# aligner indexing path
HG38_HISAT2_IDX = os.path.join(HG38_DIR, "hisat2_hg38")
C_ELEGAN_BOWTIE2_IDX = os.path.join(C_ELEGANS_DIR, "bowtie2_c_elegans")
D_SIMULANS_HISAT2_IDX = os.path.join(D_SIMULANS_DIR, "hisat2_drosophila")

# dir that store real data
SRR529095_DIR = os.path.join(REAL_DATA_DIR, "SRR529095")
SRR014966_DIR = os.path.join(REAL_DATA_DIR, "SRR014966")
SRR330569_DIR = os.path.join(REAL_DATA_DIR, "SRR330569")
SRR5000681_DIR = os.path.join(REAL_DATA_DIR, "SRR5000681")

# real data path
SRR529095_PATH = os.path.join(SRR529095_DIR, "SRR529095.1.fastq")
SRR014966_PATH = os.path.join(SRR014966_DIR, "SRR014966.2.fastq")
SRR330569_PATH = [os.path.join(SRR330569_DIR, "SRR330569.3_1.fastq"), 
                  os.path.join(SRR330569_DIR, "SRR330569.3_2.fastq")]
SRR5000681_PATH = [os.path.join(SRR5000681_DIR, "SRR5000681.1_1.fastq"), 
                   os.path.join(SRR5000681_DIR, "SRR5000681.1_2.fastq")]

### end of modifying path ###

PATH_PWD = os.path.dirname(os.path.realpath(__file__))
EARRINGS_EXE = os.path.join(PATH_PWD, "..", "build", "EARRINGS")

# prog path
RUN_TRIMMER_DIR = os.path.join(PATH_PWD, "run_trimmers")
RUN_SIM_EXE = os.path.join(PATH_PWD, "..", "simulator", "pIRS", "script.py")
SORT_PROG = os.path.join(PATH_PWD, "eval", "order.py")
SE_EVAL = os.path.join(PATH_PWD, "eval", "se_evaluation.py")
PE_EVAL = os.path.join(PATH_PWD, "eval", "pe_evaluation.py")

def build_correspoinding_dir():
    os.makedirs(HG38_DIR, exist_ok=True)
    os.makedirs(HG38_CHR1_DIR, exist_ok=True)
    os.makedirs(C_ELEGANS_DIR, exist_ok=True)
    os.makedirs(D_SIMULANS_DIR, exist_ok=True)

    os.makedirs(SRR529095_DIR, exist_ok=True)
    os.makedirs(SRR014966_DIR, exist_ok=True)
    os.makedirs(SRR330569_DIR, exist_ok=True)
    os.makedirs(SRR5000681_DIR, exist_ok=True)

if __name__ == "__main__":
    build_correspoinding_dir()
