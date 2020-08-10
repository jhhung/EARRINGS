import os

### modify the following path to the path where you put your data. ###

# reference path
# HG38_REF = ...
# HG38_CHR1_REF = ...

# earrings indexing path
# HG38_EARRINGS_IDX = ...
# HG38_CHR1_EARRINGS_IDX = ...

# aligner exe
# BOWTIE2_EXE = ...
# HISAT2_EXE = ...

### end of modifying path ###

PATH_PWD = os.path.dirname(os.path.realpath(__file__))
EARRINGS_EXE = os.path.join(PATH_PWD, "..", "build", "EARRINGS")

# prog path
RUN_TRIMMER_DIR = os.path.join(PATH_PWD, "run_trimmers")
RUN_SIM_EXE = os.path.join(PATH_PWD, "..", "simulator", "pIRS", "script.py")
SORT_PROG = os.path.join(PATH_PWD, "eval", "order.py")
SE_EVAL = os.path.join(PATH_PWD, "eval", "se_evaluation.py")
PE_EVAL = os.path.join(PATH_PWD, "eval", "pe_evaluation.py")
