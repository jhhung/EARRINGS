# EARRINGS

EARRINGS is an efficient and accurate adapter trimmer that entails no a priori adapter sequences for both single- and paired-end NGS reads.

## Information

EARRINGS is a more powerful and capable successor of [PEAT](https://github.com/jhhung/PEAT), which is now deprecated.

Like PEAT, EARRING adapts [skewer](https://github.com/relipmoc/skewer) for single-end read trimming.

## Requirement

- GUN [g++-8](https://gcc.gnu.org/gcc-8/) or higher (Gtest works with GUN [g++-8](https://gcc.gnu.org/gcc-8/) only )
- cmake [3.10.0](https://cmake.org/download/) or higher to build EARRINGS
- python [3.7](https://www.python.org/downloads/) or higher as well as [numpy](https://numpy.org/) and [pandas](https://pandas.pydata.org/) packages for the benchmarking

## Build

```sh
# In root directory of repository
> mkdir build
> cd build
> cmake .. -DBUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../stage
> cmake --build . --target install
> ./EARRINGS -h
```

### **Q&A**

Q: Failed while compiling and get the error message about ```fatal: unable to access 'https://github.com/PhoebeWangintw/hunter.git/': Could not resolve host: github.com;```.\
A: Please make sure your connection to ```github.com`` is clear, especially for the user who lives inside the GFW of China.

Q: How do I change my default GUN GCC compiler?\
A: To change default gcc and g++ version with 8.0 as example:\
```
# This environment setting works in the current session only.
export CC=/Path/Locate/To/The/gcc-8
export CXX=/Path/Locate/To/The/g++-8
```

## Execution

There are 3 modes to execute EARRINGS: **build**, **single**, and **paired**.

Build mode generates an index for the source reference sequence (e.g., the entire genome, a chromosome, a collection of panel genes, [GreenGenes](http://greengenes.secondgenome.com/) for meta-genomics etc.) of the **single-end** reads. The index is **not** used for trimming paired-end reads.

Single mode and paired mode are used for single-end reads and paired-end reads respectively. Both of these two modes can auto-detect file types. ([.fa](https://en.wikipedia.org/wiki/FASTA_format)/[.fq](https://en.wikipedia.org/wiki/FASTQ_format), and their [.gz](https://en.wikipedia.org/wiki/FASTQ_format#General_compressors) or [.bam](https://en.wikipedia.org/wiki/Binary_Alignment_Map)/[.ubam](http://129.130.90.13/ion-docs/GUID-C202F9D0-386F-412D-97F9-E4DB77F1BB6E.html))

### **Build**

Before conducting single-end adapter trimming, **one has to prebuild the index** once for a specific reference which is the source of the target reads.

```sh
# ./EARRINGS build -r [ref_path] -p [index_prefix]
> ./EARRINGS build -r hg38_chr1.fa -p earrings_hg38_chr1
> ./EARRINGS build -r 16S_rRNA_panel.fa -p 16S_rRNA_panel
```

Build mode parameters

- Required
  - -r [ --ref_path ] arg</br>
  Path to the reference sequence, which is the source of reads.
  - -p [ --index_prefix ] arg</br>
  An user-defined index prefix for index table.
- Optional
  - -h [ --help ]</br>
  Display help message and exit.

### **Single-End**

In single-end mode, EARRINGS first detects adapter then feeds the detected adapter to skewer.

```sh
# ./EARRINGS single -p [index_prefix] -1 [input_file]
> ./EARRINGS single -p path_to_index -1 ../test_data/has_adapter_1.fq
> ./EARRINGS single -p path_to_index -1 ../test_data/has_adapter_1.fq.gz
```

Single-End mode parameters

- Required
  - -p [ --index_prefix ] arg</br>
  The index prefix for pre-built index table.
  - -1 [ --input1 ] arg</br>
  The file path of Single-End reads.
- Optional
  - Utils
    - -h [ --help ]</br>
    Display help message and exit.
    - -t [ --thread ] arg (=1)</br>
    The number of threads used to run the program.
  - Input / Output
    - -o [ --output ] arg (=trimmed_se)</br>
    The file prefix of Single-End FastQ output.
  - Extract seeds / Alignment
    - -d [ --seed_len ] arg (=50)</br>
    The first ***--seed_len*** bases are seen as seed and allows 1 mismatch at most, or do not allow any mismatch if ***--no_mismatch*** is set. The sequence follows first mismatch out of the seed portion will be reported as a tail.</br>
    It's recommended to set this to 18 for very short reads like **miRNA**, otherwise, it is recommended to set to 50.
    - -e [ --no_mismatch ]</br>
    By default, EARRINGS can tolerate 1 error base at most if be set as true, this flag can disable this mismatch toleration mechenism.
    - -M [ --max_align ] arg (=0)</br>
    Maximum number of candidates used in seed finding stage, 0 means unlimited.
  - Assemble adapter
    - -f [ --prune_factor ] arg (=0.03)</br>
    Prune factor used when assembling adapters using the de-brujin graph. Kmer frequency lower than this value will be skipped.
    - --sensitive</br>
    By default, minimum number of kmers must exceed 10 during assembly adapters. However, if user have confidence that the dataset contains adapters, sensitive mode would be more suitable.</br>
    Under sensitive mode, minimum number of kmers (***--prune_factor***) would not be restricted.
  - Trimming
    - -m [ --min_length ] arg (=0)</br>
    Skip the read if the length of the read is less than ***--min_length*** after trimming.
  - Adapter setting
    - -a [ --adapter1 ] arg (=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)</br>
    Alternative adapter if auto-detect mechanism fails.
    - -u [ --UMI ]</br>
    Estimate the size of UMI sequences, results will be printed to console by default.

### **Paired-End**

```sh
# ./EARRINGS paired -1 [input1] -2 [input2]
> ./EARRINGS paired -1 ../test_data/has_adapter_1.fq -2 ../test_data/has_adapter_2.fq
> ./EARRINGS paired -1 ../test_data/has_adapter_1.fq.gz -2 ../test_data/has_adapter_2.fq.gz
```

Paired-end mode parameters

- Required
  - -1 [ --input1 ] arg</br>
    The Paired-end reads input file 1.
  - -2 [ --input2 ] arg</br>
    The Paired-end reads input file 2, note that if input file is .bam/.ubam, this parameter has no function.
- Optional
  - Utils
    - -h [ --help ]</br>
    Display help message and exit.
    - -t [ --thread ] arg (=1)</br>
    The number of threads used to run the program.
  - Input / Output
    - -o [ --output ] arg (=trimmed_pe)</br>
    The Paired-End FastQ output file prefix.
  - Assemble adapter
    - -f [ --prune_factor ] arg (=0.03)</br>
    Prune factor used when assembling adapters using the de Bruijn graph. Kmer frequency lower than the prune factor will be skipped.
    - --sensitive</br>
    By default, minimum number of kmers must exceed 10 during assembly adapters. However, if user have confidence that the dataset contains adapters, sensitive mode would be more suitable.</br>
    Under sensitive mode, minimum number of kmers (***--prune_factor***) would not be restricted.
  - Trimming
    - -m [ --min_length ] arg (=0)</br>
    Skip the read if the length of the read is less than this value.
    - -M [ --rc_thres ] arg (=0.7)</br>
    Mismatch threshold applied in reverse complement scan.
    - -s [ --ss_thres ] arg (=0.9)</br>
    Mismatch threshold applied in gene portion check.
    - -S [ --as_thres ] arg (=0.8)</br>
    Mismatch threshold applied in adapter portion check.
  - Adapter setting
    - -a [ --adapter1 ] arg (=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)</br>
    Alternative adapter 1 if auto-detect mechanism fails.
    - -A [ --adapter2 ] arg (=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA)</br>
    Alternative adapter 2 if auto-detect mechanism fails.

## Run Simulation

The simulation is carried out using a modified version of [pIRS](https://github.com/galaxy001/pirs)
(profile based Illumina pair-end Reads Simulator).</br>
Please build pIRS first before using it:

```sh
> cd simulator/pIRS
> make
```

pIRS simulates Illumina PE reads from a reference genome. In all the benchmarking except real data benchmarking, we use GRCh38 chr1 as reference genome to simulate reads.

One can run the simulation by:

```sh
# python3 script.py [seed] -r [ref_genome] -n [read_num] -m [insert_mean] -std [insert_td] -ad1 [ad1] -ad2 [ad1]
# To generate reads without adapters:
> python3 script.py 0 -r hg38_chr1.fa -n 10000 -m 150 -std 50 -ad1 "" -ad2 ""
# To generate reads with default adapters:
> python3 script.py 0 -r hg38_chr1.fa -n 10000 -m 150 -std 50
```

## Run Benchmarking

Before running benchmarking, please install all the prerequisites and set the locations of all the executions to $PATH:

1. [AdapterRemoval ver. 2.3.0](https://github.com/MikkelSchubert/adapterremoval/tree/v2.3.0)
2. [skewer ver. 0.2.2](https://github.com/relipmoc/skewer/tree/0.2.2)
3. [Cutadapt ver. 2.4](https://github.com/marcelm/cutadapt/tree/v2.4)
4. AKATrim ver. 1.3.3 (no source link now)
5. [Trimmomatic ver. 0.39](http://www.usadellab.org/cms/?page=trimmomatic)
6. [SeqPurge ver. 2019_11](https://github.com/imgag/ngs-bits/tree/2019_11)
7. [fastp ver. 0.20.0 - commit 6ff0ffa](https://github.com/OpenGene/fastp/tree/6ff0ffa5b4ee9968d4369207aca6f31804194533)
8. [atropos ver. 1.1.21 - commit 2b15c77](https://github.com/jdidion/atropos/tree/2b15c778f0ccf1d0fb753e4334fa6dc0048a9ee6)
9. [PEAT ver. 1.2.5](https://github.com/jhhung/PEAT/tree/v1.2.5)
10. [python ver. 3.7 or higher](https://www.python.org/downloads/)
11. [picard 2.23.0](https://github.com/broadinstitute/picard/tree/2.23.0)
12. our modified version of pIRS (see Run Simulation part)

For speed, memory, performance and adapter benchmarking, please download the first chromosome of GRCh38 from UCSC genome browser.
For running real data benchmarking, please download the following four datasets and their corresponding reference genome(.fa) from ncbi.

1. SRR529095: RNA IP-Seq of Argonaute proteins from Homo sapiens
2. SRR014866: miRNA of C.elegans
3. SRR330569: RNA-Seq of Gonads and Carcasses in D. simulans and D. pseudoobscura
4. SRR5000681: ATAC-Seq of Early Embryo replicate in C.elegans

Also, don't forget to install the aligners:

1. Bowtie2 v2.4.1
2. HISAT2 v2.2.0

After that, modify path-related variables in path&#46;py to the path where the downloaded datasets were placed.
Then prebuild indices for the aligners and EARRINGS by running:

```sh
> python3 build_index.py
```

Run real data benchmarking:

```sh
> python3 benchmark_real_data.py
```

Run speed benchmarking:

```sh
> python3 benchmark_speed.py
```

Run memory benchmarking:

```sh
> python3 benchmark_mem_usage.py
```

Run performance benchmarking:

```sh
> python3 benchmark_performance.py
```

Run adapter benchmarking:

```sh
> python3 benchmark_adapter.py
```

## Reference

1. Li, Y.-L., Weng, J.-C., Hsiao, C.-C., Chou, M.-T., Tseng, C.-W., & Hung, J.-H. (2015). PEAT: an intelligent and efficient paired-end sequencing adapter trimming algorithm. BMC Bioinformatics, 16(Suppl 1), S2. doi:10.1186/1471-2105-16-S1-S2</br>
[[paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-16-S1-S2)]
[[github](https://github.com/jhhung/PEAT)]
2. Jiang, H., Lei, R., Ding, S.W. and Zhu, S. (2014) Skewer: a fast and accurate adapter trimmer for next-generation sequencing paired-end reads. BMC Bioinformatics, 15, 182.</br>
[[paper](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-182)]
[[github](https://github.com/relipmoc/skewer)]
3. Hu, X., et al. (2012) pIRS: Profile-based Illumina pair-end reads simulator, Bioinformatics, 28, 1533-1535.</br>
[[paper](https://academic.oup.com/bioinformatics/article/28/11/1533/267409)]
[[github](https://github.com/galaxy001/pirs)]
4. Chou, M.-T., Han, B. W., Hsiao, C.-P., Zamore, P. D., Weng, Z., and Hung, J.-H. (2015). Tailor: a computational framework for detecting non-templated tailing of small silencing RNAs. Nucleic Acids Res. 43, e109.</br>
[[paper](https://academic.oup.com/nar/article/43/17/e109/2414301)]
[[github](https://github.com/jhhung/Tailor)]

## Citation

Wang et al. EARRINGS: An Efficient and Accurate Adapter Trimmer Entails No a Priori Adapter Sequences. Bioinformatics. Accepted.</br>
[[paper](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btab025/6103563)]
[[github](https://github.com/jhhung/EARRINGS)]

## Contact

```sh
Jui-Hung Hung <juihunghung@gmail.com>

Ting-Hsuan Wang <phoebewangintw@gmail.com>

Cheng-Ching Huang <ken5112840801@gmail.com>
```
