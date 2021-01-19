# EARRINGS

EARRINGS is an efficient and accurate adapter trimmer that entails no a priori adapter sequences for both single- and paired-end NGS reads.

## Information

EARRINGS is a more powerful and capable successor of [PEAT](https://github.com/jhhung/PEAT), which is now deprecated.

Like PEAT, EARRING adapts [skewer](https://github.com/relipmoc/skewer)for single-end read trimming.

## Requirement

- [g++-8](https://gcc.gnu.org/gcc-8/) and cmake [3.10.0](https://cmake.org/download/) or higher to build EARRINGS
- python3.7 or higher for the benchmarking

## Build

```sh
# In root directory of repository
> mkdir build
> cd build
# If default gcc or g++ version != 8.0
# -DCMAKE_C_COMPILER=/usr/bin/gcc-8 or -DCMAKE_CXX_COMPILER=/usr/bin/g++-8 is needed.
> cmake .. -DBUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../stage
> cmake --build . --target install
> ./EARRINGS -h
```

## Execution

There are 3 modes to execute EARRINGS: **build**, **single**, and **paired**.

Build mode generates an index for the source reference sequence (e.g., the entire genome, a chromosome, a collection of panel genes, etc.) of the **single-end** reads. The index is **not** used for trimming paired-end reads.

Single mode and paired mode are used for single-end reads and paired-end reads respectively.


### **Single-End**

### Build reference index

Before conducting single-end adapter trimming, **one has to prebuild the index** once for a specific reference 
which is the source of the target reads.

```sh
# ./EARRINGS build -r [ref_path] -p [index_prefix]
> ./EARRINGS build -r hg38_chr1.fa -p earrings_hg38_chr1
> ./EARRINGS build -r 16S_rRNA_panel.fa -p 16S_rRNA_penel
```

Build parameters

- Required
  - -r [ --ref_path ]                     Path to the reference genome.
  - -p [ --index_prefix ]                 The index prefix for the built index table.
- Optional
  - -h [ --help ]                         Display help message and exit.

### Execute Single-End trimming

EARRINGS first detects adapter then feeds the detected adapter to skewer.

```sh
# ./EARRINGS single -p [index_prefix] --skewer [input_file] [skewer_parameters]
> ./EARRINGS single -p path_to_index --skewer ../test_data/has_adapter_1.fq
```

Single-End mode parameters

- Required
  - -p [ --index_prefix ] The index prefix for pre-built index table.
  - -s [ --skewer ] Skewer flag, options after this would be fed to Skewer, such as input</br> 
  					file name and the **number of thread** used to run the program. These two parameters</br> 
  					will also be used by EARRINGS. 
- Optional
  - -h [ --help ] Display help message and exit.
  - -d [ --seed_len ] The first seed_len bases at 5' portion is viewed as seed when conducting alignment.</br> 
  		EARRINGS allows at most one mismatch in the seed portion if enable_mismatch is set to true.</br> 
  		Reads will be aborted if more than one mismatch is found in the seed portion. If one mismatch</br> 
  		is found outside the seed region, the remainder is reported as a tail. It is recommended to set</br> 
  		the seed_len to 18 for very short reads like miRNA, otherwise, it is recommended to set it to 50. (default: 50)
  - -m [ --max_align ] Maximum number of candidates used in seed finding stage. (default: 0, not limited)
  - -e [ --enable_mismatch ] Enable/disable mismatch when conducting seed finding. (default: true)
  - -f [ --prune_factor ] Prune factor used when assembling adapters using the de-brujin graph. </br>
        Kmer frequency lower than the prune factor will be aborted. (default: 0.03)
  - -F [ --fasta ] Specify input file type as FastA. (Default input file format: FastQ)
  - -a [ --adapter ] Alternative adapter if auto-detect mechanism fails.</br>
        (default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)
  - --sensitive Sensitive mode can be used when the user is sure that the dataset contains adapters.</br>
        Under sensitive mode, we do not restrict the minimum number of occurence of kmers when assembly</br>
        adapters. By default, the minimum number of occurence of kmers must exceed 10.
  - -b [ --bam_input ] Transform reads in a BAM file into a FastA file. Then trim off adapters from</br>
  		the FastA file. The file name of the untrimmed Fasta file is specified by the input file name</br>
  		for Skewer, while the file name of the trimmed Fasta file is also specified by the output file name for Skewer.
  - -u [ --UMI ] Estimate the size of UMI sequences. The result will be printed on the screen.

### **Paired-End**

### Execute Paired-End trimming

```sh
# ./EARRINGS paired -i [input1] -I [input2] -t [thread_num]
> ./EARRINGS paired -i ../test_data/has_adapter_1.fq -I ../test_data/has_adapter_2.fq
```

Paired-end mode parameters

- Required
  - -i [ --input1 ] The PE FastQ input file 1 (.fq)
  - -I [ --input2 ] The PE FastQ input file 2 (.fq)
- Optional
  - -h [ --help ] Display help message and exit.
  - -o [ --output1 ] The PE FastQ output file 1 (.fq) (default: EARRINGS_2.fq)
  - -O [ --output2 ] The PE FastQ output file 2 (.fq) (default: EARRINGS_2.fq)
  - -a [ --adapter1 ] Default adapter 1 when auto-detect fails.</br>
        (default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)
  - -A [ --adapter2 ] Default adapter 2 when auto-detect fails.</br>
        (default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA)
  - -t [ --thread ] The number of threads used to run the program. (default: 1)
  - -m [ --min_length ] Abort the read if the length of the read is less than m. (default: 0)
  - -M [ --rc_thres ] Mismatch threshold applied in reverse complement scan. (default: 0.7)
  - -s [ --ss_thres ] Mismatch threshold applied in gene portion check. (default: 0.9)
  - -S [ --as_thres ] Mismatch threshold applied in adapter portion check. (default: 0.8)
  - -f [ --prune_factor ] Prune factor used when assembling adapters using the de Bruijn graph. kmer</br>
        frequency lower than the prune factor will be aborted. (default: 0.03)
  - --sensitive Sensitive mode can be used when the user is sure that the dataset contains adapters.</br>
        Under sensitive mode, we do not restrict the minimum number of occurence of kmers when assembly adapters.</br>
        By default, the minimum number of occurence of kmers must exceed 10.
  - -F [ --fasta ] Specify input file type as FastA. (default input file format: FastQ)
  - -b [ --bam_input ] Transform reads in a BAM file into two FastA files. Then trim off adapters</br> 
  		from the FastA files. The file names of the untrimmed Fasta files are determined by input1 and</br> 
  		input2, while the file names of the trimmed Fasta files are determined by output1 and output2.

## Run Simulation

The simulation is carried out using a modified version of [pIRS](https://github.com/galaxy001/pirs)
(profile based Illumina pair-end Reads Simulator).
Please build pIRS first before using it:

```sh
> cd simulator/pIRS
> make
```

pIRS simulates Illumina PE reads from a reference genome. In all the benchmarking except real data </br>
benchmarking, we use GRCh38 chr1 as reference genome to simulate reads.

One can run the simulation by:

```sh
# python3 script.py [seed] -r [ref_genome] -n [read_num] -m [insert_mean] -std [insert_td] -ad1 [ad1] -ad2 [ad1]
# To generate reads without adapters:
> python3 script.py 0 -r hg38_chr1.fa -n 10000 -m 150 -std 50 -ad1 "" -ad2 ""
# To generate reads with default adapters:
> python3 script.py 0 -r hg38_chr1.fa -n 10000 -m 150 -std 50
```

## Run Benchmarking

Before running benchmarking, please install all the prerequisites and set the locations of all the </br>
executions to $PATH:

1. AdapterRemoval ver. 2.3.0
2. skewer ver. 0.2.2
3. Cutadapt ver. 2.4
4. AKATrim ver. 1.3.3
5. Trimmomatic ver. 0.39
6. SeqPurge ver. 2019_11
7. fastp ver. 0.20.0 - commit 6ff0ffa
8. atropos ver. 1.1.21 - commit 2b15c77
9. PEAT ver. 1.2.5
10. python ver. 3.6
11. our modified version of pIRS (see Run Simulation part)
12. picard 2.23.0

For speed, memory, performance and adapter benchmarking, please download the first chromosome of GRCh38 from UCSC genome browser.
For running real data benchmarking, please download the following four datasets and their corresponding reference genome(.fa) from ncbi.

1. SRR529095: RNA IP-Seq of Argonaute proteins from Homo sapiens
2. SRR014866: miRNA of C.elegans
3. SRR330569: RNA-Seq of Gonads and Carcasses in D. simulans and D. pseudoobscura
4. SRR5000681: ATAC-Seq of Early Embryo replicate in C.elegans

Also, don't forget to install the aligners:

1. Bowtie2 v2.4.1
2. HISAT2 v2.2.0

After that, modify path-related variables in path.py to the path where the downloaded datasets were placed.
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

## Contact

```sh
Jui-Hung Hung <juihunghung@gmail.com>

Ting-Husan Wang <phoebewangintw@gmail.com>

Cheng-Ching Huang <ken5112840801@gmail.com>
```
