# EARRINGS
EARRINGS is an efficient and accurate adapter trimmer that entails no a priori adapter sequences.

## Requirement
- g++8 and cmake 3.10.0 or higher to build EARRINGS
- python3.7 or higher to run benchmarking

## Build
```sh
> mkdir build
> cd build
> cmake .. -DBUILD_TESTS=OFF -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../stage
> cmake --build . --target install
> ./EARRINGS -h
```

## Execution

### SE
#### SE build
Before conducting SE adapter trimming, **one has to prebuild index** once for a specific reference genome.
#### run SE build
```sh
# ./EARRINGS build -r [ref_path] -p [index_prefix]
> ./EARRINGS build -r ref/hg38_chr1.fa -p earrings_hg38_chr1.
```
#### SE build parameters
- -h [ --help ]                         Display help message and exit.<br />
- -r [ --ref_path ]                     Path to the reference genome.<br />
- -p [ --index_prefix ]                 The index prefix for the built index table.<br />
#### run SE
EARRINGS first detects adapter then feeds the detected adapter to skewer.
```sh
# ./EARRINGS single -p [index_prefix] --skewer [input_file] [skewer_parameters]
> ./EARRINGS single -p path_to_index. --skewer ../test_data/has_adapter_1.fq
```
#### SE parameters
- -h [ --help ]                         Display help message and exit.<br />
- -p [ --index_prefix ]                 The index prefix for prebuilt index table..<br />
- -d [ --seed_len ]                     Seed length used when aligning reads. For very short reads like miRNA, 
                                        it is recommended to set seed_len to 18. (default: 50)<br />
- -m [ --max_align ]                    Control the maximum number of alignment
                                        to abort the reads. (default: 0, not limited)<br />
- -e [ --enable_mismatch ]              Enable/disable mismatch when doing seed
                                        finding. (default: true)<br />
- -f [ --prune_factor ]                 Prune factor used when assembling adapters using the de Bruijn graph.
                                        kmer occurence lower than the prune factor will be aborted.(default: 0.03)<br />
- -F [ --fasta ]                        Specify input file type as FastA. (Default input file format: FastQ)<br />
- -a [ --adapter1 ]                     Default adapter used when auto-detect fails. (default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)<br />
- --sensitive                           Sensitive mode can be used when the user is sure that the dataset 
                                        contains adapters. Under sensitive mode, we do not restrict the 
                                        minimum number of kmers when assembly adapters. By default, the 
                                        minimum number of kmers must exceed 10.<br/>
- -b [ --bam_input ] arg                Detect and trim off adapters from a BAM file.<br/>
- -u [ --UMI ]                          Estimate the size of UMI sequences.<br/>
- -s [ --skewer ]                       Skewer flag, which follows by Skewer's program options.<br />

### PE
#### run PE
```sh
# ./EARRINGS paired -i [input1] -I [input2] -t [thread_num]
> ./EARRINGS paired -i ../test_data/has_adapter_1.fq -I ../test_data/has_adapter_2.fq
```
#### PE parameters
- -h [ --help ]                         Display help message and exit.<br />
- -i [ --input1 ]                       The PE FastQ input file 1 (.fq)<br />
- -I [ --input2 ]                       The PE FastQ input file 2 (.fq)<br />
- -o [ --output1 ]                      The PE FastQ output file 1 (.fq) (default: EARRINGS_2.fq)<br />
- -O [ --output2 ]                      The PE FastQ output file 2 (.fq) (default: EARRINGS_2.fq)<br />
- -a [ --adapter1 ]                     Default adapter 1 when auto-detect
                                        fails. (default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)<br />
- -A [ --adapter2 ]                     Default adapter 2 when auto-detect
                                        fails. (default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA)<br />
- -t [ --thread ]                       The number of threads used to run the program. (default: 1)<br />
- -m [ --min_length ]                   Abort the read if the length of the read is less than m. (default: 0)<br />
- -l [ --adapter_loc ]                  Specify the location of the adapter. (default: tail)<br />
- -M [ --rc_thres ]                     Setting the threshold of reverse complement check. (default: 0.7)<br />
- -s [ --ss_thres ]                     Setting the threshold of gene portion check. (default: 0.9)<br />
- -S [ --as_thres ]                     Setting the threshold of adapter portion check. (default: 0.8)<br />
- -f [ --prune_factor ]                 Prune factor used when assembling adapters using the de Bruijn graph.
                                        kmer occurence lower than the prune factor will be aborted. (default: 0.03)<br />
- --sensitive                           Sensitive mode can be used when the user is sure that the dataset 
                                        contains adapters. Under sensitive mode, we do not restrict the 
                                        minimum number of kmers when assembly adapters. By default, the 
                                        minimum number of kmers must exceed 10.<br/>
- -F [ --fasta ]                        Specify input file type as FastA. (default input file format: FastQ)<br />
- -b [ --bam_input ] arg                Detect and trim off adapters from a BAM file.<br />


## Run Simulation
The simulation is carried out using a modified version of pIRS(profile based Illumina pair-end Reads Simulator).
Please build pIRS first before using it:
```sh
> cd simulator/pIRS
> make
```
pIRS simulates Illumina PE reads from a reference genome. In all the benchmarking except real data benchmarking,
we use GRCh38 chr1 as reference genome to simulate reads.
One can run the simulation by:
```sh
# python3 script.py [seed] -r [ref_genome] -n [read_num] -m [insert_mean] -std [insert_td] -ad1 [ad1] -ad2 [ad1]
# To generate reads without adapters:
> python3 script.py 0 -r hg38_chr1.fa -n 10000 -m 150 -std 50 -ad1 "" -ad2 ""
# To generate reads with default adapters:
> python3 script.py 0 -r hg38_chr1.fa -n 10000 -m 150 -std 50
```
## Run Benchmarking
Before running benchmarking, please install all the prerequisites
and set the locations of all the executions to $PATH:
1. AdapterRemoval ver. 2.3.0<br />
2. skewer ver. 0.2.2<br />
3. Cutadapt ver. 2.4<br />
4. AKATrim ver. 1.3.3<br />
5. Trimmomatic ver. 0.39<br />
6. SeqPurge ver. 2019_11<br />
7. fastp ver. 0.20.0 - commit 6ff0ffa<br />
8. atropos ver. 1.1.21 - commit 2b15c77<br />
9. PEAT ver. 1.2.5<br />
10. python ver. 3.6<br />
11. our modified version of pIRS (see Run Simulation part)<br />
12. picard 2.23.0<br/>

For speed, memory, performance and adapter benchmarking, please download the first chromosome of GRCh38 from UCSC genome browser.<br />
For running real data benchmarking, please download the following four datasets and their corresponding reference genome(.fa) from ncbi.
1. SRR529095: RNA IP-Seq of Argonaute proteins from Homo sapiens
2. SRR014866: miRNA of C.elegans
3. SRR330569: RNA-Seq of Gonads and Carcasses in D. simulans and D. pseudoobscura
4. SRR5000681: ATAC-Seq of Early Embryo replicate in C.elegans

Also, don't forget to install the aligners:
1. Bowtie2 v2.4.1<br />
2. HISAT2 v2.2.0<br />

After that, modify path-related variables in path.py to the path where the downloaded datasets were placed.<br />
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
```
Jui-Hung Hung <juihunghung@gmail.com>

Ting-Husan Wang <phoebewangintw@gmail.com>

Cheng-Ching Huang <ken5112840801@gmail.com>
```
