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
Since SE algorithm uses alignment-based method to find possible adapter, one have to build index first.
#### SE build
Before conducting SE adapter trimming, **one has to prebuild index** once for a specific reference genome.
```sh
# ./EARRINGS build -r [ref_path] -p [index_prefix]
> ./EARRINGS build -r ref/hg38_chr1.fa -p earrings_hg38_chr1.
```
#### run SE
EARRINGS first detects adapter then feeds detected adapter to skewer.
```sh
# ./EARRINGS single -p [index_prefix] --skewer [input_file] [skewer_parameters]
> ./EARRINGS single -p path_to_index. --skewer ../test_data/has_adapter_1.fq
```
#### SE parameters
-h [ --help ]                         Display this help message and exit.<br />
-p [ --index_prefix ] arg             The index path for single-end.<br />
-d [ --seed_len ] arg (=50)           Seed length for finding adapters for
                                      single end reads. For very short reads like miRNA, 
                                      it is recommended to set seed_len to 18. (default: 50)<br />
-m [ --max_align ] arg (=0)           Control the maximum number of alignment
                                      to abort the reads. (default: 0, not
                                      limited)<br />
-e [ --enable_mismatch ] arg (=1)     Enable/disable mismatch when doing seed
                                      finding. (default: true)<br />
-f [ --prune_factor ] arg (=0.03)     Prune factor used when assembling adapters using debruijn graph.
                                      kmer number lower than prune factor will be aborted.(default: 0.03)<br />
-F [ --fasta ]                        Specify that the input is FastA.
                                      (Default input file format: FastQ)<br />
-a [ --adapter1 ] arg (=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)
                                      Default adapter used when auto-detect fails.<br />
--sensitive                           Sensitive mode can be used when the user is sure that the dataset 
                                      contains adapters. Under sensitive mode, we do not restrict the 
                                      minimum number of kmers when assembly adapters. By default, the 
                                      minimum number of kmers must exceed 10.<br/>
-s [ --skewer ]                       skewer flag, which follows by skewer's
                                      program options<br />

### PE
#### run PE
```sh
# ./EARRINGS paired -i [input1] -I [input2] -t [thread_num]
> ./EARRINGS paired -i ../test_data/has_adapter_1.fq -I ../test_data/has_adapter_2.fq
```
#### PE parameters
-h [ --help ]                         Display this help message and exit.<br />
-i [ --input1 ] arg                   The PE input FastQ file 1 (.fq)<br />
-I [ --input2 ] arg                   The PE input FastQ file 2 (.fq)<br />
-o [ --output1 ] arg (=EARRINGS_1.fq)     The PE output FastQ file 1 (.fq)<br />
-O [ --output2 ] arg (=EARRINGS_2.fq)     The PE output FastQ file 2 (.fq)<br />
-a [ --adapter1 ] arg (=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)
                                      Default adapter 1 when auto-detect
                                      fails.<br />
-A [ --adapter2 ] arg (=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA)
                                      Default adapter 2 when auto-detect
                                      fails.<br />
-t [ --thread ] arg (=1)              The number of threads to use.<br />
-m [ --min_length ] arg (=0)          Abort reads if reads is less than m.<br />
-l [ --adapter_loc ] arg (=tail)      Adapter locates at 5'(head) or
                                      3'(tail). (default: tail)<br />
-M [ --match_rate ] arg (=0.7)        Sequence match rate when detecting
                                      possible adapter positions.(default:
                                      0.7)<br />
-s [ --seq_cmp_rate ] arg (=0.9)      Sequence similariy when comparing first
                                      strand to the reverse complement of
                                      second strand.(default: 0.9)<br />
-S [ --adapter_cmp_rate ] arg (=0.8)  Adapter similariy when comparing with
                                      detected adapter.(default: 0.8)<br />
-f [ --prune_factor ] arg (=0.03)     Prune factor used when assembling adapters using debruijn graph.
                                      kmer number lower than prune factor will be aborted.(default: 0.03)<br />
--sensitive                           Sensitive mode can be used when the user is sure that the dataset 
                                      contains adapters. Under sensitive mode, we do not restrict the 
                                      minimum number of kmers when assembly adapters. By default, the 
                                      minimum number of kmers must exceed 10.<br/>
-a [ --fasta ]                        Specify that the input is FastA.
                                      (Default input file format: FastQ)<br />


## Run Simulation
The simulation is carried out using a modified version of pIRS(profile based Illumina pair-end Reads Simulator).
Please build it first before using it:
```sh
> cd simulator/pIRS
> make
```
pIRS simulates Illumina PE reads from a reference genome. In all the benchmarking except real data benchmarking,
we use hg38 chr1 as reference genome to simulate reads.
Then one can run the simulation by:
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
Before running real data benchmarking, please download these four sets of dataset and their corresponding reference genome(.fa) from ncbi.
1. SRR529095: RNA IP-Seq of Argonaute proteins from Homo sapiens
2. SRR014866: miRNA of C.elegans
3. SRR330569: RNA-Seq of Gonads and Carcasses in D. simulans and D. pseudoobscura
4. SRR5000681: ATAC-Seq of Early Embryo replicate in C.elegans

Also, don't forget to download the aligners and build indices for them:
1. Bowtie2 v2.4.1<br />
2. HISAT2 v2.2.0<br />

After that, modify path-related variables in benchmark_real_data.py and path.py
to the path where you put all your downloaded datasets/prebuild indices.<br />
After that, you can conduct real data benchmarking by:
```sh
> python3 benchmark_real_data.py
```