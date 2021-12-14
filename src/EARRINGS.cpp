#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <boost/program_options.hpp>
#include <EARRINGS/SE/SE_auto_detect.hpp>
#include <EARRINGS/PE/PE_trimmer.hpp>
#include <EARRINGS/estimate_adapter_from_BAMS.hpp>
#include <EARRINGS/common.hpp>
#include "skewer/main.hpp"
#include "skewer/parameter.h"

void init_single(int argc, const char* argv[]);
void init_paired(int argc, const char* argv[]);
void init_build(int argc, const char* argv[]);

int main(int argc, const char* argv[])
{
    std::chrono::time_point<std::chrono::steady_clock> start, end;
    start = std::chrono::steady_clock::now();

    std::string ver = GET_EARRINGS_VERSION;
    std::string help = R"(
    *********************************************************************************
    +--------+
    |EARRINGS|
    +--------+
    EARRINGS v)" + ver + 
    R"( is an adapter trimmer with no a priori knowledge of adapter sequences.
    Usage:
    (1) Build index for reference sequence. This step is only apply to single-end ada-
        pter detection. 
    > EARRINGS build -r ref_path -p index_prefix
    (2) Adapter trimming
    > EARRINGS single -p index_prefix -1 input1.fq -t thread_num
    > EARRINGS paired -i input1.fq -I input2.fq -t thread_num 
    See EARRINGS single/paired --help for more information about the parameters.
    *********************************************************************************
    )";

    if (argc == 1)
    {
        std::cout << help << "\n";
        return 0;
    }

    if (std::string(argv[1]) == "single")
    {
        skewer::cParameter para;
        char errMsg[256];
        init_single(argc, argv);

        if (is_bam)
        {
            std::string tmp_name("/tmp/EARRINGS_bam_reads.tmp");

            std::cerr << "Processing BAM file...\n";
            auto num_records = Process_uBAMs::extract_reads_from_uBAMs(
                                                    ifs_name[0]
                                                  , tmp_name);
            ifs_name[0] = tmp_name;
            std::cerr << "Finish processing BAM file!\n";
            // check if the number of BAM records is gt than DETECT_N_READS
            if (num_records < DETECT_N_READS)
            {
                std::cerr << "Warning: Too few BAM records: " << num_records << "\n";
            }
        }

        auto adapter_info = seat_adapter_auto_detect(ifs_name[0], para.nThreads);  // auto-detect adapter 
        
        // input, output, min_len, thread, adapter, quiet flag
        std::vector<const char*> skewer_argv(10);
        skewer_argv[0] = "skewer";  // skewer is required to install beforehead.
        skewer_argv[1] = ifs_name[0].c_str(); // input
        skewer_argv[2] = "-o";  // output
        skewer_argv[3] = ofs_name[0].c_str();
        skewer_argv[4] = "-l";  // min_len
        skewer_argv[5] = std::to_string(min_length).c_str();
        skewer_argv[6] = "-t";  // thread
        skewer_argv[7] = std::to_string(thread_num).c_str();
        int32_t iRet = para.GetOpt(skewer_argv.size() - 2, skewer_argv.data(), errMsg);
        // copy from skewer's main program.
        if (iRet < 0)
        {
            const char * program = strrchr(argv[0], '/');
            program = (program == NULL) ? argv[0] : (program + 1);
            if (iRet == -1)
            {
                if(para.bEnquireVersion)
                {
                    para.PrintVersion(stdout);
                    return 0;
                }
                para.PrintUsage(program, stdout);
            }
            else
            {
                fprintf(stderr, "%s (%s): %s\n\n", program, para.version, errMsg);
                para.PrintSimpleUsage(program, stderr);
            }
            return 1;
        }
        skewer_argv[8] = "-x";
        skewer_argv[9] = std::get<0>(adapter_info).c_str();
        if (std::get<1>(adapter_info))
        {
            skewer_argv.emplace_back("-C");
        }

        skewer::main(skewer_argv.size(), skewer_argv.data());
    }
    else if (std::string(argv[1]) == "paired")
    {
        init_paired(argc, argv);

        if (is_bam)
        {
            std::string tmp_name1("/tmp/EARRINGS_bam_reads1.tmp");
            std::string tmp_name2("/tmp/EARRINGS_bam_reads2.tmp");
            // extract reads are fasta by default
            auto num_records = Process_uBAMs::extract_reads_from_uBAMs(
                                                    ifs_name[0]
                                                  , tmp_name1
                                                  , tmp_name2);
            
            // check if the number of BAM records is gt than DETECT_N_READS
            ifs_name[0] = tmp_name1;
            ifs_name[1] = tmp_name2;
            if (num_records < DETECT_N_READS)
            {
                std::cerr << "Warning: Too few BAM records: " << num_records << "\n";
            }
        }
        PE_trim();
    }
    else if (std::string(argv[1]) == "build")
    {
        init_build(argc, argv);
    }
    else
    {
        std::cout << help << "\n";
        return 0;
    }

    end = std::chrono::steady_clock::now();
    std::chrono::duration<double> sec(end - start);
    std::cout << sec.count() << " sec\n";

    return 0;
}


void init_build(int argc, const char* argv[])
{
    std::string usage = R"(
*****************************************************************************
+-----------+
|Build index|
+-----------+
Before conducting single-end adapter trimming , EARRINGS has to prebuild the 
index once for a specific reference which is the source of the target reads.

> EARRINGS build -r hg38.fa -p earrings_idx
*****************************************************************************
)";
    boost::program_options::options_description opts {usage};
    try
    {
        opts.add_options ()
        ("help,h", "Display help message and exit.")
        ("ref_path,r",
         boost::program_options::
            value<std::string>()->required(),
            "Path to the reference sequence, which is the source of reads. (required)")
        ("index_prefix,p",
         boost::program_options::
            value<std::string>()->required(),
            "An user-defined index prefix for index table. (required)");

        boost::program_options::variables_map vm;
        boost::program_options::store (
            boost::program_options::parse_command_line(
            argc, argv, opts
            ), vm
        );

        boost::program_options::notify(vm);
        if (vm.count("help"))
        {
            std::cout << usage << "\n";
            exit(0);
        }

        if (vm.count("ref_path") && vm.count("index_prefix"))
        {
            tailor::build_tables build;
            build(vm["ref_path"].as<std::string>()
                , vm["index_prefix"].as<std::string>());
        }
    } 
    catch (std::exception& e) 
    {
        std::cerr << "Error: " << e.what() << std::endl 
            << opts << std::endl;
        exit (1);
    } 
    catch (...) 
    {
        std::cerr << "Unknown error!" << std::endl 
            << opts << std::endl;
        exit (1);
    }

}

void init_single(int argc, const char* argv[])
{
    std::string usage = R"(
*********************************************************************************
+----------+
|Single End|
+----------+
Single-end adapter trimming. 
EARRINGS detects adapter using alignment-based method. Thus, it is necessary to 
prebuild the index first. For downstream adapter trimming, it is conducted using 
Skewer with adapter parameters passed by EARRINGS automatically.

> EARRINGS single -p earrings_idx -1 test_file/test_paired1.fa
> EARRINGS single -p earrings_idx -1 test_file/test_paired1.fq
> EARRINGS single -p earrings_idx -1 test_file/test_paired1.fq.gz
*********************************************************************************
    )";
    
    boost::program_options::options_description opts {usage};
    try
    {
        opts.add_options ()
        ("index_prefix,p",
         boost::program_options::
            value<std::string>()->required(),
            "The index prefix for pre-built index table. (required)")
        ("input1,1", 
         boost::program_options::
            value<std::string>(&ifs_name[0])->required(), 
            "The file path of Single-End reads. (required)")
        ("help,h", 
            "Display help message and exit.")
        ("seed_len,d",
         boost::program_options::
            value<size_t>()->default_value(50),
            "The first --seed_len bases are seen as seed and allows 1 mismatch at most, "
            "or do not allow any mismatch if --no_mismatch is set. The sequence follows "
            "first mismatch out of the seed portion will be reported as a tail.\n"
            "It's recommended to set this to 18 for very short reads like miRNA, "
            "otherwise, it is recommended to set to 50.")
        ("output,o",
         boost::program_options::
            value<std::string>(&ofs_name[0])->default_value("trimmed_se"),
            "The file prefix of Single-End FastQ output.")
        ("min_length,m",
         boost::program_options::
            value<size_t>()->default_value(0),
            "Skip the read if the length of the read is less than --min_length after trimming.")
        ("thread,t", 
         boost::program_options::
            value<size_t>()->default_value(1), 
            "The number of threads used to run the program.")
        ("max_align,M",
         boost::program_options::
            value<size_t>()->default_value(0),
            "Maximum number of candidates used in seed finding stage, 0 means unlimited.")
        ("no_mismatch,e",
         boost::program_options::
            bool_switch(&no_mismatch),
            "By default, EARRINGS can tolerate 1 error base at most if be set as true, "
            "this flag can disable this mismatch toleration mechenism.")
        ("prune_factor,f",
         boost::program_options::
            value<float>()->default_value(0.03),
            "Prune factor used when assembling adapters using the de-brujin graph. Kmer "
            "frequency lower than this value will be skipped.")
        ("adapter,a",
         boost::program_options::
            value<std::string>(&DEFAULT_ADAPTER1)->default_value(DEFAULT_ADAPTER1),
            "Alternative adapter if auto-detect mechanism fails.")
        ("sensitive", 
            "By default, minimum number of kmers must exceed 10 during assembly adapters. "
            "However, if user have confidence that the dataset contains adapters, sensitive "
            "mode would be more suitable.\n"
            "Under sensitive mode, minimum number of kmers (--prune_factor) would not be "
            "restricted.")
        ("UMI,u", 
            "Estimate the size of UMI sequences, results will be printed to console by "
            "default.");

        
        boost::program_options::variables_map vm;
        boost::program_options::store (
            boost::program_options::command_line_parser(
            argc, argv
            ).options(opts).allow_unregistered().run(), vm
        );

        boost::program_options::notify(vm);

        if (vm.count("help"))
        {
            std::cout << usage << "\n";
            exit(0);
        }

        if (vm.count("index_prefix"))
        {
            index_prefix = vm["index_prefix"].as<std::string>();
            if (!std::filesystem::exists(index_prefix + ".table") || 
                !std::filesystem::exists(index_prefix + ".rc_table"))
            {
                throw std::runtime_error(
                        "Index " + index_prefix + ".table or " + index_prefix + ".rc_table "
                        "does not exist! Please build index first."
                    );
            }
        }
        
        thread_num = vm["thread"].as<size_t>();
        if (thread_num > 32) thread_num = 32;

        min_length = vm["min_length"].as<size_t>();
        
        if (vm.count("seed_len"))
        {
            seed_len = vm["seed_len"].as<size_t>();
        }

        if (vm.count("max_align"))
        {
            min_multi = vm["max_align"].as<size_t>();
        }

        if (vm.count("prune_factor"))
        {
            prune_factor = vm["prune_factor"].as<float>();
            if (prune_factor >= 1.0 || prune_factor < 0.0)
            {
                std::cout << "prune_factor must be between [0.0, 1.0). Setting it to 0.03.\n";
                prune_factor = 0.03;
            }
        }
        if (vm.count("sensitive"))
        {
            is_sensitive = true;
        }

        if (vm.count("UMI"))
        {
            estimate_umi_len = true;
        }

        std::string fa_ext(".fa"), fasta_ext(".fasta");
        if (ifs_name[0].find(".gz") == ifs_name[0].size() - 3) {
            is_gz_input = true;
            fa_ext.append(".gz");
            fasta_ext.append(".gz");
        }
        // if (ofs_name[0].find(".gz") == ofs_name[0].size() - 3)
        //     is_gz_output = true;

        if (ifs_name[0].find(".bam") == ifs_name[0].size() - 4) {
            is_bam = true;
            is_fastq = false;
            fa_ext.append(".bam");
            fasta_ext.append(".bam");
        }
        if (ifs_name[0].find(".ubam") == ifs_name[0].size() - 5) {
            is_bam = true;
            is_fastq = false;
            fa_ext.append(".ubam");
            fasta_ext.append(".ubam");
        }
        if (ifs_name[0].find("bam.gz") == ifs_name[0].size() - 6) {
            std::cerr << "Error: Bam mode and gz mode are not compatible.\n";
            exit(1);
        }

        if (is_bam || 
            ifs_name[0].find(    fa_ext ) == ifs_name[0].length() -    fa_ext.length() || 
            ifs_name[0].find( fasta_ext ) == ifs_name[0].length() - fasta_ext.length() )
        {
            ofs_name[0] += ".fasta";
            is_fastq = false;
        }
        else
            ofs_name[0] += ".fastq";

        std::cout << std::boolalpha;
        std::cout << "Index prefix: " << index_prefix << std::endl;
        std::cout << "Input file name: " << ifs_name[0] << std::endl;
        std::cout << "Output file name: " << ofs_name[0] << std::endl;
        std::cout << "# of threads: " << thread_num << std::endl;
        std::cout << "Is fastq: " << is_fastq << ", Is gz input: " << is_gz_input << ", Is bam: " << is_bam << std::endl;
        std::cout << "Seed length: " << seed_len << ", Max alignment: " << min_multi << ", No mismatch: " << no_mismatch << std::endl;
        std::cout << "Prune factor: " << prune_factor << ", Sensitive mode: " << is_sensitive << std::endl;
        std::cout << "Min length: " << min_length << ", UMI: " << estimate_umi_len << std::endl;
        std::cout << "Default adapter: " << DEFAULT_ADAPTER1 << std::endl;
        std::cout << std::noboolalpha;
    }
    catch (std::exception& e) 
    {
        std::cerr << "Error: " << e.what() << std::endl 
            << opts << std::endl;
        exit (1);
    } 
    catch (...) 
    {
        std::cerr << "Unknown error!" << std::endl 
            << opts << std::endl;
        exit (1);
    }
}


void init_paired(int argc, const char* argv[])
{
    std::string usage = R"(
*******************************************************************************************
+----------+
|Paired End|
+----------+
Paired-end adapter trimming.
EARRINGS takes paired-end FastQ/FastA format input files (dual files), and outputs 
adapter removed FastQ/FastA format output files (dual files).

> EARRINGS paired -1 test_file/test_paired1.fa -2 test_file/test_paired2.fa -t thread
> EARRINGS paired -1 test_file/test_paired1.fq -2 test_file/test_paired2.fq -t thread
> EARRINGS paired -1 test_file/test_paired1.fq.gz -2 test_file/test_paired2.fq.gz -t thread
*******************************************************************************************
)"; 
    
    boost::program_options::options_description opts {usage};
	
    try 
    {
        opts.add_options ()
        ("input1,1", 
         boost::program_options::
            value<std::string>(&ifs_name[0])->required(), 
            "The Paired-end reads input file 1.")
        ("input2,2", 
         boost::program_options::
            value<std::string>(&ifs_name[1])->required(), 
            "The Paired-end reads input file 2.")
        ("help,h", 
            "Display help message and exit.")
        ("output,o",
         boost::program_options::
            value<std::string>(&ofs_name[0])->default_value("trimmed_pe"),
            "The Paired-End FastQ output file prefix.")
        ("adapter1,a",
         boost::program_options::
            value<std::string>(&DEFAULT_ADAPTER1)->default_value(DEFAULT_ADAPTER1),
            "Alternative adapter 1 if auto-detect mechanism fails.")
        ("adapter2,A",
         boost::program_options::
            value<std::string>(&DEFAULT_ADAPTER2)->default_value(DEFAULT_ADAPTER2),
            "Alternative adapter 2 if auto-detect mechanism fails.")
        ("thread,t", 
         boost::program_options::
            value<size_t>()->default_value(1), 
            "The number of threads used to run the program.")
        ("min_length,m",
         boost::program_options::
            value<size_t>()->default_value(0),
            "Skip the read if the length of the read is less than this value.")
        ("rc_thres,M",
         boost::program_options::
            value<float>()->default_value(0.7),
            "Mismatch threshold applied in reverse complement scan.")
        ("ss_thres,s",
         boost::program_options::
            value<float>()->default_value(0.9),
            "Mismatch threshold applied in gene portion check.")
        ("as_thres,S",
         boost::program_options::
            value<float>()->default_value(0.8),
            "Mismatch threshold applied in adapter portion check.")
        ("prune_factor,f",
         boost::program_options::
            value<float>()->default_value(0.03),
            "Prune factor used when assembling adapters using the de Bruijn graph. Kmer "
            "frequency lower than the prune factor will be skipped.")
        ("sensitive", 
            "By default, minimum number of kmers must exceed 10 during assembly adapters. "
            "However, if user have confidence that the dataset contains adapters, sensitive "
            "mode would be more suitable.\n"
            "Under sensitive mode, minimum number of kmers (--prune_factor) would not be "
            "restricted.");
        
        boost::program_options::variables_map vm;
        boost::program_options::store (
            boost::program_options::command_line_parser(
            argc, argv
            ).options(opts).allow_unregistered().run(), vm
        );

        boost::program_options::notify(vm);
        if (vm.count("help"))
        {
            std::cout << opts << std::endl;
            exit(0);
        }

        thread_num = vm["thread"].as<size_t>();
        if (thread_num > 32) thread_num = 32;

        min_length = vm["min_length"].as<size_t>();

        if (vm.count("prune_factor"))
        {
            prune_factor = vm["prune_factor"].as<float>();
            if (prune_factor >= 1.0 || prune_factor < 0.0)
            {
                std::cout << "prune_factor must be between [0.0, 1.0). Setting it to 0.03.\n";
                prune_factor = 0.03;
            }
        }

        if (vm.count("sensitive"))
        {
            is_sensitive = true;
        }

        match_rate = (vm["rc_thres"].as<float>() > 0.0 && vm["rc_thres"].as<float>() < 1.0) ? vm["rc_thres"].as<float>() : 0.7;
        seq_cmp_rate = (vm["ss_thres"].as<float>() > 0.0 && vm["ss_thres"].as<float>() < 1.0) ? vm["ss_thres"].as<float>() : 0.9;
        adapter_cmp_rate = (vm["as_thres"].as<float>() > 0.0 && vm["as_thres"].as<float>() < 1.0) ? vm["as_thres"].as<float>() : 0.8;

        std::string fa_ext(".fa"), fasta_ext(".fasta");
        if (ifs_name[0].find(".gz") == ifs_name[0].size() - 3) {
            is_gz_input = true;
            fa_ext.append(".gz");
            fasta_ext.append(".gz");
        }
        // if (ofs_name[0].find(".gz") == ofs_name[0].size() - 3)
        //     is_gz_output = true;

        if (ifs_name[0].find(".bam") == ifs_name[0].size() - 4) {
            is_bam = true;
            is_fastq = false;
            fa_ext.append(".bam");
            fasta_ext.append(".bam");
        }
        if (ifs_name[0].find(".ubam") == ifs_name[0].size() - 5) {
            is_bam = true;
            is_fastq = false;
            fa_ext.append(".ubam");
            fasta_ext.append(".ubam");
        }
        if (ifs_name[0].find("bam.gz") == ifs_name[0].size() - 6) {
            std::cerr << "Error: Bam mode and gz mode are not compatible.\n";
            exit(1);
        }

        ofs_name[1] = ofs_name[0];
        if (is_bam ||
            ifs_name[0].find(    fa_ext ) == ifs_name[0].length() -    fa_ext.length() || 
            ifs_name[0].find( fasta_ext ) == ifs_name[0].length() - fasta_ext.length() )
        {
            ofs_name[0] += "_1.fasta";
            ofs_name[1] += "_2.fasta";
            is_fastq = false;
        }
        else {
            ofs_name[0] += "_1.fastq";
            ofs_name[1] += "_2.fastq";
        }

        std::cout << std::boolalpha;
        std::cout << "Index prefix: " << index_prefix << std::endl;
        std::cout << "Input file name 1: " << ifs_name[0] << ", Input file name 2:" << ifs_name[1] << std::endl;
        std::cout << "Output file name 1: " << ofs_name[0]<< ", Output file name 2:" << ofs_name[1]  << std::endl;
        std::cout << "# of threads: " << thread_num << std::endl;
        std::cout << "Is fastq: " << is_fastq << ", Is gz input: " << is_gz_input << ", Is bam: " << is_bam << std::endl;
        std::cout << "Prune factor: " << prune_factor << ", Sensitive mode: " << is_sensitive << std::endl;
        std::cout << "Min length: " << min_length << ", UMI: " << estimate_umi_len << std::endl;
        std::cout << "Match rate: " << match_rate << ", Seq cmp rate: " << seq_cmp_rate << ", Adapter cmp rate: " << adapter_cmp_rate << std::endl;
        std::cout << "Default adapter1: " << DEFAULT_ADAPTER1 << std::endl;
        std::cout << "Default adapter2: " << DEFAULT_ADAPTER2 << std::endl;
        std::cout << std::noboolalpha;
    }
    catch (std::exception& e) 
    {
        std::cerr << "Error: " << e.what() << std::endl 
            << opts << std::endl;
        exit (1);
    } 
    catch (...) 
    {
        std::cerr << "Unknown error!" << std::endl 
            << opts << std::endl;
        exit (1);
    }
}

