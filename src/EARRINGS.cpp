#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <experimental/filesystem>
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
    > EARRINGS single -p index_prefix --skewer input1.fq
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
        size_t earrings_argc(0);
        for (size_t i(2); i < argc; ++i)
        {
            if(std::string(argv[i]) == "--skewer")
            {
                earrings_argc = i + 1;
                break;
            }
        }
        
        std::vector<const char*> skewer_argv(argc - earrings_argc + 3);
        skewer_argv[0] = "skewer";  // skewer is required to install beforehead.

        for (size_t i(1); i < skewer_argv.size() - 2; ++i)
        {
            skewer_argv[i] = argv[i + earrings_argc - 1];
        }
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

        ifs_name[0] = std::string{para.input[0]};
        std::cout << "input: " << ifs_name[0] << "\n";
        
        if (!bam_fname.empty())
        {
            std::cerr << "Processing BAM file...\n";
            auto num_records = Process_uBAMs::extract_reads_from_uBAMs(
                                                    bam_fname
                                                  , ifs_name[0]);
            std::cerr << "Finish processing BAM file!\n";
            // check if the number of unaligned BAM records is gt than DETECT_N_READS
            if (num_records < DETECT_N_READS)
            {
                std::cerr << "Warning: Too few unaligned BAM records: " << num_records << "\n";
            }
        }

		auto adapter_info = seat_adapter_auto_detect(ifs_name[0], para.nThreads);  // auto-detect adapter 
        skewer_argv[skewer_argv.size() - 2] = "-x";
        skewer_argv[skewer_argv.size() - 1] = std::get<0>(adapter_info).c_str();
        if (std::get<1>(adapter_info))
        {
            skewer_argv.emplace_back("-C");
        }

        skewer::main(skewer_argv.size(), skewer_argv.data());
    }
    else if (std::string(argv[1]) == "paired")
    {
        init_paired(argc, argv);

        if (ifs_name[0].find(".gz") == ifs_name[0].size() - 3)
            is_gz_input = true;
        if (ofs_name[0].find(".gz") == ofs_name[0].size() - 3)
            is_gz_output = true;

        if (!bam_fname.empty())
        {
            auto num_records = Process_uBAMs::extract_reads_from_uBAMs(
                                                    bam_fname
                                                  , ifs_name[0]
                                                  , ifs_name[1]);
            // check if the number of unaligned BAM records is gt than DETECT_N_READS
            if (num_records < DETECT_N_READS)
            {
                std::cerr << "Warning: Too few unaligned BAM records: " << num_records << "\n";
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
    *********************************************************************************
    +-----------+
    |Build index|
    +-----------+
    > EARRINGS build -r ref_path -p index_prefix
    > eg. EARRINGS build -r hg38.fa -p earrings_idx
    *********************************************************************************
    )";
    boost::program_options::options_description opts {usage};
    try
    {
        opts.add_options ()
        ("help,h", "Display help message and exit.")
        ("ref_path,r",
         boost::program_options::
            value<std::string>()->required(),
            "Path to the reference genome.")
        ("index_prefix,p",
         boost::program_options::
            value<std::string>()->required(),
            "The index prefix for the built index table.");

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
    
    > EARRINGS single -p earrings_idx --skewer test_file/test_paired1.fq
    *********************************************************************************
    )";
    
    boost::program_options::options_description opts {usage};
    try
    {
        opts.add_options ()
        ("index_prefix,p",
         boost::program_options::
            value<std::string>()->required(),
            "The index prefix for prebuilt index table. (required)")
        ("skewer,s", "Skewer flag, options after this would be fed to Skewer, such as input file name and the number of thread used to run the program. These two parameters will also be used by EARRINGS. Moreover, EARRINGS will pass the auto-detected adapter sequence to Skewer.")
        ("help,h", "Display help message and exit.")
        ("seed_len,d",
         boost::program_options::
            value<size_t>()->default_value(50),
            "The first seed_len bases at 5' portion is viewed as seed when conducting alignment. EARRINGS allows at most one mismatch in the seed portion if enable_mismatch is set to true. Reads will be aborted if more than one mismatch is found in the seed portion. If one mismatch is found outside the seed region, the remainder is reported as a tail. It is recommended to set the seed_len to 18 for very short reads like miRNA, otherwise, it is recommended to set it to 50. (default: 50)")
        ("max_align,m",
         boost::program_options::
            value<size_t>()->default_value(0),
            "Maximum number of candidates used in seed finding stage. (default: 0, not limited)")
        ("enable_mismatch,e",
         boost::program_options::
            value<bool>()->default_value(true),
            "Enable/disable mismatch when conducting seed finding. (default: true)")
        ("prune_factor,f",
         boost::program_options::
            value<float>()->default_value(0.03),
            "Prune factor used when assembling adapters using the de Bruijn graph. kmer frequency lower than the prune factor will be aborted.(default: 0.03)")
        ("fasta,F", "Specify input file type as FastA. (Default input file format: FastQ)")
        ("adapter,a",
        boost::program_options::
            value<std::string>(&DEFAULT_ADAPTER1)->default_value(DEFAULT_ADAPTER1),
        "Alternative adapter if auto-detect mechanism fails.. (default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)")
        ("sensitive", "Sensitive mode can be used when the user is sure that the dataset contains adapters. Under sensitive mode, we do not restrict the minimum number of occurence of kmers when assembly adapters. By default, the minimum number of occurence of kmers must exceed 10.")
        ("bam_input,b", 
        boost::program_options::
            value<std::string>(&bam_fname)->default_value(""),
            "Transform reads in a BAM file into a FastA file. Then trim off adapters from the FastA file. The file name of the untrimmed Fasta file is specified by the input file name for Skewer, while the file name of the trimmed Fasta file is also specified by the output file name for Skewer.")
        ("UMI,u", "Estimate the size of UMI sequences. The result will be printed on the screen.");

        
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
            if (!std::experimental::filesystem::exists(index_prefix + ".table") || 
                !std::experimental::filesystem::exists(index_prefix + ".rc_table"))
            {
                throw std::runtime_error(
                        "Index " + index_prefix + ".table or " + index_prefix + ".rc_table "
                        "does not exist! Please build index first."
                    );
            }
        }
        
        if (vm.count("seed_len"))
        {
            seed_len = vm["seed_len"].as<size_t>();
        }

        if (vm.count("max_align"))
        {
            min_multi = vm["max_align"].as<size_t>();
        }

        if (vm.count("enable_mismatch"))
        {
            enable_mismatch = vm["enable_mismatch"].as<bool>();
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
        
        if (vm.count("fasta") || !bam_fname.empty())
        {
            is_fastq = false;
        }

        if (vm.count("sensitive"))
        {
            is_sensitive = true;
        }

        if (vm.count("UMI"))
        {
            estimate_umi_len = true;
        }

        if (!vm.count("skewer"))
        {
            throw std::runtime_error("--skewer flag is mandatory, "
                "which follows by the arguments that are passed into skewer program.");
        }

        std::cout << "Index prefix: " << index_prefix << std::endl;
        std::cout << "Seed length: " << seed_len << ", Maximum alignment: " << min_multi << ", sensitive mode: " << is_sensitive << std::endl;
        std::cout << "Enable mismatch: " << enable_mismatch << ", Prune factor: " << prune_factor << ", is fastq: " << is_fastq << std::endl;
        std::cout << "Default adapter: " << DEFAULT_ADAPTER1 << std::endl;

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
	*********************************************************************************
	+----------+
	|Paired End|
	+----------+
	Paired-end adapter trimming.
	EARRINGS takes paired-end FastQ/FastA format input files (dual files), and outputs 
	adapter removed FastQ/FastA format output files (dual files).

	> EARRINGS paired -i test_file/test_paired1.fq -I test_file/test_paired2.fq -t thread
	*********************************************************************************
)"; 
    
    boost::program_options::options_description opts {usage};
	
    try 
    {
        opts.add_options ()
        ("input1,i", 
         boost::program_options::
            value<std::string>(&ifs_name[0])->required(), 
         "The Paired-End FastQ input file 1 (.fq) (required)")
        ("input2,I", 
         boost::program_options::
            value<std::string>(&ifs_name[1])->required(), 
         "The Paired-End FastQ input file 2 (.fq) (required)")
        ("help,h", "Display help message and exit.")
        ("output1,o",
         boost::program_options::
            value<std::string>(&ofs_name[0])->default_value("EARRINGS_1.fq"),
        "The Paired-End FastQ output file 1 (.fq) (default: EARRINGS_2.fq)")
        ("output2,O",
         boost::program_options::
            value<std::string>(&ofs_name[1])->default_value("EARRINGS_2.fq"),
        "The Paired-End FastQ output file 2 (.fq) (default: EARRINGS_2.fq)")
        ("adapter1,a",
        boost::program_options::
            value<std::string>(&DEFAULT_ADAPTER1)->default_value(DEFAULT_ADAPTER1),
        "Default adapter 1 when auto-detect fails. (default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC)")
        ("adapter2,A",
        boost::program_options::
            value<std::string>(&DEFAULT_ADAPTER2)->default_value(DEFAULT_ADAPTER2),
        "Default adapter 2 when auto-detect fails. (default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA)")
        ("thread,t", 
         boost::program_options::
            value<size_t>()->default_value(1), 
         "The number of threads used to run the program. (default: 1)")
        ("min_length,m",
         boost::program_options::
            value<size_t>()->default_value(0),
            "Abort the read if the length of the read is less than m. (default: 0)")
        ("rc_thres,M",
         boost::program_options::
            value<float>()->default_value(0.7),
            "Mismatch threshold applied in reverse complement scan. (default: 0.7)")
        ("ss_thres,s",
         boost::program_options::
            value<float>()->default_value(0.9),
            "Mismatch threshold applied in gene portion check. (default: 0.9)")
        ("as_thres,S",
         boost::program_options::
            value<float>()->default_value(0.8),
            "Mismatch threshold applied in adapter portion check. (default: 0.8)")
        ("prune_factor,f",
         boost::program_options::
            value<float>()->default_value(0.03),
            "Prune factor used when assembling adapters using the de Bruijn graph. kmer frequency lower than prune factor will be aborted.(default: 0.03)")
        ("sensitive", "Sensitive mode can be used when the user is sure that the dataset contains adapters. Under sensitive mode, we do not restrict the minimum number of occurence of kmers when assembly adapters. By default, the minimum number of occurence of kmers must exceed 10.")
        ("fasta,F", "Specify input file type as FastA. (default input file format: FastQ)")
        ("bam_input,b", 
         boost::program_options::
            value<std::string>(&bam_fname)->default_value(""),
            "Transform reads in a BAM file into two FastA files. Then trim off adapters from the FastA files. The file names of the untrimmed Fasta files are determined by input1 and input2, while the file names of the trimmed Fasta files are determined by output1 and output2.");
        
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

        if (vm.count("fasta") || !bam_fname.empty())
        {
            is_fastq = false;
        }

        if (vm.count("sensitive"))
        {
            is_sensitive = true;
        }

        match_rate = (vm["rc_thres"].as<float>() > 0.0 && vm["rc_thres"].as<float>() < 1.0) ? vm["rc_thres"].as<float>() : 0.7;
        seq_cmp_rate = (vm["ss_thres"].as<float>() > 0.0 && vm["ss_thres"].as<float>() < 1.0) ? vm["ss_thres"].as<float>() : 0.9;
        adapter_cmp_rate = (vm["as_thres"].as<float>() > 0.0 && vm["as_thres"].as<float>() < 1.0) ? vm["as_thres"].as<float>() : 0.8;

        std::cout << "# of threads: " << thread_num << ", Prune factor: " << prune_factor << ", Minimum output length: " << min_length << std::endl;
        std::cout << "Match rate: " << match_rate << ", Seq cmp rate: " << seq_cmp_rate << ", Adapter cmp rate: " << adapter_cmp_rate << std::endl; 
        std::cout << "is fastq: " << is_fastq << ", sensitive mode: " << is_sensitive << std::endl;
        std::cout << "Default adapter1: " << DEFAULT_ADAPTER1 << std::endl;
        std::cout << "Default adapter2: " << DEFAULT_ADAPTER2 << std::endl;
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

