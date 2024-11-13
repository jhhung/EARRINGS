#pragma once
#include <vector>
#include <string>
#include <filesystem>
#include <EARRINGS/PE/task.hpp>
#include <EARRINGS/PE/buffer_manager.hpp>
#include <EARRINGS/PE/rw_count.hpp>
#include <EARRINGS/PE/trimmer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <tuple>
#include <string_view>
#include <algorithm>

namespace fs = std::filesystem;
using namespace EARRINGS;

namespace EARRINGS
{
template<template<bool, class> class FORMAT, class BITSTR, typename IFS, typename OFS>
class TaskProcessor
{
private:
    template <typename T>
    using remove_cvr_t = std::remove_cv_t<std::remove_reference_t<T>>;
    using BIO_filtering_istream = boost::iostreams::filtering_istream;
    using BIO_filtering_ostream = boost::iostreams::filtering_ostream;
    using FORMAT2BIT = FORMAT<false, BITSTR>;
    fs::path _tmp_path;
    BufferManager _buf_manager;
    Trimmer<PAIRED, FORMAT2BIT> _tr;
    RWCount _rw_count;
    std::vector<fs::path> _tmp_dir;
    std::vector<IFS> _ifs;
    std::vector<OFS> _ofs;
    std::vector<std::string> _default_adapters;
    std::vector<BITSTR> _adapters;
    size_t _detect_n_reads;
    size_t _thread_num;
    size_t _record_line;
    size_t _chunk_size;
    bool _loc_tail;
    bool _is_sensitive;

    void preprocess(FORMAT2BIT&, FORMAT2BIT&, size_t);
    void detect_adapters();
    void trim_reads(Task&);
    bool read_reads(Task&);

    template<class POOL>
    bool read_task(POOL&);
    template<class POOL>
    void trim_task(Task, POOL&);
    void write_task(Task& task);
    void write_all_task();

public:
    TaskProcessor(size_t
                , size_t
                , size_t
                , std::vector<std::string>&
                , std::vector<std::string>&
                , std::tuple<float, float, float>&
                , bool
                , std::vector<std::string>&
                , size_t
                , bool);
    void process();
};

template<template<bool, class> class FORMAT, class BITSTR, typename IFS, typename OFS>
TaskProcessor<FORMAT, BITSTR, IFS, OFS>::TaskProcessor(size_t thread_num
                                                    , size_t chunk_size
                                                    , size_t record_line
                                                    , std::vector<std::string>& ifs_name
                                                    , std::vector<std::string>& ofs_name
                                                    , std::tuple<float, float, float>& trimmer_param
                                                    , bool loc_tail
                                                    , std::vector<std::string>& default_adapter
                                                    , size_t detect_n_reads
                                                    , bool is_sensitive)
    // input/output files
    : _ifs(2)
    , _ofs(2)
    , _record_line(record_line)
    , _chunk_size(chunk_size)
    , _thread_num(thread_num)
    , _detect_n_reads(detect_n_reads)
    , _default_adapters(default_adapter)
{
    _buf_manager.set_chunk_size(chunk_size, _thread_num);
        
    for (size_t i = 0; i < 2; ++i)
    {
        if constexpr (std::is_same_v<remove_cvr_t<IFS>, BIO_filtering_istream>)
        {
            _ifs[i].push(boost::iostreams::gzip_decompressor());
            auto&& src(boost::iostreams::file_source(ifs_name[i], std::ios_base::binary));
            if (!src.is_open())
                throw std::runtime_error("Can't open input gz file normally\n");
            
            _ifs[i].push(src);
            if (!_ifs[i].good())
                throw std::runtime_error("Can't open input gz stream normally\n");
        }
        else
        {
            _ifs[i].open(ifs_name[i]);
            if (!(_ifs[i].is_open() && _ifs[i].good()))
                throw std::runtime_error("Can't open input file normally\n");
        }

        if constexpr (std::is_same_v<remove_cvr_t<OFS>, BIO_filtering_ostream>)
        {
            _ofs[i].push(boost::iostreams::gzip_compressor());
            auto&& sink(boost::iostreams::file_sink(ofs_name[i], std::ios_base::binary));
            if (!sink.is_open())
                throw std::runtime_error("Can't open output gz file normally\n");

            _ofs[i].push(sink);
            if (!_ifs[i].good())
                throw std::runtime_error("Can't open input gz stream normally\n");
        }
        else
        {
            _ofs[i].open(ofs_name[i]);
            if (!(_ofs[i].is_open() && _ofs[i].good()))
                throw std::runtime_error("Can't open output file normally\n");
        }
    }

    // create tmp dir
    _tmp_path = fs::temp_directory_path();
    _tmp_dir = std::vector<fs::path>{fs::path(_tmp_path / fs::path("tmp1"))
                                   , fs::path(_tmp_path / fs::path("tmp2"))};
    fs::create_directories(_tmp_dir[0]);
    fs::create_directories(_tmp_dir[1]);

    // init adapters
    _adapters = decltype(_adapters)(2, BITSTR(""));

    // trimmer param
    _tr.trait_parm.match_rate = std::get<0>(trimmer_param);
    _tr.trait_parm.seq_cmp_rate = std::get<1>(trimmer_param);
    _tr.trait_parm.adapter_cmp_rate = std::get<2>(trimmer_param);

    // location of tails
    _loc_tail = loc_tail;

    // sensitive mode
    _is_sensitive = is_sensitive;
}

template<template<bool, class> class FORMAT, class BITSTR, typename IFS, typename OFS>
void TaskProcessor<FORMAT, BITSTR, IFS, OFS>::process()
{
    detect_adapters();
    bool eof = false;
    auto pool = nucleona::parallel::make_asio_pool(_thread_num);
	while (true)
	{
		if (!eof)
		{
			auto fut = pool.submit([this, &pool](){ return read_task(pool); } );
			eof = fut.sync();
		}
		else
		;
        
        if (_rw_count.all_finished)
        {
            std::cerr << "";
        }

        if (_rw_count.all_finished)
        {
            break;
        }
	}

    pool.flush();
}


template<template<bool, class> class FORMAT, class BITSTR, typename IFS, typename OFS>
void TaskProcessor<FORMAT, BITSTR, IFS, OFS>::preprocess(FORMAT2BIT& fm1
                                                    , FORMAT2BIT& fm2
                                                    , size_t pos)
{
    // adapter locates at 5'
    if (!_loc_tail)
    {
        fm1.seq.flip();
        fm2.seq.flip();
    }
}

template<template<bool, class> class FORMAT, class BITSTR, typename IFS, typename OFS>
void TaskProcessor<FORMAT, BITSTR, IFS, OFS>::detect_adapters()
{
    // detect adapter using the first N reads
    size_t total_lines = record_line * _detect_n_reads;
    std::vector<std::vector<std::string>> reads(2, std::vector<std::string>(total_lines));

    size_t i(0);
    for (; i < total_lines && _ifs[0].good() && _ifs[1].good(); ++i)
    {
        std::getline(_ifs[0], reads[0][i]);
        std::getline(_ifs[1], reads[1][i]);
    }
    // cutoff lines
    i = (i / record_line) * record_line;

    // reset ifstreams
    for (size_t j(0); j < 2; ++j)
    {
        if constexpr (std::is_same_v<remove_cvr_t<IFS>, BIO_filtering_istream>)
        {
            _ifs[j].pop();
            auto&& src(boost::iostreams::file_source(ifs_name[j], std::ios_base::binary));
            if (!src.is_open())
                throw std::runtime_error("Can't reopen input gz file normally\n");
            
            _ifs[j].push(src);
            if (!_ifs[j].good())
                throw std::runtime_error("Can't reopen input gz stream normally\n");
        }
        else
        {
            _ifs[j].close();
            _ifs[j].open(ifs_name[j]);
            if (!(_ifs[j].is_open() && _ifs[j].good()))
                throw std::runtime_error("Can't reopen input file normally\n");
        }
    }

    FORMAT2BIT fm1, fm2;
    // possible adapter fragments
    std::vector<std::vector<std::string>> adapter_frags(2);
    adapter_frags[0].reserve(i);
    adapter_frags[1].reserve(i);

    for (size_t j(0); j < i - 1; j += record_line) 
	{
        fm1 = FORMAT2BIT::parse_obj(
            reads[0].begin() + j);
        fm2 = FORMAT2BIT::parse_obj(
            reads[1].begin() + j);

        preprocess(fm1, fm2, j);

        _tr.cut_off_longer_seq(fm1, fm2);
        auto trim_pos = _tr.find_rc_pos(fm1, fm2, _adapters);
        if (trim_pos == fm1.seq.size())
            continue;

        adapter_frags[0].emplace_back(reads[0][j + 1].substr(trim_pos));
        adapter_frags[1].emplace_back(reads[1][j + 1].substr(trim_pos));
	}
    std::string tmp1, tmp2;
    bool is_low_complexity = false;
    if (_is_sensitive)
    {
        auto ret1 = assemble_adapters<true>(adapter_frags[0], init_kmer_size, 5);
        auto ret2 = assemble_adapters<true>(adapter_frags[1], init_kmer_size, 5);
        is_low_complexity = std::get<1>(ret1) || std::get<1>(ret2);
        tmp1 = std::get<0>(ret1);
        tmp2 = std::get<0>(ret2);
    }
    else
    {
        auto ret1 = assemble_adapters<false>(adapter_frags[0], init_kmer_size, 3);
        auto ret2 = assemble_adapters<false>(adapter_frags[1], init_kmer_size, 3);
        is_low_complexity = std::get<1>(ret1) || std::get<1>(ret2);
        tmp1 = std::get<0>(ret1);
        tmp2 = std::get<0>(ret2);
    }

    size_t min_len = 32;
    if (is_low_complexity)
    {
        min_len = std::min({tmp1.length(), tmp2.length(), (size_t)15});
    }
    else
    {
        min_len = std::min({tmp1.length(), tmp2.length(), (size_t)32});
    }
    
    _adapters[0] = tmp1.substr(0, min_len);
    _adapters[1] = tmp2.substr(0, min_len);
    
    std::cout << "detected adapter1: " << _adapters[0] << "\n";
    std::cout << "detected adapter2: " << _adapters[1] << "\n";

    if (_adapters[0] == "" || _adapters[1] == "")
    {
        std::cout << "unable to detect adapter, use default adapter\n";
        _adapters[0] = _default_adapters[0].substr(0, 32);
        _adapters[1] = _default_adapters[1].substr(0, 32);    
    }

}

template<template<bool, class> class FORMAT, class BITSTR, typename IFS, typename OFS>
template<class POOL>
bool TaskProcessor<FORMAT, BITSTR, IFS, OFS>::read_task(POOL& pool)
{
	bool eof = false;
	Task task;
    // write to buf until it is full.
	while((task.buf_idx = _buf_manager.inc_chunk_cnt()) != -1)
	{
		task.f_idx = _rw_count.rcount_fetch_add();
		eof = read_reads(task);
		
		if (eof)
		{
			_rw_count.rend_count = task.f_idx;
			pool.submit([this, task, &pool](){trim_task(task, pool);} );
			break;
		}

		pool.submit([this, task, &pool](){trim_task(task, pool);});
	}

	if (!eof)
	{
		std::this_thread::sleep_for(std::chrono::milliseconds(1));
	}

	return eof;
}

template<template<bool, class> class FORMAT, class BITSTR, typename IFS, typename OFS>
template<class POOL>
void TaskProcessor<FORMAT, BITSTR, IFS, OFS>::trim_task(Task task, POOL& pool)
{
	// trimming
	trim_reads(task);
	// add write task
	write_task(task);
	_buf_manager.dec_chunk_cnt(task.buf_idx);
	uint32_t wcount = _rw_count.wcount_fetch_add();
	
	if (_rw_count.rend_count == wcount)
	{
		pool.submit([this](){
			write_all_task();
			_rw_count.all_finished = true;
		});
	}
	else;
}

template<template<bool, class> class FORMAT, class BITSTR, typename IFS, typename OFS>
void TaskProcessor<FORMAT, BITSTR, IFS, OFS>::write_task(Task& task)
{
	std::ofstream f1(_tmp_dir[0] / fs::path(std::to_string(task.f_idx)), std::ios_base::binary);
	std::ofstream f2(_tmp_dir[1] / fs::path(std::to_string(task.f_idx)), std::ios_base::binary);

	if(!f1.good())
	    throw std::runtime_error("Can't read from tmp file1\n");

	if(!f2.good())
	    throw std::runtime_error("Can't read from tmp file2\n");

    f1.exceptions(std::ofstream::badbit | std::ofstream::failbit);
    f2.exceptions(std::ofstream::badbit | std::ofstream::failbit);

	size_t start_pos(task.buf_idx * _chunk_size);
	size_t end_pos(start_pos + _chunk_size);

	if (task.f_idx == _rw_count.rend_count)
	{
		end_pos = start_pos + _rw_count.rend_idx - 1;
	}
	else
	;

    std::string tmp1(""), tmp2("");
	tmp1.reserve(_buf_manager.buf[0][start_pos + 1].length() * 3 * _chunk_size);
	tmp2.reserve(_buf_manager.buf[0][start_pos + 1].length() * 3 * _chunk_size);

	for(size_t i(start_pos); i < end_pos; i += _record_line)
    {
         if (_buf_manager.buf[0][i].length() == 0) continue;
         for (size_t j(0); j < _record_line; ++j)
         {
             tmp1.append(_buf_manager.buf[0][i + j]);
             tmp1.append("\n");

             tmp2.append(_buf_manager.buf[1][i + j]);
             tmp2.append("\n");
         }
    }

    try
    {
        f1 << tmp1;
        f2 << tmp2;
    }
    catch(const std::exception& e)
    {
        std::cerr << "ERROR: " << e.what() << std::endl;
        std::cerr << "Failed to create/write temporary file, " << task.f_idx << "/" << _rw_count.rend_count << std::endl;
        std::abort();
    }
    
    f1.close();
    f2.close();
}

template<template<bool, class> class FORMAT, class BITSTR, typename IFS, typename OFS>
void TaskProcessor<FORMAT, BITSTR, IFS, OFS>::write_all_task()
{	
	// get all files under temporary directory 
	// that stores all the trimmed reads
	std::string tmp_fname;

    for (size_t j(0); j < 2; j++)
    {
        std::string path(_tmp_dir[j]);
        for (size_t i(0); i < _rw_count.rend_count + 1; i++)
        {
            tmp_fname = path / fs::path(std::to_string(i));
            std::ifstream tmp(tmp_fname, std::ios_base::binary);
            if(!tmp.good())
                throw std::runtime_error("Can't read from specified files\n");
    
            _ofs[j] << tmp.rdbuf();
            std::filesystem::remove(tmp_fname);
        }
        std::filesystem::remove_all(path);   
    }
}

template<template<bool, class> class FORMAT, class BITSTR, typename IFS, typename OFS>
void TaskProcessor<FORMAT, BITSTR, IFS, OFS>::trim_reads(Task& task)
{
	size_t start_pos(task.buf_idx * _chunk_size);
	size_t end_pos(start_pos + _chunk_size);
	if (task.f_idx == _rw_count.rend_count)
	{
		end_pos = start_pos + _rw_count.rend_idx - 1;
	}
	else
	;

    FORMAT2BIT fm1, fm2;

	for (size_t k(start_pos); k < end_pos; k += _record_line) 
	{
		fm1 = FORMAT2BIT::parse_obj(
            _buf_manager.buf[0].begin() + k);
        fm2 = FORMAT2BIT::parse_obj(
            _buf_manager.buf[1].begin() + k);

        preprocess(fm1, fm2, k);
        _tr.cut_off_longer_seq(fm1, fm2);
        
        auto trim_pos = _tr.template find_rc_pos<true>(fm1, fm2, _adapters);
        
        size_t seq_size(fm1.seq.size());
        if (trim_pos >= min_length)
        {
            fm1.to_iterator(_buf_manager.buf[0].begin() + k);
            fm2.to_iterator(_buf_manager.buf[1].begin() + k);
            
            if (trim_pos != seq_size)
            {
                fm1.trim(trim_pos);
                fm2.trim(trim_pos);

                for(size_t i(1); i < _record_line; i += 2)
                {
                    _buf_manager.buf[0][k + i].resize(
                                            trim_pos
                                            );

                    _buf_manager.buf[1][k + i].resize(
                                            trim_pos
                                            );
                }
             }
        }
        else
        {
            // trim_pos < min_length, abort reads
            _buf_manager.buf[0][k] = "";
            _buf_manager.buf[1][k] = "";
        }

	}
}

template<template<bool, class> class FORMAT, class BITSTR, typename IFS, typename OFS>
bool TaskProcessor<FORMAT, BITSTR, IFS, OFS>::read_reads(Task& task)
{
    // if (!_ifs[0].good())
	//     throw std::runtime_error("Can't read from input file1\n");

	// if (!_ifs[1].good())
	//     throw std::runtime_error("Can't read from input file2\n");
	size_t start_pos = task.buf_idx * _chunk_size;
	size_t i(start_pos);

	for (; i < start_pos + _chunk_size && _ifs[0].good() && _ifs[1].good(); ++i)
	{
		std::getline(_ifs[0], _buf_manager.buf[0][i]);
		std::getline(_ifs[1], _buf_manager.buf[1][i]);
	}

	bool eof1 = !_ifs[0].good();
	bool eof2 = !_ifs[1].good();
	
	if (eof1 != eof2)
	{
		throw std::runtime_error("Files aren't equal length.\n");
	}

	if (eof1)
	{
        // end buf index, and cutoff lines
		_rw_count.rend_idx = ((i / record_line) * record_line) - start_pos;
	}
	
	return eof1;
}

}