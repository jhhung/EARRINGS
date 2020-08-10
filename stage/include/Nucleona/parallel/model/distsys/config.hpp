#pragma once
#include <utility>
namespace nucleona{ namespace parallel{ namespace model{ namespace distsys{
struct Config
{
	const static std::size_t worker_buffer_size = 2000;
	const static std::size_t boost_join_time_out = 300;//millisecond.
};
}}}}
