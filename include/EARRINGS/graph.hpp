#ifndef KYUTORA_GRAPH_HPP
#define KYUTORA_GRAPH_HPP

#pragma once

#include <string_view>
#include <limits>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <EARRINGS/common.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/filtered_graph.hpp>

namespace kyutora
{

static auto PRUNE_FACTOR = 50;

struct VertexProperty {
    std::string_view kmer;
};

struct EdgeProperty {
    std::size_t count = 0;
    double score = std::numeric_limits<double>::lowest();
};

/*
 * adjacency lis:
 * OutEdgeList: vecS (std::vector)
 * VertexList: vecS 
 *
 */
using Graph = boost::adjacency_list<
    boost::vecS, boost::vecS, boost::bidirectionalS,
    VertexProperty, EdgeProperty>;
using Vertex = boost::graph_traits<Graph>::vertex_descriptor;
using Edge = boost::graph_traits<Graph>::edge_descriptor;
using Path = std::vector<Vertex>;

struct cycle_detector : public boost::dfs_visitor<>
{
    cycle_detector(bool& has_cycle) : _has_cycle(has_cycle) {}
    template <class Edge, class Graph>
    void back_edge(Edge, Graph&) {
        _has_cycle = true;
    }
protected:
    bool& _has_cycle;
};


class GraphWrapper
{
private:
    size_t kmer_size;
public:
    static constexpr auto NUM_PATHS = 32;
    //static constexpr auto KMER_SIZE = 10;
    
    static auto get_dup_kmers(std::string_view seq, std::size_t size)
    {
        std::set<std::string_view> all_kmers, dup_kmers;
        for (std::size_t i = 0; i <= seq.size()-size; i++)
        {
            auto kmer = seq.substr(i, size);
            if (auto [iter, success] = all_kmers.insert(kmer); !success)
                dup_kmers.insert(kmer);
        }
        return dup_kmers;
    }

    GraphWrapper(size_t kmer_size = init_kmer_size): kmer_size(kmer_size) {}
    size_t get_kmer_size() { return kmer_size; }
// private:
    Graph g;
    std::vector<Path> paths;
    struct EdgeFilter
    {
        bool operator()(Edge e) const
        { return (*g)[e].count >= PRUNE_FACTOR; }
        bool operator()(Vertex v) const
        {
            for (auto e : boost::make_iterator_range(boost::out_edges(v, *g)))
                if ((*g)[e].count >= PRUNE_FACTOR)
                    return true;
            for (auto e : boost::make_iterator_range(boost::in_edges(v, *g)))
                if ((*g)[e].count >= PRUNE_FACTOR)
                    return true;
            return false;
        }
        Graph* g;
    } filter{&g};

    boost::filtered_graph<Graph, EdgeFilter, EdgeFilter> fg = {g, filter, filter};
    std::vector<Vertex> sources, sinks;
    std::set<std::string_view> dup_kmers;
    std::map<std::string_view, Vertex> unique_kmers;
// private:
    auto is_source(Vertex v) const
    {
        if (boost::in_degree(v, fg) == 0)
            for (auto e : boost::make_iterator_range(boost::out_edges(v, fg)))
                if (g[e].count >= PRUNE_FACTOR)
                    return true;
        return false;
    }

    auto is_sink(Vertex v) const
    {
        if (boost::out_degree(v, fg) == 0)
            for (auto e : boost::make_iterator_range(boost::in_edges(v, fg)))
                if (g[e].count >= PRUNE_FACTOR)
                    return true;
        return false;
    }

    auto find_sources_and_sinks()
    {
        for (auto v : boost::make_iterator_range(boost::vertices(fg)))
        {
            if (is_source(v)) sources.push_back(v);
            if (is_sink(v)) sinks.push_back(v);
        }
    }

    auto create_edge(Vertex u, Vertex v)
    {
        auto [e, success] = boost::add_edge(u, v, g);
        g[e].count++;
    }

    auto create_vertex(std::string_view kmer)
    {
        auto v = boost::add_vertex(g);
        g[v].kmer = kmer; // add new vertex
        if (dup_kmers.find(kmer) == dup_kmers.end())
            unique_kmers.emplace(kmer, v);
        return v;
    }

    // create new vertex if not exist
    auto get_vertex(std::string_view kmer)
    {
        if (auto it = unique_kmers.find(kmer); it != unique_kmers.end())
            return it->second;  // return vertex given string
        return create_vertex(kmer);
    }

    auto extend_chain(Vertex u, std::string_view kmer)
    {
        // get iterator's out edges of a given vertices
        for (auto e : boost::make_iterator_range(boost::out_edges(u, g)))
        {
            // given an edge, get the targeting vertex that the edge is pointing to 
            auto v = boost::target(e, g);
            // every vertex contains only one character (except for the source vertex)
            if (g[v].kmer.back() == kmer.back())
            {
                g[e].count++;
                return v;
            }
        }

        auto v = get_vertex(kmer);
        create_edge(u, v);
        return v;
    }

    void add_seq(std::string_view seq)
    {
        // create source vertex if not exist
        auto v = get_vertex(seq.substr(0, kmer_size));
        // for every k-mer, extending it
        for (auto i = 1u; i <= seq.size()-kmer_size; i++)
            v = extend_chain(v, seq.substr(i, kmer_size));
    }

    void path_finder(Vertex from, Vertex to, Path& path)
    {
        path.push_back(from);
        if (from == to)
            paths.push_back(path);
        else
        {
            for (auto e : boost::make_iterator_range(boost::out_edges(from, fg)))
            {
                auto v = target(e, fg);
                if (std::find(path.begin(), path.end(), v) == path.end())
                    path_finder(v, to, path);
            }
        }
        path.pop_back();
    }

    void find_paths()
    {
        for (auto &i : sources)
        {
            for (auto &j : sinks)
            {
                Path path;
                path_finder(i, j, path);
            }
        }
    }

    auto get_adapters()
    {
        std::vector<std::string> adapters;
        // sort according to the count of source vertex
        std::sort(paths.begin(), paths.end(),
                [this](const auto& lhs, const auto& rhs)
                {
                    // lhs[0]/rhs[0] -> source
                    for (auto lhs_e : boost::make_iterator_range(boost::out_edges(lhs[0], g)))
                    {
                        for (auto rhs_e : boost::make_iterator_range(boost::out_edges(rhs[0], g)))
                        {
                            return g[lhs_e].count > g[rhs_e].count;
                        }
                    }

                    //return lhs.size() > rhs.size();
                }); 
        
        for (const auto& path : paths)
        {
            auto u = path[0];
            auto seq = std::string{g[u].kmer.data(), g[u].kmer.size()};
            for (std::size_t i = 1; i < path.size(); i++)
            {
                auto v = path[i];
                seq += g[v].kmer.back();
                u = v;
            }
            adapters.emplace_back(std::move(seq));
        }

        return adapters;
    }


public:
    // build graph from tails
    void build(const std::vector<std::string>& seqs)
    {
        // get non-unique kmers from all tails
        for (const auto& seq : seqs)
            for (auto kmer : get_dup_kmers(seq, kmer_size))
                dup_kmers.insert(kmer);

        for (const auto& seq : seqs)
            add_seq(seq);
    }

    bool has_cycles() const
    {
        bool has_cycle = false;
        cycle_detector vis(has_cycle);
        boost::depth_first_search(g, boost::visitor(vis));
        return has_cycle;
    }

    void print() const
    {
        std::ofstream os("graph.dot");
        os << "digraph assembly_graphs {";
        for (auto e : boost::make_iterator_range(boost::edges(fg)))
        {
            os << boost::source(e, fg) << " -> " << boost::target(e, fg) << " ";
            auto count = fg[e].count;
            if (count < PRUNE_FACTOR)
                 os << "[label=" << count << ",style=dotted,color=grey];\n";
            else os << "[label=" << count << "];\n";
        }

        for (auto v : boost::make_iterator_range(boost::vertices(fg)))
        {
            os << v << " ";
            auto kmer = fg[v].kmer;
            if (boost::in_degree(v, fg) == 0)
                 os << "[label=" << kmer << ",shape=box]\n";
            else os << "[label=" << kmer.back() << ",shape=box]\n";
        }
        os << "}";
    }
};

} // kyutora

#endif //KYUTORA_GRAPH_HPP
