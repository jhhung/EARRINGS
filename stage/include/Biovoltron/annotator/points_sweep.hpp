/**
 * @file points_sweep.hpp
 * @brief An algorithm that finds overlapped annotations in databases. 
 * This algorithm is especially efficient when user has large queries 
 * and only query once.
 *
 * @author JHH corp
 */

#pragma once
#include <iostream>
#include <variant>
#include <vector>
#include <tuple>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <map>
#include <Biovoltron/annotator/exception.hpp>

namespace biovoltron::annotator
{
    
/**
 * @class Points_Sweep
 * @brief An algorithm that finds overlapped annotations between 
 * databases and queries. User has to specify the file format types 
 * that the databases and queries have. Note that the algorithm use 
 * "0-started, fully closed" coordinates.
 * 
 * The algorithm are as follows:<br>
 * 1. Store input annotations into databases.<br>
 * 2. Break each annotations into start and end points respectively.<br>
 * 3. Add query annotations and repeat step 2.<br>
 * 4. Annotate queries using databases and retrieve results.<br>
 * 
 * Member variables:<br>
 * private:<br>
 *      std::map<std::string, UInt> chrom_str_to_uint_<br>
 *      std::vector<std::vector<Point>> pts_pool_<br>
 *      std::vector<ret_value_type> ret_pool_<br>
 *      size_t chrom_count_<br>
 *      size_t cur_id_count_<br>
 *      size_t db_interval_num_<br>
 *
 */
template <class FIRST_INPUT_TYPE, class... REST_INPUT_TYPES>
class Points_Sweep
{
public:
    ///Assume chromosome length is less than uint32_t.
    using UInt = uint32_t;
    ///Assume point id range is less than int32_t.
    using Int = int32_t;
    ///ret_value_type depends on user-provided input types.
    using ret_value_type = std::tuple<
                                std::vector<FIRST_INPUT_TYPE *>
                              , std::vector<REST_INPUT_TYPES *>...
                                >;
    
    /**
     * @class Coverage
     * @brief A class which stores gene coverages.
     *
     * Member variable:<br>
     * std::string chrom<br>
     * UInt start<br>
     * UInt end<br>
     * size_t overlaps_num<br>
     *
     */
    struct Coverage
    {
        std::string chrom;
        UInt start;
        UInt end;
        size_t overlaps_num;
    };
    
    ///Empty constructor which does nothing.
    Points_Sweep() {}

    /**
     * @brief An interface that allows user to set multiple databases.
     * set_db() can only be called once, otherwise, the program will 
     * throw exceptions. reset() must be done if the user wants to 
     * change databases.
     *
     * @param dbs One or more vectors of objects. Objects in 
     * the vector must provide functions: "std::string get_chrom()", 
     * "Uint get_start()" and "Uint get_end()".
     */
    template <class... DBs>
    void set_db(DBs &&... dbs) 
    {
        // reject setting db if db is not empty
        if (db_interval_num_ != 0) 
        {
            throw Overlap_Annotator(
                "Error: Cannot re-add db. \
                \"Reset\" must be done before adding.\n"
            );
        } 
        else 
        {
            add_pts(std::forward<DBs>(dbs)...);

            // reserve db
            for (const auto &pool : pts_pool_)
            {
                db_interval_num_ += pool.size();
            }
            db_interval_num_ = (db_interval_num_ >> 1);
        }
    }

    /**
     * @brief An interface that allows user to set multiple queries.
     * set_query() can only be called once, otherwise, the program will 
     * throw exceptions. reset() must be done if the user wants to 
     * change queries.
     *
     * @param queries One or more vectors of objects. Objects in 
     * the vector must provide functions: "std::string get_chrom()", 
     * "Uint get_start()" and "Uint get_end()".
     */
    template <class... QUERIES>
    void set_query(QUERIES &&... queries) 
    {
        // reject setting query if query is not empty
        if (ret_pool_.size() != 0) 
        {
            throw Overlap_Annotator(
                "Error: Cannot re-add query. \
                \"Reset\" must be done before adding.\n"
            );
        } 
        else 
        {
            if (db_interval_num_ == 0) 
            {
                throw Overlap_Annotator(
                    "Error: Please add db first before adding query.\n"
                );
            } 
            else
            {
                add_pts(std::forward<QUERIES>(queries)...);

                // initialize ret pool
                UInt total_query = 0;
                for (size_t i = 0; i < chrom_count_; ++i)
                {
                    total_query += (pts_pool_[i].size() >> 1);
                }
                size_t reserve_size = total_query - db_interval_num_;
                ret_pool_.resize(reserve_size);
                coverage_pool_.resize(reserve_size);
            }
        }
    }

    /**
     * @brief Reset databases and queries.
     * 
     */
    void reset() 
    {
        // reset db & query
        decltype(pts_pool_)().swap(pts_pool_);
        decltype(ret_pool_)().swap(ret_pool_);
        decltype(coverage_pool_)().swap(coverage_pool_);
        decltype(chrom_str_to_uint_)().swap(chrom_str_to_uint_);
        db_interval_num_ = 0;
        chrom_count_ = 0;
        cur_id_count_ = 1;
    }

    /**
     * @brief Get overlapped annotations.
     *
     * @return A reference vector that stores the overlapped annotations.
     */
    std::vector<ret_value_type>& annotate()
    {
        if (check_db_query())
        {
            sort_pts();
            find_overlap();
        }
        else
        ;  // check fail, check_db_query() will throw exception.

        return ret_pool_;
    }

    ///Get gene coverages.
    std::vector<std::vector<Coverage>>& get_coverage()
    {
        if (check_db_query())
        {
            sort_pts();
            cal_coverage();
        }
        else
        ;  // check fail, check_db_query() will throw exception.

        return coverage_pool_;
    }

    ~Points_Sweep(){}

private:
    ///Store the pointer pointing to original annotation.
    using variant_type = std::variant<
                            FIRST_INPUT_TYPE *
                          , REST_INPUT_TYPES *...
                            >;
    
    /**
     * @class Point
     * @brief A class which stores the splitted annotation.
     *
     * Member variable:<br>
     * private:<br>
     * variant_type variant_<br>
     * UInt pos_<br>
     * Int id_<br>
     *
     * This class provides several functions to give access to 
     * original annotation's pointer, position (either start or end 
     * position) and an unique id.
     */
    class Point
    {
    private:
        ///Pointing to original annotation.
        variant_type variant_;
        ///Storing original annotation's start or end position.
        UInt pos_;
        ///Unique id to identify each annotation.
        Int id_;

    public:
         /**
         * @brief Constructor of class Point.
         * 
         * @param var: variant that stores the pointer pointing to 
         *        the original annotation.
         *        pos: start or end position of the annotation.
         *        id: the unique id of the annotation.
         */
        template<class VAR, class POS, class ID>
        Point(VAR&& var, POS&& pos, ID&& id)
            :variant_(var), pos_(pos), id_(id){}
        
        ///Get original annotation pointer.
        variant_type get_var() const { return variant_; }

        /**
         * @brief Get annotation position, either start or end position. 
         * Start or end position is decided by the sign of id_.
         *
         * @return An UInt that represents the position.
         */
        UInt get_pos() const { return pos_; }

        /**
         * @brief Get unique id of an annotation. Positive value 
         * indicates start position, while negative value indicates 
         * end position.
         *
         * @return A Int id number.
         */
        Int get_id() const { return id_; }

        /**
         * @brief Two points are originated from the same annotation 
         * if they add up to 0.
         *
         * @return Two points are from the same annotation or not.
         */
        bool operator==(const Point &pt) const 
        {
            return (this->get_id() + pt.get_id()) == 0;
        }
    };

    ///Custom hash function for an unordered_multiset.
    struct pts_hash
    {
		size_t operator()(const Point &pt) const
		{ 
            return std::hash<UInt>()(abs(pt.get_id()));
        }
	};

    ///Store the mapping of chromosome string to UInt.
    std::map<std::string, UInt> chrom_str_to_uint_;
    /**
     * @brief Storing all internal data structure called "Point". 
     * This is stored according to chromosome.
     *
     */
    std::vector<std::vector<Point>> pts_pool_;
    ///Store the overlapped annotations.
    std::vector<ret_value_type> ret_pool_;
    ///Store gene coverages.
    std::vector<std::vector<Coverage>> coverage_pool_;

    ///Record current mapping count of chrom_str_to_uint_;
    size_t chrom_count_ = 0;
    ///Record current id number used to construct "Point".
    size_t cur_id_count_ = 1;
    ///Record total number of annotations in database.
    size_t db_interval_num_ = 0;

    /**
     * @brief Make sure that databases and queries are not empty 
     * before finding overlaps or calculating coverages.
     * 
     */
    bool check_db_query()
    {
        if (db_interval_num_ == 0) 
        {
            throw Overlap_Annotator(
                "Error: Please add db before annotate().\n"
            );
        } 
        else
        {
            if (ret_pool_.size() == 0) 
            {
                throw Overlap_Annotator(
                    "Error: Please add query before annotate().\n"
                );
            } 
            else
            {
                return true;
            }
        }

        return false;
    }

    ///Sort all the points stored in pts_pool_.
    void sort_pts()
    {
        // sort pts
        for (auto &pool : pts_pool_)
        {
            std::sort(pool.begin(), pool.end(),
                [](const Point &pts1, const Point &pts2) 
                    {
                        if (pts1.get_pos() == pts2.get_pos()) 
                        {
                            // start must precedes end.
                            return pts1.get_id()
                                    > pts2.get_id();
                        }
                        else
                        {
                            return pts1.get_pos() 
                                    < pts2.get_pos();
                        }
                    }
                );
        }
    }

    ///End of variadic parameter unpacking.
    void add_pts() {}

    /**
     * @brief Transform input annotations into internal data structure 
     * "Point" and stores them into pts_pool_.
     *
     * @param input All inputs must be stored into a std::vector.
     * 
     */
    template <class INPUT, class... INPUTS>
    void add_pts(INPUT &&input, INPUTS &&... inputs)
    {
        size_t chrom;
        UInt id;
        
        for (auto &element : input)
        {
            if (chrom_str_to_uint_.find(element.get_chrom()) 
                                != chrom_str_to_uint_.end()) 
            {
                chrom = chrom_str_to_uint_[element.get_chrom()];
            } 
            else 
            {
                chrom = chrom_count_;
                chrom_str_to_uint_[element.get_chrom()] = chrom_count_++;
                pts_pool_.emplace_back(std::vector<Point>());
            }
            id = cur_id_count_++;  // id starts from 1
            Point start_pts(&element, element.get_start(), id);
            Point end_pts(&element, element.get_end(), -id);
            pts_pool_[chrom].emplace_back(start_pts);
            pts_pool_[chrom].emplace_back(end_pts);
        }
        
        add_pts(std::forward<INPUTS>(inputs)...);
    }

    /**
     * @brief Find overlapped annotations between databases and 
     * queries.
     * 
     * The algorithm is as follows:<br>
     * 1. Two empty unordered_multiset, cache and q_cache, will be 
     * created to store points. The point's id is used to identify 
     * if two points are originated from the same annotation.<br>
     * 2. Iterate through all the points in pts_pools_. Each point 
     * will be inserted into a cache/q_cache if it is a start point 
     * (id > 0), or be deleted from cache/q_cache if it is an end 
     * point (id < 0). If a point belongs to an query, then it will be 
     * inserted/deleted to/from q_cache, otherwise, cache.<br>
     * 3. Store overlapped points if there exists.
     *
     */
    void find_overlap()
    {
        std::unordered_multiset<Point, pts_hash> cache, q_cache;
        
        bool is_query = false;

        for (size_t i = 0; i < chrom_count_; ++i)
        {
            for (const auto &pt : pts_pool_[i])
            {
                is_query = size_t(abs(pt.get_id())) > db_interval_num_;
                UInt condition = (pt.get_id() < 0) + (is_query << 1);
                
                // switch is used here to replace if/else statement 
                // in order to speed up the whole process.
                switch(condition) 
                {
                    // pt is start but not a query
                    case 0:
                        cache.insert(pt);
                        break;
                    
                    // pt is end but not a query
                    case 1:
                        for (auto &q : q_cache) 
                            store_overlaps(q, pt);
                        
                        cache.erase(pt);
                        break;

                    // pt is start and also a query
                    case 2:
                        q_cache.insert(pt);
                        break;

                    // pt is end and also a query
                    case 3:
                        for (auto &p : cache) 
                            store_overlaps(pt, p);   
                            
                        q_cache.erase(pt);
                        break;
                    // default exception
                    default:
                        throw Overlap_Annotator(
                            "Errors occure when finding overlaps.\n"
                        );
                        break;
                }
            }
            cache.clear();
            q_cache.clear();
        }
    }

    /**
     * @brief Store the pointer of the original annotation of the 
     * overlapped point using std::visit to deduce the type stored 
     * in variant.
     *
     */
    void store_overlaps(const Point &query_pts, const Point &overlap_pts)
    {
        UInt query_idx = size_t(abs(query_pts.get_id())-1) - db_interval_num_;

        std::visit([this, &query_idx, &overlap_pts](auto&& arg) 
        {
            using Q = std::decay_t<decltype(arg)>;
            store_overlap_by_type<Q>(query_idx, overlap_pts);
        }, overlap_pts.get_var());
    }

    ///Deduce variant type and store the overlapped annotations.
    template<class TYPE>
    void store_overlap_by_type(
        const UInt &query_idx
      , const Point &overlap_pts
      ) 
    {
        std::get<std::vector<TYPE>>(ret_pool_[query_idx])
                        .emplace_back(
                            std::get<TYPE>(overlap_pts.get_var())
                            );
    }


    /**
     * @brief Calculate gene coverages.
     * 
     * The algorithm is similar to find_overlap(). Coverages only 
     * change when inserting or deleting points.
     *
     */
    void cal_coverage()
    {
        std::unordered_multiset<Point, pts_hash> cache, q_cache;
        
        bool is_query = false;

        for (size_t i = 0; i < chrom_count_; ++i)
        {
            for (const auto &pt : pts_pool_[i])
            {
                is_query = size_t(abs(pt.get_id())) > db_interval_num_;
                UInt condition = (pt.get_id() < 0) + (is_query << 1);
                
                switch(condition) 
                {
                    // coverages change
                    case 0:
                        cache.insert(pt);
                        for (auto &q : q_cache) 
                            store_coverage(q, pt.get_pos(), cache.size(), 0);
                        break;
                    
                    // coverages change
                    case 1:
                        cache.erase(pt);
                        for (auto &q : q_cache) 
                            store_coverage(q, pt.get_pos(), cache.size(), 0);
                        break;

                    // start recording coverages
                    case 2:
                        q_cache.insert(pt);
                        store_coverage(pt, pt.get_pos(), cache.size(), 1);
                        break;

                    // stop recording coverages
                    case 3:
                        q_cache.erase(pt);
                        store_coverage(pt, pt.get_pos(), cache.size(), 2);
                        break;
                    // default exception
                    default:
                        throw Overlap_Annotator(
                            "Errors occure when calculating coverages.\n"
                        );
                        break;
                }
            }
            cache.clear();
            q_cache.clear();
        }
    }

    ///Store coverages when inserting or deleting points.
    void store_coverage(const Point &query_pts 
                      , const UInt &pos
                      , const size_t &overlaps_num
                      , const size_t &state
                       )
    {
        UInt query_idx = size_t(abs(query_pts.get_id())-1) - db_interval_num_;
        Coverage coverage;

        switch(state)
        {
            case 0:
                coverage_pool_[query_idx][coverage_pool_[
                                            query_idx
                                            ].size() - 1
                                         ].end
                                         = pos;
            case 1:
                coverage.end = coverage.start = pos;
                std::visit([&coverage, &query_pts](auto&& arg) 
                {
                    using Q = std::decay_t<decltype(arg)>;
                    coverage.chrom = std::get<Q>(query_pts.get_var())->get_chrom();
                }, query_pts.get_var());
                coverage.overlaps_num = overlaps_num;
                coverage_pool_[query_idx].emplace_back(coverage);
                break;
            case 2:
                coverage_pool_[query_idx][coverage_pool_[
                                            query_idx
                                            ].size() - 1
                                         ].end
                                         = pos;
                break;
            default:
                break;
        }
    }
};
} // namespace biovoltron::annotator
