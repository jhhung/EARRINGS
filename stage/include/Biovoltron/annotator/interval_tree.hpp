/**
 * @file points_sweep.hpp
 * @brief An algorithm that finds overlapped annotations in databases. 
 * This algorithm is especially efficient when user has small queries 
 * and needs to query multiple times.
 *
 * @author JHH corp
 */

#pragma once
#include <iostream>
#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <algorithm>
#include <variant>
#include <typeinfo>
#include <Biovoltron/annotator/exception.hpp>

namespace biovoltron::annotator
{

/**
 * @class Interval_Tree
 * @brief An algorithm that finds overlapped annotations in databases. 
 * User has to specify the file format types that the databases and 
 * queries have.
 * 
 * The algorithm is as follows:<br>
 * 1. Store input annotations into databases.<br>
 * 2. Transform all the annotations into internal data structure 
 * called Interval and sort the intervals according to their end 
 * positions to construct interval trees.<br>
 * 3. Find overlap annotations using the interval trees.<br>
 * 
 * Member variables:<br>
 *      std::vector<std::vector<Interval> > db_interval_<br>
 *      std::vector<ITNode*> roots_<br>
 *      std::vector<ret_value_type> overlaps_<br>
 *      std::map<std::string, UInt> chrom_str_to_uint_<br>
 *      size_t chrom_count_<br>
 *
 */
template <class FIRST_INPUT_TYPE, class... REST_INPUT_TYPES>
class Interval_Tree
{
public:
    ///Assume chromosome length is less than uint32_t.
    using UInt = uint32_t;
    ///ret_value_type depends on user-provided input types.
    using ret_value_type = std::tuple<
                                std::vector<FIRST_INPUT_TYPE *>
                              , std::vector<REST_INPUT_TYPES *>...
                                >;

    ///Empty constructor which does nothing.
    Interval_Tree() {}

    /**
     * @brief An interface that allows user to set multiple databases.
     * set_db() can only be called once, otherwise, the program will 
     * throw exceptions. reset() must be done if user wants to change 
     * databases.
     *
     * @param dbs One or more vectors of objects. Objects in 
     * the vector must provide functions: "std::string get_chrom()", 
     * "Uint get_start()" and "Uint get_end()".
     */
    template <class... DBs>
    void set_db(DBs &&... dbs) 
    {
        if (db_interval_.size() != 0)
        {
            throw Overlap_Annotator(
                "Interval_Tree Error: Cannot re-add db."
                " \"Reset\" must be done before adding.\n"
            );
        }
        else
        {
            construct_db(std::forward<DBs>(dbs)...);

            UInt med_idx = 0;
            for (auto &v : db_interval_)
            {
                med_idx = v.size() - 1;
                med_idx = ((med_idx & 1) == 0) 
                            ? UInt(med_idx >> 1) 
                            : UInt(med_idx >> 1) + 1;
                
                std::vector<bool> is_used(v.size(), false);
                std::sort(v.begin(), v.end());

                roots_.emplace_back(
                    construct_tree(v, is_used, med_idx, 0, v.size()));
            }
        }
    }
    
    // end of variadic
    void set_query() {}

    /**
     * @brief An interface that allows user to set multiple queries.
     * set_query() can only be called once, otherwise, the program will 
     * throw exceptions. reset() or reset_query() must be done if user 
     * wants to change queries.
     *
     * @param queries One or more vectors of objects. Objects in 
     * the vector must provide functions: "std::string get_chrom()", 
     * "Uint get_start()" and "Uint get_end()".
     */
    template <class... QUERIES>
    void set_query(QUERIES &&... queries) 
    {
        if (overlaps_.size() != 0)
        {
            throw Overlap_Annotator(
                "Interval_Tree Error: Cannot re-add query."
                " \"Reset\" must be done before adding.\n"
            );
        }
        else
        {
            if (db_interval_.size() == 0)
            {
                throw Overlap_Annotator(
                    "Interval_Tree Error:"
                    "Please add db first before adding query.\n"
                );
            }
            else
            {
                construct_ret(std::forward<QUERIES>(queries)...);   
            }
        }
    }

    /**
     * @brief Reset databases and queries.
     * 
     */
    void reset() 
    {
        decltype(db_interval_)().swap(db_interval_);
        decltype(roots_)().swap(roots_);
        decltype(overlaps_)().swap(overlaps_);
        decltype(chrom_str_to_uint_)().swap(chrom_str_to_uint_);
        chrom_count_ = 0;
    }

    /**
     * @brief Reset only queries.
     * 
     */
    void reset_query()
    {
        decltype(overlaps_)().swap(overlaps_);
    }

    /**
     * @brief Get overlapped annotations.
     *
     * @return A reference vector that stores the overlapped annotations.
     */
    std::vector<ret_value_type>& annotate() 
    {
        if (db_interval_.size() == 0)
        {
            throw Overlap_Annotator(
                "Interval_Tree Error: Please add db before annotate().\n"
            );
        }
        else 
        {
            if (overlaps_.size() == 0)
            {
                throw Overlap_Annotator(
                    "Interval_Tree: Please add query before annotate().\n"
                );
            }
            else
            ;  // do nothing.
        }
        return overlaps_;
    }

private:
    ///Store the pointer pointing to the original annotation.
    using variant_type = std::variant<
                            FIRST_INPUT_TYPE *
                          , REST_INPUT_TYPES *...
                            >;

    /**
     * @class Interval
     * @brief A class which turns annotation into interval.
     *
     * Member variable:<br>
     * private:<br>
     *      variant_type data_ptr<br>
     *
     * This class provides several functions to give access to 
     * original interval's pointer, start and end position.
     */
    class Interval
    {
    private:
        // variant which stores data pointer
        // eg. variant<GFF*, VCF*, WIG*> w = &gff;
        variant_type data_ptr;

    public:
        /**
         * @brief class Interval's constructor which takes pointer to 
         * annotation object as input.
         *
         * @param ptr pointer to annotation object.
         */
        Interval(variant_type&& ptr): data_ptr(ptr) { }
        
        ///Return pointer to original annotation object.
        variant_type get_data_ptr() const { return data_ptr; }

        ///Get start position of a pointer to an annotation object.
        static UInt get_start(const variant_type& var)
        {
            UInt start;
            std::visit([&var, &start](auto&& arg) 
            {
                using Q = std::decay_t<decltype(arg)>;
                start = std::get<Q>(var)->get_start();
            }, var);
            
            return start;
        }

        ///Get end position of a pointer to an annotation object.
        static UInt get_end(const variant_type& var)
        {
            UInt end;
            std::visit([&var, &end](auto&& arg) 
            {
                using Q = std::decay_t<decltype(arg)>;
                end = std::get<Q>(var)->get_end();
            }, var);
            
            return end;
        }

        ///Compare the end position of two interval objects.
        bool operator<(const Interval &interval) const
        {
            return Interval::get_end(this->get_data_ptr()) 
                    < Interval::get_end(interval.get_data_ptr());
        }
    };

    /**
     * @class INode
     * @brief A tree's node class.
     *
     * Member variable:<br>
     *      UInt median<br>
     *      std::vector<Interval> vec_interval<br>
     *      ITNode *left<br>
     *      ITNode *right<br>
     *
     */
    struct ITNode
    {
        ///Store the median value of the end position of intervals.
        UInt median;
        ///Store all the interval that crosses median.
        std::vector<Interval> vec_interval;
        ///Pointer that points to the left node.
        ITNode *left;
        ///Pointer that points to the right node.
        ITNode *right;
    };

    ///Store all the input intervals.
    std::vector<std::vector<Interval> > db_interval_;
    ///Store all the roots of the interval trees.
    std::vector<ITNode*> roots_;
    ///Store the overlapped annotations.
    std::vector<ret_value_type> overlaps_;
    ///Store the mapping of chromosome string to UInt.
    std::map<std::string, UInt> chrom_str_to_uint_;
    ///Record current mapping count of chrom_str_to_uint_;
    size_t chrom_count_ = 0;

    ///End of variadic parameter unpacking.
    void construct_db() {}

    /**
     * @brief Transform input annotations into internal data structure 
     * "Interval" and stores them in db_interval_.
     *
     * @param input All inputs must be stored into a std::vector.
     * 
     */
    template<class INPUT, class... INPUTS>
    void construct_db(INPUT &&input, INPUTS &&... inputs) 
    {
        UInt chrom;

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
                db_interval_.emplace_back(std::vector<Interval>());
            }
            
            Interval interval(std::move(&element));
            db_interval_[chrom].emplace_back(interval);
        }

        construct_db(std::forward<INPUTS>(inputs)...);
    }

    ///End of variadic parameter unpacking.
    void construct_ret() {}

    /**
     * @brief Transform input annotations into internal data structure 
     * "Interval" and find interval overlaps using db_interval_.
     *
     * @param input All inputs must be stored into a std::vector.
     * 
     */
    template<class INPUT, class... INPUTS>
    void construct_ret(INPUT &&input, INPUTS &&... inputs) 
    {
        UInt chrom;
        ret_value_type ret;
        for (auto& element : input) 
        {
            if (chrom_str_to_uint_.find(element.get_chrom()) 
                        != chrom_str_to_uint_.end())
            {
                chrom = chrom_str_to_uint_[element.get_chrom()];
                Interval interval(std::move(&element));
                find_overlap(roots_[chrom], interval, ret);
            } 
            else 
            ;
            
            overlaps_.emplace_back(ret);
            decltype(ret)().swap(ret);
        }
        
        construct_ret(std::forward<INPUTS>(inputs)...); 
    }

    /**
     * @brief Calculate the median index of bool vector that is not 
     * marked as true.
     *
     * @param is_used A vector which indicates if an interval is used.
     *        unused_count The number of unused intervals between 
     *                     begin_idx and end_idx.
     *        begin_idx The index that the search begins.
     *        end_idx The index that the search ends.
     */
    UInt cal_med_idx(
        std::vector<bool> &is_used
      , size_t &unused_count
      , const UInt &begin_idx
      , const UInt &end_idx
    )
    {
        if ((unused_count & 1) == 0) 
        { 
            unused_count = size_t(unused_count >> 1); 
        }
        else 
        {
            unused_count = size_t(unused_count >> 1) + 1; 
        }
        
        size_t cur_count = 0;
        for (size_t i = begin_idx; i < end_idx; ++i ) 
        {
            if (!is_used[i]) cur_count++;
            else
            ;

            if (cur_count == unused_count) return i;
            else
            ;
        }

        return end_idx;
    }

    /**
     * @brief Construct the interval trees. 
     * The algorithm is as follows:
     * 1. Find the median of the intervals' end position that are not 
     * marked as used.<br>
     * 2. The intervals which crosses the median are stores in current 
     * ITNode. Also, mark those intervals as used.<br>
     * 3. For intervals whose start position are smaller than the 
     * median value of current ITNode goes to ITNode->left.<br>
     * 4. For intervals whose end position are greater than the median 
     * value of current ITNode goes to ITNode->right.<br>
     *
     * @param intervals The sorted intervals used to construct the trees.
     *        is_used A vector which indicates if an interval is used.
     *        unused_count The number of unused intervals between 
     *                     begin_idx and end_idx.
     *        begin_idx The index that the search begins.
     *        end_idx The index that the search ends.
     */
    ITNode* construct_tree(
        std::vector<Interval> &intervals
      , std::vector<bool> &is_used
      , const UInt& med_idx
      , const UInt& begin_idx
      , const UInt& end_idx
    )
    {
        ITNode* node = new ITNode();
        node->median = Interval::get_end(
                                intervals[med_idx].get_data_ptr());
        // use to count the idx of median
        size_t unused_count_left = 0, unused_count_right = 0;

        variant_type var;
        // Intervals which pass the median are kept in the node
        for (size_t i = begin_idx; i < end_idx; ++i) 
        {
            if (!is_used[i])
            {
                var = intervals[i].get_data_ptr();
                if ((Interval::get_start(var) <= node->median) && 
                      (Interval::get_end(var) >= node->median)) 
                {
                    is_used[i] = true;
                    node->vec_interval.emplace_back(
                                            std::move(intervals[i]));
                }
                else 
                {
                    (i < med_idx) ? unused_count_left++ 
                                  : unused_count_right++;
                }
            }
            else 
            ;  // the interval is used in previous node, 
               //or the interval does not cross with the median
            
        }

        // calculate new median
        // left
        if (unused_count_left != 0) 
        {
            size_t med_idx_left = cal_med_idx(is_used
                                            , unused_count_left
                                            , begin_idx
                                            , med_idx
                                            );
            node->left = construct_tree(intervals
                                      , is_used
                                      , med_idx_left
                                      , begin_idx
                                      , med_idx
                                      );
        }
        else 
        {
            node->left = nullptr;
        }
        
        // right
        if (unused_count_right != 0)
        {
            size_t med_idx_right = cal_med_idx(is_used
                                             , unused_count_right
                                             , size_t(med_idx + 1)
                                             , end_idx
                                             );
            node->right = construct_tree(intervals
                                       , is_used
                                       , med_idx_right
                                       , size_t(med_idx + 1)
                                       , end_idx
                                       );
        }
        else 
        {
            node->right = nullptr;
        }
        
        return node;
    }


    /**
     * @brief Find overlaps between databases and queries.
     * The algorithm is as follows:<br>
     * 1. If the query interval crosses current node's median value, 
     * search overlaps for the node's vec_interval. Then, traverse 
     * to both left and right nodes.<br>
     * 2. If the end position of the query interval is smaller than 
     * current node's median value, then search overlaps for the 
     * node's vec_interval. Then, traverse only to the left node.<br>
     * 3. If the start position of the query interval is greater than 
     * current node's median value, then search overlaps for the 
     * node's vec_interval. Then, traverse only to the right node.<br>
     *
     * @param node Node which is currently traversing.
     *        query_interval The interval that the user wishes to 
     *                       find overlaps.
     *        ret Pointers which point to overlapped annotations.
     */
    void find_overlap(
        ITNode* &node
      , Interval &query_interval
      , ret_value_type& ret
    )
    {
        if (node != nullptr) 
        {
            /*
                0 & default: This situation could not happed.
                1: The interval is smaller than node.median.
                Keep traversing to node->left.
                2: The interval is larger than node.median.
                Keep traversing to node->right.
                3: The interval is Overlapped with node.median.
                Keep traversing to node->left & node->right.
            */
            variant_type var = query_interval.get_data_ptr();
            UInt condition = 
                    (Interval::get_start(var) <= node->median) + 
                    ((Interval::get_end(var) >= node->median) << 1);
            
            switch(condition) 
            {
                case 0:
                    throw Overlap_Annotator(
                        "Interval_Tree Error:"
                        " Errors occure when finding overlaps.\n"
                    );
                    break;
                case 1:
                    store_overlap(node->vec_interval, query_interval, ret);
                    find_overlap(node->left, query_interval, ret);
                    break;
                case 2:
                    store_overlap(node->vec_interval, query_interval, ret);
                    find_overlap(node->right, query_interval, ret);
                    break;
                case 3:
                    store_overlap(node->vec_interval, query_interval, ret);
                    find_overlap(node->left, query_interval, ret);
                    find_overlap(node->right, query_interval, ret);
                    break;
                default:
                    throw Overlap_Annotator(
                        "Interval_Tree Error:"
                        " Errors occure when finding overlaps.\n"
                    );
                    break;
            }
        } 
        else 
        ;  // node is empty, do nothing
        
    }

    /**
     * @brief Find overlaps between databases and queries.
     *
     * @param intervals Intervals that are going to be searched.
     *        query_interval The interval that the user wishes to 
     *                       find overlaps.
     *        ret Pointers which point to overlapped annotations.
     */
    void store_overlap(
        std::vector<Interval>& intervals
      , Interval &query_interval
      , ret_value_type& ret
    ) 
    {
        variant_type q_var, var;
        for (auto& interval : intervals) 
        {   
            q_var = query_interval.get_data_ptr();
            var = interval.get_data_ptr();
            if (Interval::get_start(q_var) <= Interval::get_end(var) && 
                Interval::get_end(q_var) >= Interval::get_start(var)) 
            {
                std::visit([&var, &ret](auto&& arg) 
                {
                    using Q = std::decay_t<decltype(arg)>;
                    std::get<std::vector<Q>>(ret)
                                .emplace_back(std::get<Q>(var));
                }, var);
            } 
            else
            {
                // no overlap
                continue;
            }
        }
    }

};
}
