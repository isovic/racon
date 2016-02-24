/**
 * @file alignment.hpp
 * @author Marko Culinovic <marko.culinovic@gmail.com>
 * @brief Alignment class header file 
 * @details Header file with declaration of Alignment class used for
 * calculating local alignment between graph and sequence. Class is
 * based on https://github.com/ljdursi/poapy/blob/master/seqgraphalignment.py
 * python implementation
 */
#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <vector>
#include <string>
#include <tuple>
#include <deque>
#include <unordered_map>

#include "./graph.hpp"

using std::deque;
using std::pair;
using std::string;
using std::tuple;
using std::unordered_map;
using std::vector;

namespace POA {

    typedef tuple<int, int, int, char> move;
    typedef pair<int, int> limit;

    /**
     * @brief dp_el struct 
     * @details stores all information related to particular element in DP matrix.
     */
    struct dp_el {
        int score;
        int prev_graph_idx;
        int prev_seq_idx;
        char insert_cost;
        char delete_cost;

        dp_el(const int score, const int prev_graph_idx, const int prev_seq_idx, const int insert_cost, const int delete_cost)
          : score(score), prev_graph_idx(prev_graph_idx), prev_seq_idx(prev_seq_idx), insert_cost(insert_cost), delete_cost(delete_cost)
        {}
    };

    /**
     * @brief Alignment class
     * @details Calculates local alignment between graph nodes and sequence
     * using modified smith-waterman algorithm.
     * 
     */
    class Alignment {
     public:
        /**
         * @brief Alignment constructor
         * @details Constructor needs two parameters - sequence and graph
         * for calculating alignment between those two.
         * 
         * @param sequence reference to string sequence
         * @param graph reference to graph
         */
        Alignment(const string& sequence, const Graph& graph);


        /**
         * @brief Aligns sequence to graph
         * @details Method performs local alignment between
         * string sequence and directed acyclic graph (DAG).
         * Algorithm is variation of smith-waterman algorithm 
         * (http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm).
         */
        void align();


        /**
         * @brief Aligns sequence to part of graph which represents part of sequence
         * on position >= pos.
         * @details Method is created to improve efficieny while aligning sequence
         * that should occur in some later parts of the graph.
         * Example:
         * AAAAAAAAAAABBBBBBFBBCCCCDC (seq1)
         *                  BBBCCCDCC (seq2)
         * seq1.indexOf('F') - THRESHOLD would be a nice 'pos' to pass since more than
         * a half of string won't be part od DP.
         *
         * Method performs local alignment between
         * string sequence and directed acyclic graph (DAG).
         * Algorithm is variation of smith-waterman algorithm
         * (http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm).
         * @param pos minimum position where to align
         */
        void align_starting_at(const uint32_t pos);


        /**
         * @brief Aligns sequence to part of graph which represents part of sequence
         * on position >= pos, with given band width.
         * @details Method is created to improve efficieny while aligning sequence
         * that should occur in some later parts of the graph.
         * Example:
         * AAAAAAAAAAABBBBBBFBBCCCCDC (seq1)
         *                  BBBCCCDCC (seq2)
         * seq1.indexOf('F') - THRESHOLD would be a nice 'pos' to pass since more than
         * a half of string won't be part od DP.
         *
         * Banded alignment explained in comment on https://www.biostars.org/p/2071/.
         * Essentialy, just the "diagonal" of 2*band_width will be considered for DP.
         *
         * Method performs local alignment between
         * string sequence and directed acyclic graph (DAG).
         * Algorithm is variation of smith-waterman algorithm
         * (http://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm).
         * @param pos minimum position where to align
         * @param band_width the band_width, -1 for no limit
         */
        void align_banded_starting_at(const uint32_t pos, const int band_width);


        /**
         * @brief Getter for sequence
         * @return sequence
         */
        const string& sequence() const { return sequence_; }


        /**
         * @brief Getter for sequence indices in alignment
         * @details Method returns sequence indices representing
         * alignment. Index value is -1 if node in graph isn't aligned
         * to any position/character in sequence, otherwise it is position
         * in sequence string.
         * @return sequence indices
         */
        const deque<int>& seq_idxs() const { return seq_idxs_; }


        /**
         * @brief Getter for node ids in alignment
         * @details Methor returns node ids representing
         * alignment. Id value is -1 if some character in sequence
         * isn't aligned to any node in graph, otherwise it is id of node
         * in graph.
         * @return node ids
         */
        const deque<int>& node_ids() const { return node_ids_; }


     private:
        // value of match in smith-waterman score matrix
        static int match_score_;
        // value of mismatch in smith-waterman score matrix
        static int mismatch_score_;
        // value of gap opening in smith-waterman score matrix
        static int open_gap_score_;
        // value of extending gap in smith-waterman score matrix
        static int extend_gap_score_;

        // sequence to be aligned to graph
        const string& sequence_;
        // poa graph
        const Graph& graph_;

        // smith-waterman dynamic programming "matrix" storage (condensed)
        vector<dp_el> dp_;

        // sequence indices representing alignment
        deque<int> seq_idxs_;
        // node ids representing alignment
        deque<int> node_ids_;

        // map node ID to index in dp matrix
        vector<int> nodeID_to_index_;
        // map index in dp matrix to node ID
        vector<uint32_t> index_to_nodeID_;
        // maps node_index -> sequence.lo, sequence.hi that has to be aligned to node on dp[node_index]
        vector<limit> seq_limits_;

        // max score in matrix
        int max_score_;
        // row index of maximum score
        int max_i_;
        // column index of maximum score
        int max_j_;

        // number of graph nodes that could be aligned with the sequence.
        // All nodes that occur before given start position for sure will not be accounted here.
        uint32_t valid_nodes_num_;

        // minimal valid node_id
        // used for compressing nodeID_to_index_ field
        int max_valid_node_id_;

        // maximal valid node_id
        // used for compressing nodeID_to_index_ field
        int min_valid_node_id_;

        // maximum dp row width.
        uint32_t dp_width_;


        /**
         * @brief Returns index of given node_id in DP table. If node is
         * not part of DP, returns -1.
         *
         * @param node_id node id
         * @return index in DP table if exists, otherwise -1.
         */
        int index_from_node_id(const uint32_t node_id) const;


        /**
         * @brief Initializes matrices and all needed values for smith-waterman
         * @details Initializes dynamic programming matrix, and sequence substring
         * coordinates for each node.
         * Since some nodes are out of DP (if aligning starts from pos > 0),
         * matrix is "sparse".
         * @param starting_pos minimum starting position in result sequence
         * where the query sequence is going to be aligned.
         * @band_width aligning band width
         */
        void init_dp(const int starting_pos, const int band_width);


        /**
         * @brief Backtracks best alignment path
         */
        void backtrack();


        /**
         * @brief returns DP element
         * @param i - dp row
         * @param j - dp column
         * @returns dp[i][j]
         */
        inline dp_el& dp(const uint32_t i, const uint32_t j);
    };
}

#endif  // ALIGNMENT_H
