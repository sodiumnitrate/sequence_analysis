#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>
#include <unordered_set>
#include <iostream>
#include <cctype>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <sstream>
#include "includes.hpp"
#include "pair_score.hpp"

namespace py = pybind11;

enum algorithms {local, global};

PairwiseAligner::PairwiseAligner(){
    PairScore ps("blosum50");
};
void PairwiseAligner::set_algorithm(std::string alg) {algorithm = alg;}
std::string PairwiseAligner::get_algorithm() {return algorithm;}
void PairwiseAligner::set_query(std::string query_) {query = query_;}
std::string PairwiseAligner::get_query(){return query;}
void PairwiseAligner::set_target(std::string target_) {target = target_;}
std::string PairwiseAligner::get_target(){return target;}
float PairwiseAligner::get_score(){return score;}

void PairwiseAligner::align(){
    if (query.size() == 0 || target.size() == 0) throw "query or target not set";
    const int n = query.size() + 1;
    const int m = target.size() + 1;

    // allocate memory for F and pointers
    F.resize(n, std::vector<int>(m));
    pointers.resize(n, std::vector<direction>(m));

    if (algorithm.compare("global")==0) needleman_wunsch();
}

// Implementation of Needleman-Wunsch
void PairwiseAligner::needleman_wunsch(){
    const int n = query.size() + 1;
    const int m = target.size() + 1;

    direction pointers[n][m];

    // init matrices
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            F[i][j] = 0;
            pointers[n][m] = none;
        }
    }

    // boundaries
    for (int i = 0; i < n; i++)
    {
        F[i][0] = i * gap_penalty;
        pointers[i][0] = up;
    }
    for (int i = 0; i < m; i++)
    {
        F[0][i] = i * gap_penalty;
        pointers[0][i] = left;
    }

    pointers[0][0] = none;

    int down, right, lower_right;
    for (int i = 1; i < n; i++)
    {
        for (int j = 1; j < m; j++)
        {
            // x_i and y_j are aligned
            lower_right = F[i - 1][j - 1] + pair_score(query[i - 1], target[i - 1]);
            // x_i is aligned to a gap
            down = F[i][j - 1] + gap_penalty;
            // y_j is aligned to a gap
            right = F[i - 1][j] + gap_penalty;

            if (lower_right > down){
                if (lower_right > right){
                    F[i][j] = lower_right;
                    pointers[i][j] = upper_left;
                }
                else{
                    F[i][j] = right;
                    pointers[i][j] = left;
                }
            }
            else{
                if (down > right){
                    F[i][j] = down;
                    pointers[i][j] = up;
                }
                else{
                    F[i][j] = right;
                    pointers[i][j] = left;
                }
            }
        }
    }

    score = F[n-1][m-1];
}

void PairwiseAligner::traceback_nw(){
    int row_ptr = n - 1;
    int col_ptr = m - 1;

    query_aligned = "";
    target_aligned = "";

    while (row_ptr > 0 || col_ptr > 0){
        switch (pointers[row_ptr][col_ptr]){
            case up:
                query_aligned = "-" + query_aligned;
                target_aligned = target[row_ptr - 1] + target_aligned;
                alignment_string = "-" + alignment_string;
                row_ptr--;
                break;
            case left:
                query_aligned = query[col_ptr - 1] + query_aligned;
                target_aligned = "-" + target_aligned;
                alignment_string = "-" + alignment_string;
                col_ptr--;
                break;
            case upper_left:
                query_aligned = query[col_ptr - 1] + query_aligned;
                target_aligned = target[row_ptr - 1] + target_aligned;
                if (query[col_ptr - 1] == target[row_ptr -1]) alignment_string = "|" + alignment_string;
                else alignment_string = " " + alignment_string;
                row_ptr--;
                col_ptr--;
                break;
        }
    }
}


void init_pairwise_aligner(py::module_ &m){
    py::class_<PairwiseAligner>(m, "PairwiseAligner", py::dynamic_attr())
        .def(py::init<>())
        .def("set_algorithm", &PairwiseAligner::set_algorithm)
        .def("get_algorithm", &PairwiseAligner::get_algorithm)
        .def_property("algorithm", &PairwiseAligner::get_algorithm, &PairwiseAligner::set_algorithm)
        .def("set_query", &PairwiseAligner::set_query)
        .def("get_query", &PairwiseAligner::get_query)
        .def_property("query", &PairwiseAligner::get_query, &PairwiseAligner::set_query)
        .def("set_target", &PairwiseAligner::set_target)
        .def("get_target", &PairwiseAligner::get_target)
        .def_property("target", &PairwiseAligner::get_target, &PairwiseAligner::set_target)
        .def("get_score", &PairwiseAligner::get_score)
    ;
}
