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

namespace py = pybind11;

PairwiseAligner::PairwiseAligner(){};
void PairwiseAligner::set_algorithm(std::string alg) {algorithm = alg;}
std::string PairwiseAligner::get_algorithm() {return algorithm;}
void PairwiseAligner::set_query(std::string query_) {query = query_;}
std::string PairwiseAligner::get_query(){return query;}
void PairwiseAligner::set_target(std::string target_) {target = target_;}
std::string PairwiseAligner::get_target(){return target;}
float PairwiseAligner::get_score(){return score;}


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
