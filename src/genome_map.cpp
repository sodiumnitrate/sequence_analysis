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

GenomeMap::GenomeMap(){};
void GenomeMap::set_chromosome_name(std::string ch_name){chromosome_name = ch_name;}
std::string GenomeMap::get_chromosome_name(){return chromosome_name;}
std::vector<unsigned int> GenomeMap::get_heatmap(int start, int end){
    // start and end are *not* relative to vector indices. They are absolute nucleotide indices
    if (start < heatmap_start || end > heatmap_end){
        std::cout << "WARNING: requested range=(" << start <<','<<end<<") is outside the heatmap range=(" << heatmap_start << ',' << heatmap_end << ")." << std::endl;
        std::vector<unsigned int> blank;
        return blank;
    }

    auto first = heatmap.begin() + (start - heatmap_start);
    auto last = heatmap.begin() + (end + 1 - heatmap_start);

    std::vector<unsigned int> sub(first, last);
    return sub;
}
void GenomeMap::set_sample_name(std::string samp){sample_name = samp;}
std::string GenomeMap::get_sample_name(){return sample_name;}
// the following method sets heatmap -- to use from SamFile and not expose it to python
void GenomeMap::set_heatmap(std::vector<unsigned int> heatmap_, int heatmap_start_, int heatmap_end_){
    heatmap = heatmap_;
    heatmap_start = heatmap_start_;
    heatmap_end = heatmap_end_;
}

// add info from another genome map
void GenomeMap::add_map(GenomeMap* new_gm){
    int new_gm_start = new_gm->heatmap_start;
    int new_gm_end = new_gm->heatmap_end;
}

void init_genome_map(py::module_ &m){
    py::class_<GenomeMap>(m, "GenomeMap", py::dynamic_attr())
        .def(py::init<>())
        .def("get_chromosome_name", &GenomeMap::get_chromosome_name)
        .def("set_chromosome_name", &GenomeMap::set_chromosome_name)
        .def_property("chromosome_name", &GenomeMap::get_chromosome_name, &GenomeMap::set_chromosome_name)
        .def("get_sample_name", &GenomeMap::get_sample_name)
        .def("set_sample_name", &GenomeMap::set_sample_name)
        .def_property("sample_name", &GenomeMap::get_sample_name, &GenomeMap::set_sample_name)
        .def("get_heatmap", &GenomeMap::get_heatmap)
        ;
}
