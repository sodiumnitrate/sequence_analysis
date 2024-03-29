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
#include "include/genome_map.hpp"

namespace py = pybind11;

GenomeMap::GenomeMap(){};
void GenomeMap::set_chromosome_name(std::string ch_name){chromosome_name = ch_name;}
std::string GenomeMap::get_chromosome_name(){return chromosome_name;}
std::vector<unsigned int> GenomeMap::get_heatmap(int start, int end){
    if(end < start){
        int dummy = end;
        end = start;
        start = end;
    }

    if (start > heatmap_end || end < heatmap_start){
        std::cout << "ranges don't overlap" << std::endl;
        throw;
    }

    // start and end are *not* relative to vector indices. They are absolute nucleotide indices
    if (start < heatmap_start || end > heatmap_end){
        std::cout << "WARNING: requested range=(" << start <<','<<end<<") is outside the heatmap range=(" << heatmap_start << ',' << heatmap_end << ")." << std::endl;
        start = heatmap_start;
        end = heatmap_end;
    }

    auto first = heatmap.begin() + (start - heatmap_start);
    auto last = heatmap.begin() + (end + 1 - heatmap_start);

    std::vector<unsigned int> sub(first, last);
    return sub;
}

void GenomeMap::set_mapped_entries(std::vector<SamEntry*> mapped_entries_){
    mapped_entries = mapped_entries_;
}

std::vector<std::string> GenomeMap::get_mapped_read_names(int start, int end){
    // start and end are *not* relative to vector indices. They are absolute nucleotide indices
    std::vector<std::string> result;

    if (start < heatmap_start || end > heatmap_end){
        std::cout << "WARNING: requested range=(" << start <<','<<end<<") is outside the heatmap range=(" << heatmap_start << ',' << heatmap_end << ")." << std::endl;
        return result;
    }

    // pick the mapped entries that overlap with the requested range
    int start_idx, end_idx;
    std::string read_name;
    for (auto& t : mapped_entries){
        start_idx = t->get_start_pos();
        end_idx = t->get_end_pos();
        read_name = t->get_read_name();
        if (start_idx < end && end_idx > start){
            result.push_back(read_name);
        } 
    }

    return result;
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
    // TODO: complete!
    int new_gm_start = new_gm->heatmap_start;
    int new_gm_end = new_gm->heatmap_end;
}

void GenomeMap::set_from_list(std::vector<int>& starts, std::vector<int>& ends){
    if (starts.size() != ends.size()){
        std::cout << "number of start indices not the same as the number of end indices" << std::endl;
        throw;
    }
    heatmap_start = * std::min_element(starts.begin(), starts.end());
    heatmap_end = *  std::max_element(ends.begin(), ends.end());

    const unsigned int length = heatmap_end - heatmap_start + 1;

    heatmap.resize(length);
    std::fill(heatmap.begin(), heatmap.end(), 0);

    for (auto i = 0; i < starts.size(); i++){
        for (int j = starts[i]; j <= ends[i]; j++){
            heatmap[j - heatmap_start] += 1;
        }
    }
}

void GenomeMap::set_from_list(std::vector<int>& starts, std::vector<int>& ends, std::vector<int>& multiplicities){
    if (starts.size() == 0){
        std::cout << "empty input" << std::endl;
        throw;
    }
    if (starts.size() != ends.size()){
        std::cout << "number of start indices not the same as the number of end indices" << std::endl;
        throw;
    }
    if (starts.size() != multiplicities.size()){
        std::cout << "number of indices does not match the number of multiplicities";
        throw;
    }

    heatmap_start = * std::min_element(starts.begin(), starts.end());
    heatmap_end = *  std::max_element(ends.begin(), ends.end());

    const unsigned int length = heatmap_end - heatmap_start + 1;

    heatmap.resize(length);
    std::fill(heatmap.begin(), heatmap.end(), 0);

    for (auto i = 0; i < starts.size(); i++){
        for (int j = starts[i]; j <= ends[i]; j++){
            heatmap[j - heatmap_start] += multiplicities[i];
        }
    }
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
        .def("get_mapped_read_names", &GenomeMap::get_mapped_read_names)
        .def("set_from_list", py::overload_cast<std::vector<int>&, std::vector<int>&>(&GenomeMap::set_from_list))
        .def("set_from_list", py::overload_cast<std::vector<int>&, std::vector<int>&, std::vector<int>&>(&GenomeMap::set_from_list))
        ;
}
