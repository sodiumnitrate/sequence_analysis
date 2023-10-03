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
#include "include/fasta_iterator.hpp"

namespace py = pybind11;

// constructor
FastaIterator::FastaIterator(std::string file_name_){
    file_name = file_name_;
    file_handle.open(file_name);

    // check if file is open
    if(!file_handle.is_open()){
        throw "Failed to open file.";
    }
};

// get next sequence
Sequence FastaIterator::get_next(){
    while(std::getline(file_handle, line)){
        if(line[0] == '>')
        {
            if (!(name.empty())){
                Sequence new_seq(seq);
                new_seq.set_name(name);
                seq = "";
                // get new name
                name = line.substr(1, line.length() - 1);
                seq_ct++;
                return new_seq;
            }
            // get new name
            name = line.substr(1, line.length() - 1);
            seq_ct++;
        }
        else{
            seq += line;
        }
    }
    Sequence new_seq("");
    return new_seq;
}

void init_fasta_iterator(py::module_ &m){
    py::class_<FastaIterator>(m, "FastaIterator")
        .def(py::init<std::string>())
        .def("get_next", &FastaIterator::get_next)
        ;
}