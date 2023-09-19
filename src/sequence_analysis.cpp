/*
C++ module with pybind11 binding to do sequence analysis.

NOTE: you will need this: https://stackoverflow.com/a/36056804/2112406
TODO: make things pickleable and add copy & deepcopy support (https://pybind11.readthedocs.io/en/stable/advanced/classes.html#deepcopy-support)
TODO: sequence logo class by itself?
*/

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

void init_sequence(py::module_ &);
void init_orf(py::module_ &);
void init_seq_set(py::module_ &);

class SamFile{
    // path to the .sam file
    std::string file_name;

    // filter options
    std::vector<int> start = {0}; //for each name in mapped onto, give start and end
    std::vector<int> end = {1};
    std::vector<std::string> mapped_onto = {""};
    float min_score = 0;

    // if AS is not normalized, we need to use sequence lengths
    bool normalized_score = true;
    std::unordered_map<std::string, unsigned int> lengths;

    // headers
    std::vector<std::string> headers;

    // vector that will hold the read lines
    std::vector<std::string> entries;

    // beginning and end for the reads 
    int seq_start = std::numeric_limits<int>::max();
    int seq_end = 0;
public:
    SamFile(){}
    void set_filter_options(std::vector<int> start_, std::vector<int> end_, std::vector<std::string> mapped_onto_, float min_score_)
    {
        // check that start, end, and mapped_onto have the same size
        if (!(start_.size() == end_.size() && end_.size() == mapped_onto_.size())){
            std::cout << "WARNING: start, end, and mapped_onto lists must have the same size." << std::endl;
            return;
        }
        start = start_;
        end = end_;
        mapped_onto = mapped_onto_;
        min_score = min_score_;
    }
    int size(){
        return entries.size();
    }
    void set_file_name(std::string file_name_) { file_name = file_name_;}
    std::string get_file_name(){return file_name;}
    void set_normalized_true(){normalized_score = true;}
    void set_normalized_false(){normalized_score = false;}
    bool get_normalized(){return normalized_score;}
    int get_seq_start() { return seq_start;}
    int get_seq_end() { return seq_end;}
    void get_lengths_from_fasta(std::string fasta_file_name){
        SeqSet sset;
        lengths = sset.get_names_and_lengths_from_fasta(fasta_file_name);
    }
    void read(){
        std::ifstream in_file(file_name);
        std::string line;

        // if scores are not normalized, we need to have length info
        if ( !normalized_score ){
            if (lengths.size() == 0){
                std::cout << "WARNING: scores are not normalized but no length info is given." << std::endl;
                return;
            }
        }

        // make sure file is open
        if(!in_file.is_open()){
            std::cout << "WARNING: Failed to open file with name " << file_name << std::endl;
            return;
        }

        std::string dummy, ref_name, seq, flag, seq_name;
        int pos, length;
        float score;
        std::string as_flag = "AS:i:";

        bool skip = false;

        while(std::getline(in_file, line)){
            skip = false;

            if (line[0] == '@'){
                // do stuff with headers
                continue;
            }
            std::istringstream ss(line);
            ss >> seq_name >> dummy >> ref_name >> pos >> dummy >> dummy >> dummy >> dummy >> dummy >> seq >> dummy;

            // apply first set of filters
            for(unsigned int i = 0; i < mapped_onto.size(); i++){
                if ( ref_name.compare(mapped_onto[i]) != 0 && mapped_onto[i].compare("") != 0){
                    skip = true;
                }
                length = seq.size();
                if ( pos + length < start[i] || (pos > end[i] && end[i] != -1) ){
                    skip = true;
                }
            }
            if (skip) continue;

            // check alignment score
            while ( ss >> flag){
                // if flag is AS
                if (flag.compare(0, 5, as_flag) == 0){
                    score = std::stof(flag.substr(5, flag.length()));
                    if (!normalized_score){
                        auto it = lengths.find(seq_name);
                        score = score / (float) it->second;
                    }
                    if (score < min_score){
                        skip = true;
                        break;
                    }
                }
            }
            if (skip) continue;

            // passed all the tests
            entries.push_back(line);
            seq_start = fmin(seq_start, pos);
            seq_end = fmax(seq_end, pos + length);
        }
    }
};



PYBIND11_MODULE(sequence_analysis_cpp, m){
    init_sequence(m);
    init_orf(m);
    init_seq_set(m);
    
    py::class_<SamFile>(m, "SamFile", py::dynamic_attr())
        .def(py::init<>())
        .def("set_filter_options", &SamFile::set_filter_options)
        .def("size", &SamFile::size)
        .def("__len__",
             [](SamFile &a){
                 return a.size();
             })
        .def("set_file_name", &SamFile::set_file_name)
        .def("get_file_name", &SamFile::get_file_name)
        .def_property("file_name", &SamFile::get_file_name, &SamFile::set_file_name)
        .def("set_normalized_true", &SamFile::set_normalized_true)
        .def("set_normalized_false", &SamFile::set_normalized_false)
        .def("get_lengths_from_fasta", &SamFile::get_lengths_from_fasta)
        .def("read", &SamFile::read)
        .def("get_normalized", &SamFile::get_normalized)
        .def("get_seq_start", &SamFile::get_seq_start)
        .def("get_seq_end", &SamFile::get_seq_end)
        .def("__repr__",
             [](SamFile &a){
                 if ( a.size() > 0){
                    return "<nonempty sequence_analysis.SamFile object>";
                 }else{
                     return "<empty sequence_analysis.SamFile object>";
                 }
             })
        ;
}