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

// constructor
SamFile::SamFile(){}

// set filter options for reading sam files
void SamFile::set_filter_options(std::vector<int> start_, std::vector<int> end_, std::vector<std::string> mapped_onto_, float min_score_)
{
    // check that start, end, and mapped_onto have the same size
    if (!(start_.size() == end_.size() && end_.size() == mapped_onto_.size())){
        throw "ERROR: start, end, and mapped_onto lists must have the same size.";
    }
    start_indices = start_;
    end_indices = end_;
    mapped_onto = mapped_onto_;
    min_score = min_score_;
}

// get size of sam_file
int SamFile::size(){
    return entries.size();
}

void SamFile::set_file_name(std::string file_name_) { file_name = file_name_;}

std::string SamFile::get_file_name(){return file_name;}

void SamFile::set_normalized_true(){normalized_score = true;}

void SamFile::set_normalized_false(){normalized_score = false;}

bool SamFile::get_normalized(){return normalized_score;}

int SamFile::get_seq_start() { return seq_start;}

int SamFile::get_seq_end() { return seq_end;}

// read in sequence lengths from a fasta file (for when normalization of alignment scores is needed)
void SamFile::get_lengths_from_fasta(std::string fasta_file_name){
    SeqSet sset;
    lengths = sset.get_names_and_lengths_from_fasta(fasta_file_name);
}

// read info from sam file
void SamFile::read(){
    std::cout << min_score << std::endl;
    // release GIL to be able to run this parallel from python side
    py::gil_scoped_release release;
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

    // some preprocessing to facilitate filtering
    SamFilter filter_obj(mapped_onto, start_indices, end_indices);

    while(std::getline(in_file, line)){
        //skip = false;
        if (line[0] == '@'){
            // do stuff with headers
            if(line[1] == 'S' && line[2] == 'Q'){
                std::istringstream ss(line);
                std::string target_name, dummy;
                ss >> dummy >> target_name >> dummy;
                target_name = target_name.substr(3, target_name.length());
                if (filter_obj.check_name(target_name)){
                    headers.push_back(line);
                }
            }
            continue;
        }
        std::istringstream ss(line);
        ss >> seq_name >> dummy >> ref_name >> pos >> dummy >> dummy >> dummy >> dummy >> dummy >> seq >> dummy;

        auto length = seq.length();

        // apply filters----------------
        if (!(filter_obj.query(ref_name, pos, pos + length))) continue;
        // -----------------------------

        // check alignment score
        bool skip = false;
        while ( ss >> flag){
            // if flag is AS
            if (flag.compare(0, 5, as_flag) == 0){
                score = std::stof(flag.substr(5, flag.length()));
                if (!normalized_score){
                    auto it = lengths.find(seq_name);
                    if ( it == lengths.end()){
                        std::cout << "ERROR: sequence with name " << seq_name << " was not found in the lengths list." << std::endl;
                    }
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

    std::cout << "Done reading. The final range is: (" << seq_start << ", " << seq_end << ")." << std::endl;

    // re-acquire GIL
    py::gil_scoped_acquire acquire;
}

std::vector<int> SamFile::get_starts(){return start_indices;}
std::vector<int> SamFile::get_ends(){return end_indices;}
std::vector<std::string> SamFile::get_names(){return mapped_onto;}

GenomeMap SamFile::get_genome_map(std::string mapped_name, std::string sample_name){
    GenomeMap result;
    result.set_chromosome_name(mapped_name);
    result.set_sample_name(sample_name);
    
    // check entries exist
    if(entries.size() == 0){
        return result;
    }
    // check that mapped_name exists in mapped_onto
    bool exists = false;
    for(unsigned int i = 0; i < mapped_onto.size(); i++){
        if(mapped_name.compare(mapped_onto[i]) == 0){
            exists = true;
            break;
        }
    }
    if (!exists){ std::cout << "ERROR: requested name not found in reference names that were mapped.\n"; return result;}

    // TODO: check values of seq_end and seq_start?
    // we are ready to create our heatmap
    int length = seq_end - seq_start + 1;
    std::vector<unsigned int> heatmap;
    heatmap.resize(length);
    std::fill(heatmap.begin(), heatmap.end(), 0);

    std::string dummy, seq, ref_name;
    int pos;
    for(auto t : entries){
        std::istringstream ss(t);
        ss >> dummy >> dummy >> ref_name >> pos >> dummy >> dummy >> dummy >> dummy >> dummy >> seq;
        // skip if name doesn't match
        if(ref_name.compare(mapped_name) != 0) continue;

        for(unsigned int i = pos - seq_start; i < pos + seq.length() - seq_start; i++){
            heatmap[i] += 1;
        }
    }
    result.set_heatmap(heatmap, seq_start, seq_end);
    return result;
}

// figure out how to make it available as a common util
template<typename T>
bool are_vector_contents_equal(std::vector<T> &first, std::vector<T> &second){
    if (first.size() != second.size()){
        return false;
    }

    std::sort(first.begin(), first.end());
    std::sort(second.begin(), second.end());

    return first == second;
}

bool SamFile::are_filters_equal(SamFile* other){
    if (min_score != other->min_score) return false;
    if (!are_vector_contents_equal(mapped_onto, other->mapped_onto)) return false;
    if (!are_vector_contents_equal(start_indices, other->start_indices)) return false;
    if (!are_vector_contents_equal(end_indices, other->end_indices)) return false;

    return true;
}

void SamFile::copy_filters_from_another(SamFile* other){
    min_score = other->min_score;
    mapped_onto = other->mapped_onto;
    start_indices = other->start_indices;
    end_indices = other->end_indices;
}

// function to add info from another sam file
void SamFile::add_sam_file(SamFile* other){
    // assume filtering options are the same, or compatible
    if (!this->are_filters_equal(other)){
        std::cout << "FILTERS ARE NOT EQUAL" << std::endl;
        throw "Filters of the SamFile objects to be added are not equal.";
    }

    for (auto entry : other->entries){
        entries.push_back(entry);
    }

    std::cout << "merged entries" << std::endl;

    seq_start = fmin(seq_start, other->seq_start);
    seq_end = fmax(seq_end, other->seq_end);

    std::cout << "merged ranges" << std::endl;

    return;
}

std::vector<std::string> SamFile::get_entries(){return entries;}
std::vector<std::string> SamFile::get_headers(){return headers;}

void init_sam_file(py::module_ &m){
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
        .def("get_genome_map", &SamFile::get_genome_map)
        .def("get_entries", &SamFile::get_entries)
        .def("get_headers", &SamFile::get_headers)
        .def("__repr__",
             [](SamFile &a){
                 if ( a.size() > 0){
                    return "<nonempty sequence_analysis.SamFile object>";
                 }else{
                     return "<empty sequence_analysis.SamFile object>";
                 }
             })
        .def("add_sam_file", &SamFile::add_sam_file)
        .def("copy_filters_from_another", &SamFile::copy_filters_from_another)
        .def("get_starts", &SamFile::get_starts)
        .def("get_ends", &SamFile::get_ends)
        .def("get_names", &SamFile::get_names)
        ;
}
