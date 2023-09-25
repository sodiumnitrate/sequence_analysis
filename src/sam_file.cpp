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
        std::cout << "WARNING: start, end, and mapped_onto lists must have the same size." << std::endl;
        return;
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
        //skip = false;
        if (line[0] == '@'){
            // do stuff with headers
            if(line[1] == 'S' && line[2] == 'Q'){
                std::istringstream ss(line);
                std::string target_name, dummy;
                ss >> dummy >> target_name >> dummy;
                target_name = target_name.substr(3, target_name.length());
                for(unsigned int i = 0; i < mapped_onto.size(); i++){
                    if (target_name.compare(mapped_onto[i]) == 0){
                        headers.push_back(line);
                    }
                }
            }
            continue;
        }
        std::istringstream ss(line);
        ss >> seq_name >> dummy >> ref_name >> pos >> dummy >> dummy >> dummy >> dummy >> dummy >> seq >> dummy;

        // apply filters----------------
        skip = true;
        for(unsigned int i = 0; i < mapped_onto.size(); i++){
            // if name matches or is empty, don't skip
            if ( ref_name.compare(mapped_onto[i]) == 0 || mapped_onto[i].compare("") == 0){
                // if read overlaps with the region requested, don't skip
                length = seq.size();
                if ( pos + length > start_indices[i] && (end_indices[i] == -1 || pos < end_indices[i])){
                    skip = false;
                    break;
                }
            }
        }
        // if the conditions above are not met, skip
        if (skip) continue;
        // -----------------------------

        // check alignment score
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
}

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

// function to add info from another sam file
void SamFile::add_sam_file(SamFile* other){
    // assume filtering options are the same, or compatible
    // TODO: check for this

    for (auto entry : other->entries){
        entries.push_back(entry);
    }

    for (auto header : other->headers){
        headers.push_back(header);
    }

    seq_start = fmin(seq_start, other->seq_start);
    seq_end = fmax(seq_end, other->seq_end);

    return;
}

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
        .def("__repr__",
             [](SamFile &a){
                 if ( a.size() > 0){
                    return "<nonempty sequence_analysis.SamFile object>";
                 }else{
                     return "<empty sequence_analysis.SamFile object>";
                 }
             })
        .def("add_sam_file", &SamFile::add_sam_file)
        ;
}
