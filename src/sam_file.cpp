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
#include "include/sam_file.hpp"

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
void SamFile::read(std::string file_name){
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
    
    int idx = 0;

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
        if (!(filter_obj.query(ref_name, pos, pos + length - 1))) continue;
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
                    score = 100 * score / (2*(float) it->second);
                }
                if (score < min_score){
                    skip = true;
                    break;
                }
            }
        }
        if (skip) continue;
        normalized_scores.push_back(score);
        unique_names_to_entry_idx[seq_name].push_back(idx);
        idx++;

        // passed all the tests
        SamEntry entry(seq_name, seq, ref_name, pos, pos + length - 1, score);
        entries.push_back(entry);
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

std::vector<float> SamFile::get_normalized_scores(){return normalized_scores;}

void SamFile::set_multiplicity(std::unordered_map<std::string, int> mult){
    multiplicity = mult;
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
    int length = seq_end - seq_start;
    std::vector<unsigned int> heatmap;
    heatmap.resize(length);
    std::fill(heatmap.begin(), heatmap.end(), 0);

    std::string dummy, seq, ref_name, seq_name;
    int pos, end;
    int multi;
    bool overlap;
    std::vector<SamEntry*> mapped;
    for(auto& t : entries){
        overlap = false;
        ref_name = t.get_mapped_onto();
        seq_name = t.get_seq_str();
        pos = t.get_start_pos();
        end = t.get_end_pos();
        // skip if name doesn't match
        if(ref_name.compare(mapped_name) != 0) continue;

        if (multiplicity.empty()){
            multi = 1;
        }
        else multi = multiplicity[seq_name];

        if (multiplicity.find(seq_name) == multiplicity.end()) {
            std::cout << "seq_name " << seq_name << " not found in multiplicity dict." << std::endl;
            throw;
        }

        if (multi == 0) throw;

        for(unsigned int i = pos - seq_start; i < end - seq_start; i++){
            std::cout << "adding " << multi << " to heatmap at pos " << i << std::endl;
            heatmap[i] += multi;
            overlap = true;
        }

        if (overlap) mapped.push_back(&t);
    }
    result.set_heatmap(heatmap, seq_start, seq_end);
    result.set_mapped_entries(mapped);
    return result;
}

// TODO: make it available as a common util
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

std::vector<SamEntry> SamFile::get_entries(){return entries;}
std::vector<std::string> SamFile::get_headers(){return headers;}

void SamFile::generate_multimapping_stats(){
    if (entries.size() == 0) {
        std::cout << "WARNING: no entries found." << std::endl;
        return;
    }

    std::string read_name;
    int idx = 0;
    for (auto& t : entries){
        read_name = t.get_read_name();
        unique_names_to_entry_idx[read_name].push_back(idx);
        idx++;
    }
}

std::unordered_map<std::string, std::vector<int> > SamFile::get_multimapping_stats(){
    //generate_multimapping_stats();
    return unique_names_to_entry_idx;
}

std::unordered_map<std::string, int> SamFile::get_multiplicity(){
    return multiplicity;
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
        .def("set_normalized_true", &SamFile::set_normalized_true)
        .def("set_normalized_false", &SamFile::set_normalized_false)
        .def("get_lengths_from_fasta", &SamFile::get_lengths_from_fasta)
        .def("read", &SamFile::read)
        .def("get_normalized", &SamFile::get_normalized)
        .def("get_seq_start", &SamFile::get_seq_start)
        .def("get_seq_end", &SamFile::get_seq_end)
        .def("get_genome_map", &SamFile::get_genome_map)
        //.def("get_entries", &SamFile::get_entries)
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
        .def("get_normalized_scores", &SamFile::get_normalized_scores)
        .def("get_multimapping_stats", &SamFile::get_multimapping_stats)
        .def("set_multiplicity", &SamFile::set_multiplicity)
        .def("get_multiplicity", &SamFile::get_multiplicity)
        ;
}
