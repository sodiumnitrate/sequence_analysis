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


SeqSet::SeqSet() {};
void SeqSet::set_name(std::string &name_) {name = name_;}
std::string SeqSet::get_name() {return name;}
int SeqSet::size() { return n_seqs; }
void SeqSet::set_type(std::string seq_type=""){
    if ( n_seqs == 0){
        std::cout << "WARNING: there are no sequences in the set to determine type." << std::endl;
        return;
    }
    std::string prev;
    for (auto& t : records){
        t.set_type(seq_type);
        if(prev.length() > 0 && t.get_type().compare(prev) != 0){
            std::cout << "WARNING: set contains sequences of different types! The assigned type is probably wrong." << std::endl;
            prev = t.get_type();
        }
    }
    // set the type based on the first sequence
    // (Assumes all sequences have the same type -- TODO: add check.)
    type = records[0].get_type();
}

std::string SeqSet::get_type(){ return type; };

void SeqSet::set_records(std::vector<Sequence> &records_) {
    records = records_;
    n_seqs = records.size();
    this->set_type();
}
std::vector<Sequence> SeqSet::get_records() { return records;}

// function to be able to add sequence to the set
void SeqSet::add_sequence(Sequence seq){
    // TODO: add check for type
    records.push_back(seq);
    n_seqs++;
    return;
}

// dealign
void SeqSet::dealign(){
    std::string seq;
    for(auto& t : records){
        seq = t.get_seq();
        seq.erase(remove(seq.begin(), seq.end(), '-'), seq.end());
        t.set_seq(seq);
    }
    return;
}

// alphabetize
void SeqSet::alphabetize(){
    std::sort(records.begin(), records.end());
}

// remove duplicates
void SeqSet::remove_duplicates(){
    int to_remove = 0;
    if(records.size() < 2){
        return;
    }
    // sort
    this->alphabetize();
    // flag
    std::string flag = "**";
    unsigned int idx=0;
    std::string name = records[0].get_name();
    // label duplicates, append their names to the first occurrence
    for (unsigned int i=1; i < records.size(); i++){
        if (records[i] == records[i-1]){
            name += '_' + records[i].get_name();
            records[i].set_name(flag);
            to_remove++;
        }
        else{
            if(i - idx > 1){
                records[idx].set_name(name);
            }
            name = records[i].get_name();
            idx=i;
        }
        std::cout << i << " " << name << std::endl;
    }
    // remove the labeled ones
    records.erase(std::remove_if(
        records.begin(), records.end(), 
        [](auto t){return t.get_name().compare("**")==0;}),
        records.end());

    // update count
    n_seqs -= to_remove;
}

// add elements from another SeqSet instance
void SeqSet::add_set(SeqSet* sset){
    for (auto t : sset->get_records()){
        records.push_back(t);
        n_seqs++;
    }
}

// write fasta
void SeqSet::write_fasta(std::string file_name){
    std::ofstream out_file;
    out_file.open(file_name);
    for (auto t : records){
        out_file << ">" << t.get_name() << std::endl;
        for (int i = 0; i < t.length(); i++){
            if (i != 0 && i % 79 == 0){
                out_file << std::endl;
            }
            out_file << t.get_seq()[i];
        }
        out_file << std::endl;
    }
    out_file.close();
}

// read fasta
void SeqSet::read_fasta(std::string file_name){
    std::ifstream in_file (file_name);
    std::string line;
    std::string name;
    std::string seq;
    // make sure file is opened
    if (!in_file.is_open()){
        std::cout << "WARNING: Failed to open file with name " << file_name << std::endl;
        return;
    }
    int seq_ct = 0;
    while(std::getline(in_file, line)){
        if (line[0] == '>')
        {
            if (seq_ct > 0){
                Sequence new_seq(seq);
                new_seq.set_name(name);
                records.push_back(new_seq);
                n_seqs++;
                seq = "";
            }
            // get new name
            name = line.substr(1, line.length() - 1);
            seq_ct++;
        }
        else{
            // add new line to seq
            seq += line;
        }
    }
    // add the final one
    Sequence new_seq(seq);
    new_seq.set_name(name);
    records.push_back(new_seq);
    n_seqs++;
    in_file.close();
}

// write phylip
void SeqSet::write_phylip(std::string file_name, std::string mode){
    // in phylip format, all sequences must be of the same length
    int n_seqs = records.size();
    std::unordered_set<int> n_chars;
    unsigned int max_name_chars = 0;
    for (auto t : records){
        n_chars.insert(t.length());
        // get the longest name for later use
        if(t.get_name().length() > max_name_chars) max_name_chars = t.get_name().length();
    }
    if (n_chars.size() != 1){
        std::cout << "ERROR: sequences must be of the same length." << std::endl;
        return;
    }
    std::ofstream out_file;
    out_file.open(file_name);
    auto num_chars = n_chars.begin();
    // first lines differ depending on whether we want the sequential or interleaved versions
    if (mode.compare("sequential") == 0){
        out_file << n_seqs << "\t" << *num_chars << std::endl;
    }
    else if (mode.compare("interleaved") == 0){
        out_file << n_seqs << "\t" << *num_chars << "\tI" <<std::endl; 
    }
    else{
        std::cout << "ERROR: mode not understood." << std::endl;
        return;
    }
    // sequential is pretty straightforward
    if (mode.compare("sequential") == 0){
        for (auto t : records){
            out_file << t.get_name() << "    " << t.get_seq() << std::endl;
        }
    }
    // interleaved -> write in 6 blocks of 10
    else if (mode.compare("interleaved") == 0){
        int total_lines = *num_chars / 60 + 1;
        for (int i=0; i < total_lines; i++){
            for (auto t : records){
                if (i == 0){
                    out_file << t.get_name();
                    // proper number of spaces
                    for (unsigned int j = 0; j < (max_name_chars - t.get_name().length() + 2); j++) out_file << " ";
                }
                else
                {
                    // proper number of spaces
                    for (unsigned int j = 0; j < (max_name_chars + 2); j++) out_file << " ";
                }
                // write seq_str in chunks
                int start = i * 60;
                int len = fmin(60, t.length() - start);
                std::string relevant = t.get_seq().substr(start, len);
                int count = 0;
                for(auto c : relevant){
                    if (count != 0 && count % 10 == 0){
                        out_file << " ";
                    }
                    out_file << c;
                    count++;
                }
                out_file << std::endl;
            }
            out_file << std::endl;
        }
    }
    out_file.close();
}

// read only names and lengths
std::unordered_map<std::string, unsigned int> SeqSet::get_names_and_lengths_from_fasta(std::string file_name){
    std::unordered_map<std::string, unsigned int> lengths;
    std::ifstream in_file(file_name);
    std::string line;
    std::string name;
    std::string seq;
    // TODO: refactor -- this is nearly a duplicate of read_fasta
    // make sure file is open
    if(!in_file.is_open()){
        std::cout << "WARNING: Failed to open file with name " << file_name << std::endl;
        return lengths;
    }
    int seq_ct = 0;
    while(std::getline(in_file, line)){
        if (line[0] == '>')
        {
            if (seq_ct > 0){
                lengths[name] = seq.length();
                seq = "";
            }
            name = line.substr(1, line.length() - 1);
            seq_ct++;
        }
        else{
            seq += line;
        }
    }
    // add the final one
    lengths[name] = seq.length();


    in_file.close();
    return lengths;
}

void init_seq_set(py::module_ &m){
    py::class_<SeqSet>(m, "SeqSet", py::dynamic_attr())
        .def(py::init<>())
        .def("set_name", &SeqSet::set_name)
        .def("get_name", &SeqSet::get_name)
        .def_property("name", &SeqSet::get_name, &SeqSet::set_name)
        .def("set_type", &SeqSet::set_type, py::arg("seq_type") = "")
        .def("get_type", &SeqSet::get_type)
        .def_property("type", &SeqSet::get_type, &SeqSet::set_type)
        .def("add_sequence", &SeqSet::add_sequence)
        .def("size", &SeqSet::size)
        .def("get_records", &SeqSet::get_records)
        .def("set_records", &SeqSet::set_records)
        .def_property("records", &SeqSet::get_records, &SeqSet::set_records)
        .def("add_set", &SeqSet::add_set)
        .def("write_fasta", &SeqSet::write_fasta)
        .def("read_fasta", &SeqSet::read_fasta)
        .def("write_phylip", &SeqSet::write_phylip)
        .def("alphabetize", &SeqSet::alphabetize)
        .def("dealign", &SeqSet::dealign)
        .def("remove_duplicates", &SeqSet::remove_duplicates)
        .def("__repr__",
             [](SeqSet &a){
                 return "<sequence_analysis.SeqSet of size " + std::to_string(a.size()) + " >";
             })
        .def("__len__",
             [](SeqSet &a){
                 return a.size();
             })
        ;
}