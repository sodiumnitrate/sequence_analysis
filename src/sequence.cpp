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
#include "include/sequence.hpp"

namespace py = pybind11;

std::unordered_map<std::string, char> codon_to_aa_map = {
    {"UUU", 'F'},
    {"UUC", 'F'},
    {"UUA", 'L'},
    {"UUG", 'L'},
    {"CUU", 'L'},
    {"CUC", 'L'},
    {"CUA", 'L'},
    {"CUG", 'L'},
    {"AUU", 'I'},
    {"AUC", 'I'},
    {"AUA", 'I'},
    {"AUG", 'M'},
    {"GUU", 'V'},
    {"GUC", 'V'},
    {"GUA", 'V'},
    {"GUG", 'V'},
    {"UCU", 'S'},
    {"UCC", 'S'},
    {"UCA", 'S'},
    {"UCG", 'S'},
    {"CCU", 'P'},
    {"CCC", 'P'},
    {"CCA", 'P'},
    {"CCG", 'P'},
    {"ACU", 'T'},
    {"ACC", 'T'},
    {"ACA", 'T'},
    {"ACG", 'T'},
    {"GCU", 'A'},
    {"GCC", 'A'},
    {"GCA", 'A'},
    {"GCG", 'A'},
    {"UAU", 'Y'},
    {"UAC", 'Y'},
    {"UAA", '*'},
    {"UAG", '*'},
    {"CAU", 'H'},
    {"CAC", 'H'},
    {"CAA", 'Q'},
    {"CAG", 'Q'},
    {"AAU", 'N'},
    {"AAC", 'N'},
    {"AAA", 'K'},
    {"AAG", 'K'},
    {"GAU", 'D'},
    {"GAC", 'D'},
    {"GAA", 'E'},
    {"GAG", 'E'},
    {"UGU", 'C'},
    {"UGC", 'C'},
    {"UGA", '*'},
    {"UGG", 'W'},
    {"CGU", 'R'},
    {"CGC", 'R'},
    {"CGA", 'R'},
    {"CGG", 'R'},
    {"AGU", 'S'},
    {"AGC", 'S'},
    {"AGA", 'R'},
    {"AGG", 'R'},
    {"GGU", 'G'},
    {"GGC", 'G'},
    {"GGA", 'G'},
    {"GGG", 'G'},
    {"TTT", 'F'},
    {"TTC", 'F'},
    {"TTA", 'L'},
    {"TTG", 'L'},
    {"CTT", 'L'},
    {"CTC", 'L'},
    {"CTA", 'L'},
    {"CTG", 'L'},
    {"ATT", 'I'},
    {"ATC", 'I'},
    {"ATA", 'I'},
    {"ATG", 'M'},
    {"GTT", 'V'},
    {"GTC", 'V'},
    {"GTA", 'V'},
    {"GTG", 'V'},
    {"TCT", 'S'},
    {"TCC", 'S'},
    {"TCA", 'S'},
    {"TCG", 'S'},
    {"CCT", 'P'},
    {"ACT", 'T'},
    {"GCT", 'A'},
    {"TAT", 'Y'},
    {"TAC", 'Y'},
    {"TAA", '*'},
    {"TAG", '*'},
    {"CAT", 'H'},
    {"AAT", 'N'},
    {"GAT", 'D'},
    {"TGT", 'C'},
    {"TGC", 'C'},
    {"TGA", '*'},
    {"TGG", 'W'},
    {"CGT", 'R'},
    {"AGT", 'S'},
    {"GGT", 'G'}
};

bool is_stop(std::string &codon){
    if (codon[0] != 'U' && codon[0] != 'T') return false;
    if (codon[1] != 'A' && codon[1] != 'G') return false;
    if (codon[2] != 'A' && codon[2] != 'G') return false;
    if (codon[1] == 'G' && codon[2] == 'G') return false;
    return true;
}

bool is_start(std::string &codon){
    if (codon[0] != 'A') return false;
    if (codon[1] != 'T' && codon[1] != 'U') return false;
    if (codon[2] != 'G') return false;
    return true;
}

char codon_to_aa(std::string &codon){
    return codon_to_aa_map[codon];
}


Sequence::Sequence(std::string seq_str_) {
    std::transform(seq_str_.begin(), seq_str_.end(), seq_str_.begin(), ::toupper);
    seq_str = seq_str_;
}

void Sequence::set_name(std::string name_) { name = name_; }
std::string Sequence::get_name() { return name; }
void Sequence::set_seq(std::string seq_str_) { seq_str = seq_str_; }
std::string Sequence::get_seq() { return seq_str; }
std::string Sequence::get_type() { return type; }

// return length
int Sequence::length() { return seq_str.length();}

// checking equality
bool const Sequence::operator==(const Sequence& other) const{
    if (seq_str.compare(other.seq_str) == 0) return true;
    else return false;
}

    // checking <
bool const Sequence::operator<(const Sequence& other) const{
    if (seq_str < other.seq_str) return true;
    else return false;
}

// function to set type
void Sequence::set_type(std::string seq_type=""){
    if (seq_str.length() == 0)
    {
        std::cout << "WARNING: empty sequence. Can't set type." << std::endl;
        return;
    }
    if (seq_type.compare("") != 0){
        type = seq_type;
    }
    else
    {
        // assign type
        std::string ch;
        std::unordered_set<std::string> uniq_chars;
        int nuc_ct = 0;
        for (unsigned int i = 0; i < seq_str.length(); i++)
        {
            ch = seq_str[i];
            uniq_chars.insert(ch);

            // letters unique to protein sequences
            if (ch == "Q" || ch == "E" || ch == "H" || ch == "I" ||ch == "L" || ch == "F" || ch == "P" || ch == "O"){
                type = "protein";
                return;
            }
            if (ch == "*" || ch == "-"){
                continue;
            }
            if (ch == "A" || ch == "G" || ch == "T" || ch == "U" ||ch == "C" || ch == "N"){
                nuc_ct++;
            }
        }
        if (nuc_ct == 0){
            // highly unlikely, but check anyway
            type = "protein";
            return;
        }
        // frequency of nucleotide-associated letters
        float nuc_freq = ((float) nuc_ct) / (float) seq_str.length();
        if ( uniq_chars.size() < 10 && nuc_freq > 0.9){
            if (uniq_chars.contains("U")){
                type = "rna";
                return;
            }
            else{
                type = "dna";
                return;
            }
        }
    }
}

    // find codon
int Sequence::find_codon(std::string codon){
    // return the index of the first match
    if (type.compare("rna") != 0 && type.compare("dna") != 0){
        std::cout << "ERROR: cannot search for codons if the sequence isn't rna or dna." << std::endl;
        return -1;
    }
    if (codon.length() != 3){
        std::cout << "ERROR: codons must have 3 letters." << std::endl;
        return -1;
    }
    unsigned int ptr = 0;
    std::string curr_codon;
    while ( ptr < (seq_str.length() / 3) * 3){
        curr_codon = seq_str.substr(ptr, 3);
        if (curr_codon.compare(codon) == 0){
            return ptr;
        }
        ptr += 3;
    }
    // not found
    return -1;
}

// reverse
Sequence Sequence::reverse(){
    // reverse sequence and return as new obj
    std::string reversed;
    for (int i = seq_str.length() - 1; i >= 0; i--){
        reversed.push_back(seq_str[i]);
    }
    Sequence new_seq(reversed);
    new_seq.set_type(type);
    new_seq.set_name(name);
    return new_seq;
}

// complement
Sequence Sequence::complement(){
    // check type
    if (type.compare("rna") != 0 && type.compare("dna") != 0){
        std::cout << "ERROR: can't complement if type not rna or dna." << std::endl;
        std::string empty;
        return Sequence(empty);
    }
    std::string complemented;
    for (unsigned int i = 0; i < seq_str.length(); i++){
        if (seq_str[i] == 'A'){
            complemented.push_back(type.compare("rna") == 0 ? 'U' : 'T');
        }
        else if (seq_str[i] == 'C'){
            complemented.push_back('G');
        }
        else if (seq_str[i] == 'G'){
            complemented.push_back('C');
        }
        else if (seq_str[i] == 'U' || seq_str[i] == 'T'){
            complemented.push_back('A');
        }
        else{
            // if we have N, for instance
            complemented.push_back(seq_str[i]);
        }
    }
    Sequence new_seq(complemented);
    new_seq.set_type(type);
    new_seq.set_name(name);
    return new_seq;
}

// reverse complement
Sequence Sequence::reverse_complement(){
    Sequence new_seq = Sequence(seq_str).reverse().complement();
    new_seq.set_name(name);
    new_seq.set_type(type);
    return new_seq;
}
// frame shift
Sequence Sequence::frame_shift(int frame){
    // frame shift and return new obj
    std::string new_str;
    if (frame <= 0 || frame > 2){
        new_str = seq_str;
    }
    else if (frame == 1){
        new_str = seq_str.substr(1, seq_str.length() - 1);
    }
    else{
        new_str = seq_str.substr(2, seq_str.length() - 2);
    }
    Sequence to_return(new_str);
    to_return.set_type(type);
    to_return.set_name(name);
    return to_return;
}

// translate
Sequence Sequence::translate(){
    if (type.compare("rna") != 0 && type.compare("dna") != 0){
        std::cout << "ERROR: can't complement if type not rna or dna." << std::endl;
        std::string empty;
        return Sequence(empty);
    }
    std::string protein_sequence;
    unsigned int ptr = 0;
    std::string curr_codon;
    while (ptr < (seq_str.length() / 3) * 3){
        curr_codon = seq_str.substr(ptr, 3);
        protein_sequence.push_back(codon_to_aa(curr_codon));
        ptr += 3;
    }
    Sequence new_seq = Sequence(protein_sequence);
    new_seq.set_name(name);
    new_seq.set_type("protein");
    return new_seq;
}

// write fasta
void Sequence::write_fasta(std::string file_name){
    std::ofstream out_file;
    out_file.open(file_name);
    out_file << ">" << name << std::endl;
    for (unsigned int i=0; i < seq_str.length(); i++){
        if( i != 0 && i % 79 == 0){
            out_file << std::endl;
        }
        out_file << seq_str[i];
    }
    out_file << std::endl;
    out_file.close();
}

// get orfs
std::vector<OpenReadingFrame> Sequence::get_open_reading_frames(unsigned int min_len=80){

    std::vector<OpenReadingFrame> result;
    if (type.compare("dna") != 0 && type.compare("rna") != 0){
        std::cout << "ERROR: open reading frames can be found for DNA or RNA sequences only." << std::endl;
        return result;
    }
    // remove gaps
    std::string seq = seq_str;
    unsigned int seq_length = seq_str.length();
    seq.erase(remove(seq.begin(), seq.end(), '-'), seq.end());

    // declarations
    Sequence modded(seq);
    Sequence shifted(seq);
    int order_vals[2] = {-1, 1};
    int len_str, start_ind, stop_ind;
    std::string sequence, curr_fragment, curr_codon;
    bool recording;
    OpenReadingFrame orf;
    Sequence fragment(" "), prot(" ");

    for (int frame = 0; frame < 3; frame++){
        for (int order : order_vals ){
            // looping over 6 frames (3 frames * 2 order)
            if (order == 1){
                shifted = modded.frame_shift(frame);
            }
            else{
                shifted = modded.reverse_complement().frame_shift(frame);
            }

            len_str = shifted.length() / 3;
            sequence = shifted.get_seq();
            start_ind = 0;
            stop_ind = 0;
            recording = false;

            for (int ptr = 0; ptr < 3 * len_str; ptr += 3)
            {
                curr_codon = sequence.substr(ptr, 3);
                if (is_start(curr_codon) && !recording){
                    start_ind = ptr;
                    recording = true;
                }

                if (recording){
                    curr_fragment += curr_codon;
                }
                if (is_stop(curr_codon) || ptr == 3*(len_str-1)){
                    if (recording){
                        if (min_len > curr_fragment.length()){
                            recording = false;
                            curr_fragment = "";
                            continue;
                        }
                        stop_ind = ptr + 3;
                        start_ind += frame;
                        stop_ind += frame;
                        if (order == -1){
                            start_ind = seq_length - start_ind;
                            stop_ind = seq_length - stop_ind;
                        }
                        if (start_ind > stop_ind){
                            int dummy = stop_ind;
                            stop_ind = start_ind;
                            start_ind = dummy;
                        }
                        
                        fragment.set_seq(curr_fragment);
                        fragment.set_type();
                        prot = fragment.translate();
                        //orf.set_props(shifted.get_seq(), modded.get_seq(), prot.get_seq(), start_ind, stop_ind, order, frame);
                        //orf.set_props(shifted.get_seq(), prot.get_seq(), start_ind, stop_ind, order, frame);
                        orf.set_props(fragment.get_seq(), prot.get_seq(), start_ind, stop_ind, order, frame);
                        result.push_back(orf);
                        
                        recording = false;
                        curr_fragment = "";
                    }
                }
            }
            
        }
    }

    return result;
}

void init_sequence(py::module_ &m){
    py::class_<Sequence>(m, "Sequence", py::dynamic_attr())
        .def(py::init<std::string &>())
        .def("set_name", &Sequence::set_name)
        .def("get_name", &Sequence::get_name)
        .def_property("name", &Sequence::get_name, &Sequence::set_name)
        .def("set_seq", &Sequence::set_seq)
        .def("get_seq", &Sequence::get_seq)
        .def_property("seq_str", &Sequence::get_seq, &Sequence::set_seq)
        .def("get_type", &Sequence::get_seq)
        .def("set_type", &Sequence::set_type, py::arg("seq_type") = "")
        .def_property("type", &Sequence::get_type, &Sequence::set_type)
        .def("find_codon", &Sequence::find_codon)
        .def("reverse", &Sequence::reverse)
        .def("complement", &Sequence::complement)
        .def("reverse_complement", &Sequence::reverse_complement)
        .def("frame_shift", &Sequence::frame_shift)
        .def("translate", &Sequence::translate)
        .def("write_fasta", &Sequence::write_fasta)
        .def("get_open_reading_frames", &Sequence::get_open_reading_frames, py::arg("min_len")=80)
        .def("__repr__",
            [](Sequence &a){
                return "<sequence_analysis.Sequence " + a.get_seq().substr(0,5) + "... >";
            })
        .def("__len__", 
             [](Sequence &a){
                return a.length();
             })
        .def("__str__",
             [](Sequence &a){
                 return a.get_seq();
             })
        ;
}
