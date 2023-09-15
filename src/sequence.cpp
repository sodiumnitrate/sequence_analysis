/*
C++ class and functions for a sequence object.
*/

#include <pybind11/pybind11.h>
#include <vector>
#include <string>
#include <pybind11/stl.h>
#include <unordered_set>
#include "utils.h"
#include <iostream>
#include <cctype>

namespace py = pybind11;

class Sequence {
    std::string name;
    std::string seq_string;
    std::string type;

    // phred style quality string, if available
    std::string quality;

public:
    // constructor
    Sequence(std::string seq_string_, std::string name_){
        name = name_;
        seq_string = seq_string_;
    }
    // function to get length
    int length(){
        return seq_string.length();
    }
    // set type
    void set_type(std::string seq_type=""){
        if (seq_string.length() == 0){
            std::cout << "WARNING: empty sequence. Can't set type." << std::endl;
            return;
        }
        if (seq_type != ""){
            type = seq_type;
        }
        else
        {
            // assign type
            std::string ch;
            int nuc_ct = 0;
            std::unordered_set<std::string> uniq_chars;
            for (int i = 0; i < seq_string.length(); i++){
                ch = seq_string[i];
                uniq_chars.insert(ch);
                // letters unique to protein sequences
                if (ch == "Q" || ch == "E" || ch == "H" || ch == "I" || ch == "L" || ch == "F" || ch == "P" || ch == "O"){
                    type = "protein";
                    return;
                }
                if (ch == "*" || ch == "-"){
                    continue;
                }
                if (ch == "A" || ch == "G" || ch == "T" || ch == "U" || ch == "C" || ch == "N"){
                    nuc_ct++;
                }
            }
            if (nuc_ct == 0){
                // highly unlikely, but check anyway
                type = "protein";
                return;
            }
            float nuc_freq = ((float) nuc_ct) / (float) seq_string.length();
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

    // get type
    std::string get_type(){
        return type;
    }

    // get sequence
    std::string get_sequence(){
        return seq_string;
    }

    // set sequence
    void set_sequence(std::string seq){
        seq_string = seq;
    }

    // get name
    std::string get_naame(){
        return name;
    }

    // set name
    void set_name(std::string name_){
        name = name_;
    }

    // translate
    Sequence translate(){
        std::string protein_sequence = "";
        int ptr = 0;
        std::string curr_codon;
        while (ptr < seq_string.length() - 3){
            curr_codon = seq_string.substr(ptr, 3);
            protein_sequence += codon_to_aa(curr_codon);
        }
        Sequence result(protein_sequence, name);
        result.type = "protein";
        return result;
    }

    // find codon
    int find_codon(std::string codon){
        // return the index of the first match
        if (type.compare("rna") != 0 || type.compare("dna") != 0 ){
            // TODO: better way to handle exception?
            std::cout << "ERROR: cannot search for codonds if the sequence isn't rna or dna." << std::endl;
            return -1;
        }
        if (codon.length() != 3){
            std::cout << "ERROR: codon must have 3 letters." << std::endl;
            return -1;
        }
        int ptr = 0;
        std::string curr_codon;
        while ( ptr < seq_string.length() - 3){
            curr_codon = seq_string.substr(ptr, 3);
            if (curr_codon.compare(codon) == 0){
                return ptr;
            }
            ptr += 3;
        }
        // not found
        return -1;
    }

    // reverse complement
    Sequence reverse_complement(){
        // reverse complement and return new obj
        std::string reversed = "";
        for (int i = seq_string.length(); i >= 0; i--){
            reversed += seq_string[i];
        }

        Sequence new_seq(reversed, name);
        new_seq.set_type(type);
        return new_seq;
    }

    // frame shift
    Sequence frame_shift(int frame){
        // frame shift and return new obj
        if (frame == 0 || frame < 0 || frame > 2){
            Sequence to_return(seq_string, name);
            to_return.type = type;
            return to_return;
        }
        else if (frame == 1){
            Sequence to_return(seq_string.substr(1, seq_string.length()-1), name);
            to_return.type = type;
            return to_return;
        }
        else if (frame == 2){
            Sequence to_return(seq_string.substr(2, seq_string.length()-1), name);
            to_return.type = type;
            return to_return;
        }
    }

};

class SeqSet {
    // props
    std::vector<Sequence> records;
    std::string type;
public:
    // methods
    SeqSet (std::vector<Sequence> records_){
        records = records_;
    }

    int size(){
        return records.size();
    }
};


PYBIND11_MODULE(sequence_module, module_handle){
    module_handle.doc() = "Sequence class.";

    py::class_<Sequence>(
        module_handle, "PySequence"
    ).def(py::init<std::string, std::string>())
    .def("length", &Sequence::length)
    .def("set_name", &Sequence::set_name)
    .def("get_name", &Sequence::get_name)
    .def("set_type", &Sequence::set_type)
    .def("get_type", &Sequence::get_type)
    .def("get_sequence", &Sequence::get_sequence)
    .def("set_sequence", &Sequence::set_sequence)
    .def("translate", &Sequence::translate)
    .def("find_codon", &Sequence::find_codon)
    .def("reverse_complement", &Sequence::reverse_complement)
    .def("frame_shift", &Sequence::frame_shift)
    ;

    py::class_<SeqSet>(
        module_handle, "PySeqSet"
    ).def(py::init<std::vector<Sequence>>())
    .def("size", &SeqSet::size)
    ;

}