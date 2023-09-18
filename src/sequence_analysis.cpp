/*
C++ module with pybind11 binding to do sequence analysis.
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
    {"GGG", 'G'}
};

char codon_to_aa(std::string codon){
    // requires c++20
    std::ranges::replace(codon, 'T', 'U');
    return codon_to_aa_map[codon];
}

class Sequence {
    std::string name;
    std::string seq_str;
    std::string type;
public:
    Sequence(std::string &seq_str) : seq_str(seq_str) {}
    void set_name(std::string &name_) { name = name_; }
    std::string get_name() { return name; }
    void set_seq(std::string &seq_str_) { seq_str = seq_str_; }
    std::string get_seq() { return seq_str; }
    //void set_type(std::string &type_) { type = type_; }
    std::string get_type() { return type; }

    // return length
    int length() { return seq_str.length();}

    // function to set type
    void set_type(std::string seq_type=""){
        if (seq_str.length() == 0)
        {
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
            std::unordered_set<std::string> uniq_chars;
            int nuc_ct = 0;
            for (int i = 0; i < seq_str.length(); i++)
            {
                ch = seq_str[i];
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
    int find_codon(std::string codon){
        // return the index of the first match
        if (type.compare("rna") != 0 && type.compare("dna") != 0){
            std::cout << "ERROR: cannot search for codons if the sequence isn't rna or dna." << std::endl;
            return -1;
        }
        if (codon.length() != 3){
            std::cout << "ERROR: codons must have 3 letters." << std::endl;
            return -1;
        }
        int ptr = 0;
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
    Sequence reverse(){
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
    Sequence complement(){
        // check type
        if (type.compare("rna") != 0 && type.compare("dna") != 0){
            std::cout << "ERROR: can't complement if type not rna or dna." << std::endl;
            std::string empty;
            return Sequence(empty);
        }

        std::string complemented;
        for (int i = seq_str.length() - 1; i >= 0; i--){
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
                complemented.push_back('T');
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
    Sequence reverse_complement(){
        Sequence new_seq = Sequence(seq_str).reverse().complement();
        new_seq.set_name(name);
        new_seq.set_type(type);
        return new_seq;
    }

    // frame shift
    Sequence frame_shift(int frame){
        // frame shift and return new obj
        std::cout << frame << std::endl;
        std::string new_str;
        std::cout << seq_str << std::endl;
        if (frame <= 0 || frame > 2){
            new_str = seq_str;
        }
        else if (frame == 1){
            std::cout << "now im here" << std::endl;
            std::cout << seq_str[0] << std::endl;
            std::cout << seq_str[1] << std::endl;
            new_str = seq_str.substr(1, seq_str.length() - 2);
        }
        else{
            new_str = seq_str.substr(2, seq_str.length() - 3);
        }
        std::cout << "now im there" << std::endl;
        std::cout << new_str << std::endl;
        Sequence to_return(new_str);
        to_return.set_type(type);
        to_return.set_name(name);
        return to_return;
    }

    // translate
    Sequence translate(){
        if (type.compare("rna") != 0 && type.compare("dna") != 0){
            std::cout << "ERROR: can't complement if type not rna or dna." << std::endl;
            std::string empty;
            return Sequence(empty);
        }
        std::string protein_sequence;
        int ptr = 0;
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
    void write_fasta(std::string file_name){
        std::ofstream out_file;
        out_file.open(file_name);
        out_file << ">" << name << std::endl;
        for (int i=0; i < seq_str.length(); i++){
            if( i != 0 && i % 79 == 0){
                out_file << std::endl;
            }
            out_file << seq_str[i];
        }
        out_file << std::endl;

        out_file.close();
    }

    // get orfs
    std::vector<OpenReadingFrame> get_open_reading_frames(){
        std::vector<OpenReadingFrame> result;
        if (type.compare("dna") != 0 && type.compare("rna") != 0){
            std::cout << "ERROR: open reading frames can be found for DNA or RNA sequences only." << std::endl;
            return result;
        }

        // remove gaps
        std::string seq = seq_str;
        std::ranges::replace(seq, '-', '');
        Sequence modded(seq);
        Sequence shifted(seq);
        for (int frame = 0; frame < 3; frame++){
            for (int [-1, 1] : order ){
                if (order == 1){
                    shifted = modded.frame_shift(frame);
                }
                else{
                    shifted = modded.reverse_complement().frame_shift(frame);
                }

                int len_str = shifted.length() / 3;
                int start_ind = 0;
                int stop_ind = 0;
                std::string curr_fragment;
                bool recording = false;
                for (int ptr = 0; ptr < 3*len_str; ptr+=3)
                {
                    curr_codon = shifted.get_seq().substr(ptr, 3);
                    if (codon_to_aa[curr_codon] == 'M' && !recording){
                        start_ind = ptr;
                        recording = true;
                    }
                    if (recording){
                        curr_fragment += curr_codon;
                    }
                    if (codon_to_aa[curr_codon] == '*' || ptr == 3*(len_str-1)){
                        if (recording){
                            stop_ind = ptr + 3;
                            start_ind += frame;
                            stop_ind += frame;
                            if (order == -1){
                                start_ind = seq_str.length() - start_ind;
                                stop_ind = seq_str.length() - end_ind;
                            }
                            if (start_ind > stop_ind){
                                int dummy = stop_ind;
                                stop_ind = start_ind;
                                start_ind = dummy;
                            }
                            OpenReadingFrame orf;
                            Sequence prot = shifted.translate();
                            orf.set_props(shifted, modded, prot, start, stop, order, frame);
                            result.push_back(orf);
                        }

                        recording = false;
                        curr_fragment = "";
                    }
                }
            }
        }
    return result;
    }

    /*
    TODO:
    - get composition (freqs of letters)
    - interface affinity
    - weight
    - int affinity period
    - find orfs
    - read single seq from file
    */
};

class SeqSet {
    std::string name;
    // pointers to sequences!
    std::vector<Sequence> records;
    std::string type;
public:
    SeqSet() {};
    void set_name(std::string &name_) {name = name_;}
    std::string get_name() {return name;}
    int size() { return records.size(); }
    void set_type(std::string seq_type=""){
        if ( records.size() == 0){
            std::cout << "WARNING: there are no sequences in the set to determine type." << std::endl;
            return;
        }
        for (auto& t : records){
            t.set_type(seq_type);
        }
        // set the type based on the first sequence
        // (Assumes all sequences have the same type -- TODO: add check.)
        type = records[0].get_type();
    }
    std::string get_type(){ return type; };
    void set_records(std::vector<Sequence> &records_) {records = records_;}
    std::vector<Sequence> get_records() { return records;}

    // function to be able to add sequence to the set
    void add_sequence(Sequence seq){
        records.push_back(seq);
        return;
    }

    // add elements from another SeqSet instance
    void add_set(SeqSet* sset){
        for (auto t : sset->get_records()){
            records.push_back(t);
        }
    }

    // write fasta
    void write_fasta(std::string file_name){
        std::ofstream out_file;
        out_file.open(file_name);
        for (auto t : records){
            out_file << ">" << name << std::endl;
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

    void read_fasta(std::string file_name){
        std::ifstream in_file (file_name);
        std::string line;
        std::string name;
        std::string seq;

        // make sure file is opened
        if (!in_file.is_open()){
            std::cout << "WARNING: Failed to open file with name " << file_name << std::endl;
            return;
        }

        int seq_st = 0;
        while(std::getline(file, line)){
            if (line[0] == '>')
            {
                if (seq_ct > 0){
                    Sequence new_seq(seq);
                    new_seq.name = name;
                    records.push_back(new_seq);
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
        new_seq.name = name;
        records.push_back(new_seq);
    }

    /*
    TODO:
    - filter support for read_fasta
    - __getitem__
    - __add__
    - write_phylip
    - alphabetize
    - remove duplicates
    - get frequencies
    - read fastq
    - read phylip?
    - write fastq?
    - dealign (remove all gaps)
    - filter by frequency ?
    - find kmers
    */
};

// TODO: sequence logo class by itself?

class OpenReadingFrame{
    Sequence rna_sequence;
    Sequence parent_sequence;
    Sequence protein_sequence;
    int start;
    int stop;
    int strand;
    int frame;
public:
    OpenReadingFrame() {}
    void set_props(Sequence rna_sequence_, Sequence parent_sequence_, Sequence protein_sequence_, int start_, int stop_, int strand_, int frame_){
        rna_sequence = rna_sequence_;
        parent_sequence = parent_sequence_;
        protein_sequence = protein_sequence_;
        start = start_;
        stop = stop_;
        strand = strand_;
        frame = frame_;
    }
    void set_start(int start_) {start = start_;}
    int get_start() { return start;}
    void set_stop(int stop_) { stop = stop_;}
    int get_stop() {return stop;}
    void set_strand(int strand_) { strand = strand_;}
    int get_strand() {return strand;}
    // TODO: finish rest
};

PYBIND11_MODULE(sequence_analysis_cpp, m){
    py::class_<Sequence>(m, "Sequence", py::dynamic_attr())
        .def(py::init<std::string &>())
        .def("set_name", &Sequence::set_name)
        .def("get_name", &Sequence::get_name)
        .def_property("name", &Sequence::get_name, &Sequence::set_name)
        .def("set_seq", &Sequence::set_seq)
        .def("get_seq", &Sequence::get_seq)
        .def_property("seq_str", &Sequence::get_seq, &Sequence::set_seq)
        .def("get_type", &Sequence::get_seq)
        .def("set_type", &Sequence::set_type)
        .def_property("type", &Sequence::get_type, &Sequence::set_type)
        .def("find_codon", &Sequence::find_codon)
        .def("reverse", &Sequence::reverse)
        .def("complement", &Sequence::complement)
        .def("reverse_complement", &Sequence::reverse_complement)
        .def("frame_shift", &Sequence::frame_shift)
        .def("translate", &Sequence::translate)
        .def("write_fasta", &Sequence::write_fasta)
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

    py::class_<SeqSet>(m, "SeqSet", py::dynamic_attr())
        .def(py::init<>())
        .def("set_name", &SeqSet::set_name)
        .def("get_name", &SeqSet::get_name)
        .def_property("name", &SeqSet::get_name, &SeqSet::set_name)
        .def("set_type", &SeqSet::set_type)
        .def("get_type", &SeqSet::get_type)
        .def_property("type", &SeqSet::get_type, &SeqSet::set_type)
        .def("add_sequence", &SeqSet::add_sequence)
        .def("size", &SeqSet::size)
        .def("get_records", &SeqSet::get_records)
        .def("set_records", &SeqSet::set_records)
        .def_property("records", &SeqSet::get_records, &SeqSet::set_records)
        .def("add_set", &SeqSet::add_set)
        .def("write_fasta", &SeqSet::write_fasta)
        .def("__repr__",
             [](SeqSet &a){
                 return "<sequence_analysis.SeqSet of size " + a.size() + '>';
             })
        .def("__len__",
             [](SeqSet &a){
                 return a.size();
             })
        ;
        //TODO: make it pickleable and add copy & deepcopy support (https://pybind11.readthedocs.io/en/stable/advanced/classes.html#deepcopy-support)
}