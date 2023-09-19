/*
C++ module with pybind11 binding to do sequence analysis.

NOTE: you will need this: https://stackoverflow.com/a/36056804/2112406
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

// ORF class
class OpenReadingFrame{
    std::string rna_sequence;
    std::string parent_sequence;
    std::string protein_sequence;
    int start;
    int stop;
    int strand;
    int frame;
public:
    OpenReadingFrame() {};
    void set_props(std::string rna_sequence_, std::string parent_sequence_, std::string protein_sequence_, int start_, int stop_, int strand_, int frame_){
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
    void set_frame(int frame_) {frame = frame_;}
    int get_frame() { return frame; }
    std::string get_rna_sequence() { return rna_sequence;}
    void set_rna_sequence(std::string rs) {rna_sequence = rs;}
    std::string get_parent_sequence() {return parent_sequence;}
    void set_parent_sequence(std::string ps) { parent_sequence = ps;}
    std::string get_protein_sequence() { return protein_sequence;}
    void set_protein_sequence(std::string ps) {protein_sequence = ps;}
};


class Sequence {
    std::string name;
    std::string seq_str;
    std::string type;
public:
    Sequence(std::string &seq_str_) {
        std::transform(seq_str_.begin(), seq_str_.end(), seq_str_.begin(), ::toupper);
        seq_str = seq_str_;
    }
    void set_name(std::string &name_) { name = name_; }
    std::string get_name() { return name; }
    void set_seq(std::string &seq_str_) { seq_str = seq_str_; }
    std::string get_seq() { return seq_str; }
    //void set_type(std::string &type_) { type = type_; }
    std::string get_type() { return type; }

    // return length
    int length() { return seq_str.length();}

    // checking equality
    bool const operator==(const Sequence& other) const{
        if (seq_str.compare(other.seq_str) == 0) return true;
        else return false;
    }

    // checking <
    bool const operator<(const Sequence& other) const{
        if (seq_str < other.seq_str) return true;
        else return false;
    }

    // function to set type
    void set_type(std::string seq_type=""){
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
    Sequence reverse_complement(){
        Sequence new_seq = Sequence(seq_str).reverse().complement();
        new_seq.set_name(name);
        new_seq.set_type(type);
        return new_seq;
    }

    // frame shift
    Sequence frame_shift(int frame){
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
    Sequence translate(){
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
    void write_fasta(std::string file_name){
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
    std::vector<OpenReadingFrame> get_open_reading_frames(unsigned int min_len=80){
        std::vector<OpenReadingFrame> result;
        if (type.compare("dna") != 0 && type.compare("rna") != 0){
            std::cout << "ERROR: open reading frames can be found for DNA or RNA sequences only." << std::endl;
            return result;
        }

        // remove gaps
        std::string seq = seq_str;
        seq.erase(remove(seq.begin(), seq.end(), '-'), seq.end());
        Sequence modded(seq);
        Sequence shifted(seq);
        int order_vals[2] = {-1, 1};
        for (int frame = 0; frame < 3; frame++){
            for (int order : order_vals ){
                // looping over 6 frames (3 frames * 2 order)

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
                std::string curr_codon;
                bool recording = false;
                for (int ptr = 0; ptr < 3 * len_str; ptr += 3)
                {
                    curr_codon = shifted.get_seq().substr(ptr, 3);
                    if (codon_to_aa(curr_codon) == 'M' && !recording){
                        start_ind = ptr;
                        recording = true;
                    }
                    if (recording){
                        curr_fragment += curr_codon;
                    }
                    if (codon_to_aa(curr_codon) == '*' || ptr == 3*(len_str-1)){
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
                                start_ind = seq_str.length() - start_ind;
                                stop_ind = seq_str.length() - stop_ind;
                            }
                            if (start_ind > stop_ind){
                                int dummy = stop_ind;
                                stop_ind = start_ind;
                                start_ind = dummy;
                            }
                            OpenReadingFrame orf;
                            Sequence fragment(curr_fragment);
                            fragment.set_type();
                            Sequence prot = fragment.translate();
                            orf.set_props(shifted.get_seq(), modded.get_seq(), prot.get_seq(), start_ind, stop_ind, order, frame);
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
    std::string get_type(){ return type; };
    void set_records(std::vector<Sequence> &records_) {records = records_;
    this->set_type();}
    std::vector<Sequence> get_records() { return records;}

    // function to be able to add sequence to the set
    void add_sequence(Sequence seq){
        // TODO: add check for type
        records.push_back(seq);
        return;
    }

    // dealign
    void dealign(){
        std::string seq;
        for(auto& t : records){
            seq = t.get_seq();
            seq.erase(remove(seq.begin(), seq.end(), '-'), seq.end());
            t.set_seq(seq);
        }
        return;
    }

    // alphabetize
    void alphabetize(){
        std::sort(records.begin(), records.end());
    }

    // remove duplicates
    void remove_duplicates(){
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

    // read fasta
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

        int seq_ct = 0;
        while(std::getline(in_file, line)){
            if (line[0] == '>')
            {
                if (seq_ct > 0){
                    Sequence new_seq(seq);
                    new_seq.set_name(name);
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
        new_seq.set_name(name);
        records.push_back(new_seq);
    }

    // write phylip
    void write_phylip(std::string file_name, std::string mode){
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

    /*
    TODO:
    - filter support for read_fasta
    - __getitem__
    - __add__
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
        //TODO: make it pickleable and add copy & deepcopy support (https://pybind11.readthedocs.io/en/stable/advanced/classes.html#deepcopy-support)

    py::class_<OpenReadingFrame>(m, "OpenReadingFrame", py::dynamic_attr())
        .def(py::init<>())
        .def("set_start", &OpenReadingFrame::set_start)
        .def("get_start", &OpenReadingFrame::get_start)
        .def("set_stop", &OpenReadingFrame::set_stop)
        .def("get_stop", &OpenReadingFrame::get_stop)
        .def("set_strand", &OpenReadingFrame::set_strand)
        .def("set_strand", &OpenReadingFrame::get_strand)
        .def("set_frame", &OpenReadingFrame::set_frame)
        .def("get_frame", &OpenReadingFrame::get_frame)
        .def("set_rna_sequence", &OpenReadingFrame::set_rna_sequence)
        .def("get_rna_sequence", &OpenReadingFrame::get_rna_sequence)
        .def("set_parent_sequence", &OpenReadingFrame::set_parent_sequence)
        .def("get_parent_sequence", &OpenReadingFrame::get_parent_sequence)
        .def("set_protein_sequence", &OpenReadingFrame::set_protein_sequence)
        .def("get_protein_sequence", &OpenReadingFrame::get_protein_sequence)
        .def_property("start", &OpenReadingFrame::get_start, &OpenReadingFrame::set_start)
        .def_property("stop", &OpenReadingFrame::get_stop, &OpenReadingFrame::set_stop)
        .def_property("strand", &OpenReadingFrame::get_strand, &OpenReadingFrame::set_strand)
        .def_property("frame", &OpenReadingFrame::get_frame, &OpenReadingFrame::set_frame)
        .def_property("rna_sequence", &OpenReadingFrame::get_rna_sequence, &OpenReadingFrame::set_rna_sequence)
        .def_property("parent_sequence", &OpenReadingFrame::get_parent_sequence, &OpenReadingFrame::set_parent_sequence)
        .def_property("protein_sequence", &OpenReadingFrame::get_protein_sequence, &OpenReadingFrame::set_protein_sequence)
        .def("__repr__",
             [](OpenReadingFrame &a){
                return "<sequence_analysis.OpenReadingFrame object>";
             })
        ;        
}