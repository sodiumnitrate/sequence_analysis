#pragma once
#include <string>
#include <vector>
#include "orf.hpp"

// Sequence class
class Sequence {
    std::string name;
    std::string seq_str;
    std::string type;
public:
    Sequence();
    Sequence(std::string seq_str_);
    void set_name(std::string name_);
    std::string get_name();
    void set_seq(std::string seq_str_);
    std::string get_seq();
    std::string get_type();

    // return length
    int length();

    // checking equality
    bool const operator==(const Sequence& other) const;

    // checking <
    bool const operator<(const Sequence& other) const;

    // function to set type
    void set_type(std::string seq_type);

    // find codon
    int find_codon(std::string codon);

    // reverse
    Sequence reverse();

    // complement
    Sequence complement();

    // reverse complement
    Sequence reverse_complement();

    // frame shift
    Sequence frame_shift(int frame);

    // translate
    Sequence translate();

    // write fasta
    void write_fasta(std::string file_name_);

    // get orfs
    std::vector<OpenReadingFrame> get_open_reading_frames(unsigned int min_len);
};