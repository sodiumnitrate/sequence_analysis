#pragma once
#include <string>
#include <vector>
#include "sequence.hpp"

class SamEntry{
    std::string read_name;
    std::string sequence;
    std::string mapped_onto;
    int start_pos;
    int end_pos;

    float alignment_score;
    int bin_flag;
public:
    SamEntry(std::string& rn, std::string& s, std::string& m, int st, int ed, float as, int bin_flag);
    std::string get_read_name();
    std::string get_seq_str();
    std::string get_mapped_onto();
    int get_start_pos();
    int get_end_pos();
    int get_length();
    int get_binary_flag();
    Sequence to_sequence();
};
