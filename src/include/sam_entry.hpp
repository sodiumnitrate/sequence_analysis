#pragma once
#include <string>
#include <vector>

class SamEntry{
    std::string read_name;
    std::string sequence;
    std::string mapped_onto;
    int start_pos;
    int end_pos;

    float alignment_score;
public:
    SamEntry(std::string& rn, std::string& s, std::string& m, int st, int ed, float as);
    std::string get_read_name();
    std::string get_seq_str();
    std::string get_mapped_onto();
    int get_start_pos();
    int get_end_pos();
    int get_length();
};