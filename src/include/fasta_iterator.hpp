#pragma once
#include <fstream>
#include "sequence.hpp"

// TODO: use this in other fasta reading ops
class FastaIterator{
    std::string file_name;
    std::ifstream file_handle;

    //  things to use to keep track of reading seqs
    int seq_ct = 0;
    std::string name;
    std::string seq;
    std::string line;
public:
    FastaIterator(std::string file_name);
    Sequence get_next();
};