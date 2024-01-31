#pragma once
#include <unordered_map>
#include "sequence.hpp"

class SeqSet {
    std::string name;
    // pointers to sequences!
    std::vector<Sequence> records;
    std::string type;
    int n_seqs = 0;
public:
    SeqSet();
    void set_name(std::string &name_);
    std::string get_name();
    int size();
    void set_type(std::string seq_type);
    std::string get_type();
    void set_records(std::vector<Sequence> &records_);
    std::vector<Sequence> get_records();

    // function to be able to add sequence to the set
    void add_sequence(Sequence seq);

    // dealign
    void dealign();

    // alphabetize
    void alphabetize();

    // remove duplicates
    void remove_duplicates();

    // add elements from another SeqSet instance
    void add_set(SeqSet* sset);

    // write fasta
    void write_fasta(std::string file_name);

    // read fasta
    void read_fasta(std::string file_name);

    // write phylip
    void write_phylip(std::string file_name, std::string mode);
    // read only names and lengths
    std::unordered_map<std::string, unsigned int> get_names_and_lengths_from_fasta(std::string file_name);

    void trim_gaps();
    void trim_gaps(float threshold);

    std::vector<std::vector<int>> pairwise_distance();
    SeqSet find_subset_with_names(std::vector<std::string> names);

    // remove columns at given indices
    void remove_columns(std::vector<int> indices);

    // check if all seqs have the same length
    bool are_lengths_identical();
};
