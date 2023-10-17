#pragma once
#include <string>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <limits>
#include "genome_map.hpp"
#include "sam_filter.hpp"
#include "seq_set.hpp"
#include "sam_entry.hpp"

class SamFile{
    // filter options
    std::vector<int> start_indices = {0}; //for each name in mapped onto, give start and end
    std::vector<int> end_indices = {-1};
    std::vector<std::string> mapped_onto = {""};
    float min_score = 0;

    // if AS is not normalized, we need to use sequence lengths
    bool normalized_score = true;
    std::unordered_map<std::string, unsigned int> lengths;
    std::vector<float> normalized_scores;

    // headers
    std::vector<std::string> headers;

    // vector that will hold the read lines
    std::vector<SamEntry> entries;

    // unordered_map to hold info about multimappers
    std::unordered_map<std::string, std::vector<int> > unique_names_to_entry_idx;

    // beginning and end for the reads 
    int seq_start = std::numeric_limits<int>::max();
    int seq_end = 0;

    // multiplicities
    std::unordered_map<std::string, int> multiplicity;
public:
    SamFile();
    void set_filter_options(std::vector<int> start_, std::vector<int> end_, std::vector<std::string> mapped_onto_, float min_score_);
    int size();
    void set_normalized_true();
    void set_normalized_false();
    bool get_normalized();
    int get_seq_start();
    int get_seq_end();
    void get_lengths_from_fasta(std::string fasta_file_name);
    void read(std::string file_name);
    GenomeMap get_genome_map(std::string mapped_name, std::string sample_name);
    void add_sam_file(SamFile* other);
    std::vector<SamEntry> get_entries();
    std::vector<std::string> get_headers();
    bool are_filters_equal(SamFile* other);
    void copy_filters_from_another(SamFile* other);
    std::vector<int> get_starts();
    std::vector<int> get_ends();
    std::vector<std::string> get_names();
    std::vector<float> get_normalized_scores();
    void generate_multimapping_stats();
    std::unordered_map<std::string, std::vector<int> > get_multimapping_stats();
    
    void set_multiplicity(std::unordered_map<std::string, int> mult);
    std::unordered_map<std::string, int> get_multiplicity();

    std::unordered_set<std::string> get_seq_names();

    SeqSet to_seq_set();
};


