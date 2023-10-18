#pragma once
#include <string>
#include <vector>
#include "sam_entry.hpp"

// genome map, to be created and utilized by SamFile
// TODO: functionality to get names of reads that map to each pos?
class GenomeMap{
    std::string chromosome_name;
    std::string sample_name;

    std::vector<unsigned int> heatmap;
    int heatmap_start = 0;
    int heatmap_end = -1;

    std::vector<SamEntry*> mapped_entries;
public:
    GenomeMap();
    std::vector<unsigned int> get_heatmap(int start, int end);
    std::string get_chromosome_name();
    void set_chromosome_name(std::string ch_name);
    std::string get_sample_name();
    void set_sample_name(std::string samp);
    void set_heatmap(std::vector<unsigned int> heatmap_, int heatmap_start_, int heatmap_end_);
    void add_map(GenomeMap* new_gm);
    std::vector<std::string> get_mapped_read_names(int start, int end);
    void set_mapped_entries(std::vector<SamEntry*> mapped_entries_);

    void set_from_list(std::vector<int>& starts, std::vector<int>& ends);
    void set_from_list(std::vector<int>& starts, std::vector<int>& ends, std::vector<int>& multiplicities);
};