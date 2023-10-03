#pragma once
#include <string>
#include <vector>

// genome map, to be created and utilized by SamFile
// TODO: functionality to get names of reads that map to each pos?
class GenomeMap{
    std::string chromosome_name;
    std::string sample_name;

    std::vector<unsigned int> heatmap;
    int heatmap_start = 0;
    int heatmap_end = -1;
public:
    GenomeMap();
    std::vector<unsigned int> get_heatmap(int start, int end);
    std::string get_chromosome_name();
    void set_chromosome_name(std::string ch_name);
    std::string get_sample_name();
    void set_sample_name(std::string samp);
    void set_heatmap(std::vector<unsigned int> heatmap_, int heatmap_start_, int heatmap_end_);
    void add_map(GenomeMap* new_gm);
};