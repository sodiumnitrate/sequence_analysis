#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <fstream>
#include "pair_score.hpp"


// ORF class
class OpenReadingFrame{
    std::string rna_sequence;
    //std::string parent_sequence;
    std::string protein_sequence;
    int start;
    int stop;
    int strand;
    int frame;
public:
    OpenReadingFrame();
    //void set_props(std::string rna_sequence_, std::string parent_sequence_, std::string protein_sequence_, int start_, int stop_, int strand_, int frame_);
    void set_props(std::string rna_sequence_, std::string protein_sequence_, int start_, int stop_, int strand_, int frame_);
    void set_start(int start_);
    int get_start();
    void set_stop(int stop_);
    int get_stop();
    void set_strand(int strand_);
    int get_strand();
    void set_frame(int frame_);
    int get_frame();
    std::string get_rna_sequence();
    void set_rna_sequence(std::string rs);
    //std::string get_parent_sequence();
    //void set_parent_sequence(std::string ps);
    std::string get_protein_sequence();
    void set_protein_sequence(std::string ps);
};

// Sequence class
class Sequence {
    std::string name;
    std::string seq_str;
    std::string type;
public:
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

};

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

class SamFilter{
    std::unordered_set<std::string> nameset;
    std::unordered_map<std::string, std::vector<std::tuple<int, int> > > nuc_ranges;
    bool contains_empty = false;
public:
    SamFilter(std::vector<std::string>& names, std::vector<int>& starts, std::vector<int>& ends);
    bool query(std::string& name, int start, int end);
    bool check_name(std::string& name);
};

class SamFile{
    // path to the .sam file
    std::string file_name;

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
    std::vector<std::string> entries;

    // beginning and end for the reads 
    int seq_start = std::numeric_limits<int>::max();
    int seq_end = 0;
public:
    SamFile();
    void set_filter_options(std::vector<int> start_, std::vector<int> end_, std::vector<std::string> mapped_onto_, float min_score_);
    int size();
    void set_file_name(std::string file_name_);
    std::string get_file_name();
    void set_normalized_true();
    void set_normalized_false();
    bool get_normalized();
    int get_seq_start();
    int get_seq_end();
    void get_lengths_from_fasta(std::string fasta_file_name);
    void read();
    GenomeMap get_genome_map(std::string mapped_name, std::string sample_name);
    void add_sam_file(SamFile* other);
    std::vector<std::string> get_entries();
    std::vector<std::string> get_headers();
    bool are_filters_equal(SamFile* other);
    void copy_filters_from_another(SamFile* other);
    std::vector<int> get_starts();
    std::vector<int> get_ends();
    std::vector<std::string> get_names();
    std::vector<float> get_normalized_scores();
};


enum direction {up, left, upper_left, none};

class PairwiseAligner{
    std::string algorithm = "global"; //TODO: switch to enum
    std::string query;
    std::string target;

    int gap_penalty = -8;

    float score;

    std::string query_aligned;
    std::string target_aligned;
    std::string alignment_string;

    std::vector<std::vector<int>> F;
    std::vector<std::vector<direction>> pointers;

    PairScore* ps;
public:
    PairwiseAligner();
    void set_algorithm(std::string alg);
    std::string get_algorithm();
    void set_query(std::string query_);
    std::string get_query();
    void set_target(std::string target_);
    std::string get_target();

    // alloc memory based on query and target size, choose correct algs, etc.
    void align();

    void needleman_wunsch();
    void traceback_nw();

    float get_score();
    std::string get_query_aligned();
    std::string get_target_aligned();
    std::string get_match_string();
    std::vector<std::vector<int>> get_F();
    std::vector<std::vector<direction>> get_pointers();
};

