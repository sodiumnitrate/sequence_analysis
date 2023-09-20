#include <string>
#include <vector>

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
    OpenReadingFrame();
    void set_props(std::string rna_sequence_, std::string parent_sequence_, std::string protein_sequence_, int start_, int stop_, int strand_, int frame_);
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
    std::string get_parent_sequence();
    void set_parent_sequence(std::string ps);
    std::string get_protein_sequence();
    void set_protein_sequence(std::string ps);
};

// Sequence class
class Sequence {
    std::string name;
    std::string seq_str;
    std::string type;
public:
    Sequence(std::string &seq_str_);
    void set_name(std::string &name_);
    std::string get_name();
    void set_seq(std::string &seq_str_);
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
    void write_fasta(std::string file_name);

    // get orfs
    std::vector<OpenReadingFrame> get_open_reading_frames(unsigned int min_len);
};

class SeqSet {
    std::string name;
    // pointers to sequences!
    std::vector<Sequence> records;
    std::string type;
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

class SamFile{
    // path to the .sam file
    std::string file_name;

    // filter options
    std::vector<int> start = {0}; //for each name in mapped onto, give start and end
    std::vector<int> end = {1};
    std::vector<std::string> mapped_onto = {""};
    float min_score = 0;

    // if AS is not normalized, we need to use sequence lengths
    bool normalized_score = true;
    std::unordered_map<std::string, unsigned int> lengths;

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
};

class PairwiseAligner{
    std::string algorithm = "local"; //TODO: switch to enum
    std::string query;
    std::string target;

    float score;
public:
    PairwiseAligner();
    void set_algorithm(std::string alg);
    std::string get_algorithm();
    void set_query(std::string query_);
    std::string get_query();
    void set_target(std::string target_);
    std::string get_target();

    float get_score();
};

class GenomeMap{
    std::string chromosome_name;
    std::vector<unsigned int> heatmap;
    int heatmap_start;
    int heatmap_end;
public:
    GenomeMap();
}
