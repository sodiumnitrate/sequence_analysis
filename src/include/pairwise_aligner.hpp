#include <vector>
#include <string>
#include "pair_score.hpp"

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
