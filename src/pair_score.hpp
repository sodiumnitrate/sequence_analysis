#include <unordered_map>
#include <string>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>

std::unordered_map<std::string, int> read_scores(std::string& name);

class PairScore{
    std::string name;
    std::unordered_map<std::string, int> scores;
public:
    PairScore(std::string name_);

    void init_blosum50();
    void init_blosum62();
    int query(char res1, char res2);
};
