#include "pair_score.hpp"

PairScore::PairScore(std::string name_)
{
    name = name_;
    if (name.compare("blosum50") == 0) init_blosum50();
    else if (name.compare("blosum62") == 0) init_blosum62();
    else throw "unrecognized name";
}

std::unordered_map<std::string, int> read_scores(std::string& name){
    std::ifstream file("data/"+ name +".dat");
    std::string line;
    char res1, res2;
    std::string cols;

    std::unordered_map<std::string, int> result;

    if(!file.is_open()) throw "score matrix file not found";

    // header line
    std::getline(file, line);

    // cols
    std::getline(file, line);
    std::istringstream ss(line);
    char curr;
    for ( int i = 0; i < 25; i++){
        ss >> curr;
        cols.push_back(curr);
    }
    std::string pair;
    char second;
    int sij;
    while(std::getline(file, line)){
        ss.str(line);
        ss >> second;
        for(int i = 0; i < 25; i++){
            pair = second + cols[i];
            ss >> sij;
            result[pair] = sij;
        }
    }
    return result;
}

void PairScore::init_blosum50(){
    scores.clear();
    // load data and populate map
    scores = read_scores("blosum50");
}

void PairScore::init_blosum62(){
    scores.clear();
    scores = read_scores("blosum62");
}

int PairScore::query(char res1, char res2){
    if (scores.find(res1 + res2) == scores.end()) throw "char not found";
    return scores[res1 + res2];
}

