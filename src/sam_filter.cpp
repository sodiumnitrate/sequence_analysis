#include <vector>
#include <string>
#include <tuple>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <cmath>
#include <iostream>
#include "include/sam_filter.hpp"

SamFilter::SamFilter(std::vector<std::string>& names, std::vector<int>& starts, std::vector<int>& ends){
    // check input sizes
    if (!(names.size() == starts.size() && starts.size() == ends.size())) throw "size mismatch in input";

    unsigned int len = names.size();

    // process names
    for (const auto& name : names){
        // if (name.compare("") != 0){
            nameset.insert(name);
            std::vector<std::tuple<int, int> > empty;
            nuc_ranges[name] = empty;
        // }
        if (name.compare("") == 0) contains_empty = true;
    }

    std::string name;
    int s_idx, e_idx;
    for (unsigned i = 0; i < len; i++){
        name = names[i];
        s_idx = starts[i];
        e_idx = ends[i];

        if (e_idx == -1){
            e_idx = std::numeric_limits<int>::max();
        }

        if (name.compare("") == 0){
            for (auto t : nameset){
                nuc_ranges[t].push_back(std::make_tuple(s_idx, e_idx));
            }
        }
        else{
            nuc_ranges[name].push_back(std::make_tuple(s_idx, e_idx));
        }
    }

    // merge intervals
    for (auto name : nameset){
        unsigned int len = nuc_ranges.at(name).size();
        std::vector<int> visited;
        for (unsigned int i = 0; i < len; i++) visited.push_back(0);

        // adjacency matrix
        std::vector<std::vector<int> > adjacency(len);
        int s1, e1, s2, e2;
        for (unsigned int i = 0; i < len; i++){
            s1 = std::get<0>(nuc_ranges.at(name)[i]);
            e1 = std::get<1>(nuc_ranges.at(name)[i]);
            for ( unsigned int j = i + 1; j < len; j++){
                s2 = std::get<0>(nuc_ranges.at(name)[j]);
                e2 = std::get<1>(nuc_ranges.at(name)[j]);
                if(!(e1 < s2 || e2 < s1)){
                    adjacency[i].push_back(j);
                    adjacency[j].push_back(i);
                }
            }
        }

        // bfs
        int curr_idx, min_start, max_end;
        for (auto& t : visited) t = 0;
        std::vector<std::tuple<int, int> > curr;
        for(int i = 0; i < len; i++){
            if(visited[i] == 1) continue;
            min_start = std::get<0>(nuc_ranges.at(name)[i]);
            max_end = std::get<1>(nuc_ranges.at(name)[i]);

            std::queue<int> Q;
            Q.push(i);
            while(!Q.empty()){
                curr_idx = Q.front();
                Q.pop();
                visited[curr_idx] = 1;
                min_start = fmin(min_start, std::get<0>(nuc_ranges.at(name)[curr_idx]));
                max_end = fmax(max_end, std::get<1>(nuc_ranges.at(name)[curr_idx]));
                for(auto t: adjacency[curr_idx]){
                    if (visited[t] == 0) Q.push(t);
                }
            }
            curr.push_back(std::make_tuple(min_start, max_end));
        }
        nuc_ranges[name] = curr;
    }
}

bool SamFilter::query(std::string& name, int start, int end){
    std::string curr_name = name;
    if(nameset.find(name) == nameset.end()){
        if (!contains_empty) return false;
        else curr_name = "";
    }
    int s, e;
    bool works = false;
    for (auto t : nuc_ranges.at(curr_name)){
        s = std::get<0>(t);
        e = std::get<1>(t);

        if(!(start > e || end < s)) works = true;
    }
    return works;
}

bool SamFilter::check_name(std::string& name){
    if (contains_empty) return true;
    if (nameset.find(name) == nameset.end()) return false;
    return true;
}
