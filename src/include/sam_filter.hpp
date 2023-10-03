#include <unordered_map>
#include <unordered_set>
#include <string>
#include <vector>

class SamFilter{
    std::unordered_set<std::string> nameset;
    std::unordered_map<std::string, std::vector<std::tuple<int, int> > > nuc_ranges;
    bool contains_empty = false;
public:
    SamFilter(std::vector<std::string>& names, std::vector<int>& starts, std::vector<int>& ends);
    bool query(std::string& name, int start, int end);
    bool check_name(std::string& name);
};