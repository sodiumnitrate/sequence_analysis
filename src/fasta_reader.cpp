/*
C++ function to make reading fasta faster, especially when filtering by sequence name.
*/
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

py::list fasta_reader(std::string file_name,
                      std::vector<std::string> sequence_names)
{
    std::ifstream file (file_name);
    std::string line;
    std::string name;
    std::string seq;
    std::vector<std::string> lines;
    
    if (!file.is_open())
    {
        // print error message
        std::cout << "WARNING: Failed to open file with name " << file_name << std::endl;    
        // return empty list
        py::list py_list = py::cast(lines);
        return py_list;
    }

    // create unordered set from sequence names (O(N))
    std::unordered_set<std::string> seq_names;
    for(auto i = 0; i < sequence_names.size(); i++)
    {
        seq_names.insert(sequence_names[i]);
    }

    // whether or not any sequence names are provided for filtering
    bool filter = (sequence_names.size() > 0) ? true : false;

    int seq_ct = 0;
    while(std::getline(file, line))
    {
        if (line[0] == '>')
        {
            // we are at a new sequence
            if (seq_ct > 0)
            {
                // add only if filter is empty or name in the filter list
                if (!filter || seq_names.find(name) != seq_names.end())
                {
                    lines.push_back(name);
                    lines.push_back(seq);
                }
                // reset seq string
                seq = "";               
            }
            // get new name
            name = line.substr(1, line.length()-1);
            seq_ct++;
        }
        else
        {
            // add new line to seq
            seq += line;
        }
    }

    // add the final one
    if (!filter || seq_names.find(name) != seq_names.end())
    {
        lines.push_back(name);
        lines.push_back(seq);
    }

    // return the result to python
    py::list py_list = py::cast(lines);
    return py_list;
}

PYBIND11_MODULE(fasta_reader_cpp, m)
{
    m.doc() = "C++ function to read fasta files while filtering sequence names.";
    m.def("fasta_reader", &fasta_reader, "filters and reads fasta files.");
}