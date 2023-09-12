/*
C++ function to make filtering while reading SAM files much faster.
(Hopefully)
*/

#include <cmath>
#include <pybind11/pybind11.h>
#include "pybind11/numpy.h"
#include <pybind11/stl.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <limits>

namespace py = pybind11;

py::list sam_reader(std::string file_name,
                    int start,
                    int end,
                    std::string mapped_onto,
                    float min_score)
{
    std::ifstream file (file_name);
    std::string line;
    std::vector<std::string> lines;

    // if no end is given, make it (practically) infinite
    if (end == -1){
        end = std::numeric_limits<int>::max();
    }

    if (!file.is_open())
    {
        // print error message
        std::cout << "WARNING: Failed to open file with name " << file_name << std::endl;
        // return empty list
        py::list py_list = py::cast(lines);
        return py_list;
    }

    std::string dummy, ref_name, seq, flag;
    int pos, length;
    float alignment_score;
    std::string as_flag = "AS:i:";

    int seq_start, seq_end;
    seq_start = std::numeric_limits<int>::max();
    seq_end = 0;

    std::cout << "Mapped onto: " << mapped_onto << std::endl;

    while(std::getline(file, line))
    {
        // skip headers
        if (line[0] == '@')
        {
            continue;
        }
        std::istringstream ss(line);

        ss >> dummy >> dummy >> ref_name >> pos >> dummy >> dummy >> dummy >> dummy >> dummy >> seq >> dummy;

        // check the name of reference read is mapped onto
        if (ref_name.compare(mapped_onto) != 0 && mapped_onto.size() != 0)
        {
            continue;
        }

        // check that there is any overlap between the requested range and the read
        length = seq.size();

        if (pos + length < start || pos > end)
        {
            continue;
        }

        while( ss >> flag)
        {
            // if flag is AS
            if ( as_flag.compare(0, 5, flag) )
            {
                alignment_score = std::stof(flag.substr(5, flag.length()));
                if (alignment_score < min_score)
                {
                    continue;
                }
            }
        }

        // passed all the tests, append to vector
        lines.push_back(line);
        seq_start = fmin(seq_start, pos);
        seq_end = fmax(seq_end, pos + length);

        if (seq_end > 20053868))
        {
            std::cout << ref_name << "," << pos << "," << length << std::endl;
        }
    }

    // add seq_start and seq_end 
    lines.push_back(std::to_string(seq_start) + "-" + std::to_string(seq_end));
    
    // return as python list
    py::list py_lines = py::cast(lines);
    return py_lines;
}

PYBIND11_MODULE(sam_reader_cpp, m)
{
    m.doc() = "C++ function to filter and read SAM file.";
    m.def("sam_reader", &sam_reader, "filters and reads SAM files.");
}