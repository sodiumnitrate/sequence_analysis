#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "include/sam_entry.hpp"

namespace py = pybind11;

// constructor
SamEntry::SamEntry(std::string& rn, std::string& s, std::string& m, int st, int ed, float as)
{
    read_name = rn;
    sequence = s;
    mapped_onto = m;
    start_pos = st;
    end_pos = ed;
    alignment_score = as;
}

std::string SamEntry::get_read_name(){ return read_name; }
std::string SamEntry::get_seq_str(){ return sequence; }
std::string SamEntry::get_mapped_onto(){ return mapped_onto; }
int SamEntry::get_start_pos(){ return start_pos; }
int SamEntry::get_end_pos(){ return end_pos; }
int SamEntry::get_length(){ return end_pos - start_pos + 1; }

Sequence SamEntry::to_sequence() {
    Sequence result(sequence);
    result.set_name(read_name);
    return result;
}

void init_sam_entry(py::module_ &m){
    py::class_<SamEntry>(m, "SamEntry", py::dynamic_attr())
        .def(py::init<std::string&, std::string&, std::string&, int, int, float >())
        .def("get_length", &SamEntry::get_length)
        .def("get_read_name", &SamEntry::get_read_name)
        .def("get_seq_str", &SamEntry::get_seq_str)
        .def("get_mapped_onto", &SamEntry::get_mapped_onto)
        .def("get_start_pos", &SamEntry::get_start_pos)
        .def("get_end_pos", &SamEntry::get_end_pos)
        .def("__len__", &SamEntry::get_length)
        ;
}
