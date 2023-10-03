#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <cctype>
#include <unordered_map>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <sstream>
#include "include/orf.hpp"

namespace py = pybind11;

// ORF class
OpenReadingFrame::OpenReadingFrame() {};

void OpenReadingFrame::set_props(std::string rna_sequence_, std::string protein_sequence_, int start_, int stop_, int strand_, int frame_){
    rna_sequence = rna_sequence_;
    protein_sequence = protein_sequence_;
    start = start_;
    stop = stop_;
    strand = strand_;
    frame = frame_;
}
void OpenReadingFrame::set_start(int start_) {start = start_;}
int OpenReadingFrame::get_start() { return start;}
void OpenReadingFrame::set_stop(int stop_) { stop = stop_;}
int OpenReadingFrame::get_stop() {return stop;}
void OpenReadingFrame::set_strand(int strand_) { strand = strand_;}
int OpenReadingFrame::get_strand() {return strand;}
void OpenReadingFrame::set_frame(int frame_) {frame = frame_;}
int OpenReadingFrame::get_frame() { return frame; }
std::string OpenReadingFrame::get_rna_sequence() { return rna_sequence;}
void OpenReadingFrame::set_rna_sequence(std::string rs) {rna_sequence = rs;}
std::string OpenReadingFrame::get_protein_sequence() { return protein_sequence;}
void OpenReadingFrame::set_protein_sequence(std::string ps) {protein_sequence = ps;}

void init_orf(py::module_ &m){
    py::class_<OpenReadingFrame>(m, "OpenReadingFrame", py::dynamic_attr())
        .def(py::init<>())
        .def("set_start", &OpenReadingFrame::set_start)
        .def("get_start", &OpenReadingFrame::get_start)
        .def("set_stop", &OpenReadingFrame::set_stop)
        .def("get_stop", &OpenReadingFrame::get_stop)
        .def("set_strand", &OpenReadingFrame::set_strand)
        .def("set_strand", &OpenReadingFrame::get_strand)
        .def("set_frame", &OpenReadingFrame::set_frame)
        .def("get_frame", &OpenReadingFrame::get_frame)
        .def("set_rna_sequence", &OpenReadingFrame::set_rna_sequence)
        .def("get_rna_sequence", &OpenReadingFrame::get_rna_sequence)
        .def("set_protein_sequence", &OpenReadingFrame::set_protein_sequence)
        .def("get_protein_sequence", &OpenReadingFrame::get_protein_sequence)
        .def_property("start", &OpenReadingFrame::get_start, &OpenReadingFrame::set_start)
        .def_property("stop", &OpenReadingFrame::get_stop, &OpenReadingFrame::set_stop)
        .def_property("strand", &OpenReadingFrame::get_strand, &OpenReadingFrame::set_strand)
        .def_property("frame", &OpenReadingFrame::get_frame, &OpenReadingFrame::set_frame)
        .def_property("rna_sequence", &OpenReadingFrame::get_rna_sequence, &OpenReadingFrame::set_rna_sequence)
        .def_property("protein_sequence", &OpenReadingFrame::get_protein_sequence, &OpenReadingFrame::set_protein_sequence)
        .def("__repr__",
             [](OpenReadingFrame &a){
                return "<sequence_analysis.OpenReadingFrame object>";
             })
        ;        
}
