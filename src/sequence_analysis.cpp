/*
C++ module with pybind11 binding to do sequence analysis.

NOTE: you will need this: https://stackoverflow.com/a/36056804/2112406
TODO: make things pickleable and add copy & deepcopy support (https://pybind11.readthedocs.io/en/stable/advanced/classes.html#deepcopy-support)
TODO: sequence logo class by itself?
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "include/seq_set.hpp"
#include "include/sam_file.hpp"
#include "include/sam_entry.hpp"
#include "include/pairwise_aligner.hpp"
#include "include/genome_map.hpp"
#include "include/fasta_iterator.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

void init_sequence(py::module_ &);
void init_orf(py::module_ &);
void init_seq_set(py::module_ &);
void init_sam_file(py::module_ &);
void init_sam_entry(py::module_ &);
void init_pairwise_aligner(py::module_ &);
void init_genome_map(py::module_ &);
void init_fasta_iterator(py::module_ &);

PYBIND11_MODULE(sequence_analysis_cpp, m){
    init_sequence(m);
    init_orf(m);
    init_seq_set(m);
    init_sam_file(m);
    init_sam_entry(m);
    init_pairwise_aligner(m);
    init_genome_map(m);
    init_fasta_iterator(m);
    
#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
