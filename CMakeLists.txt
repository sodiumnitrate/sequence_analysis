cmake_minimum_required(VERSION 3.15...3.26)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

project(
    ${SKBUILD_PROJECT_NAME}
    VERSION ${SKBUILD_PROJECT_VERSION}
    LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

python_add_library(sequence_analysis_cpp MODULE src/genome_map.cpp src/orf.cpp src/pairwise_aligner.cpp src/sam_file.cpp src/seq_set.cpp src/sequence.cpp src/sequence_analysis.cpp WITH_SOABI)
target_link_libraries(sequence_analysis_cpp PRIVATE pybind11::headers)
target_compile_definitions(sequence_analysis_cpp PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS sequence_analysis_cpp DESTINATION sequence_analysis)
