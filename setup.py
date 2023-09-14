import sys
import os
from pybind11 import get_cmake_dir
from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages
import sequence_analysis

cxx_std = int(os.environ.get("CMAKE_CXX_STANDARD", "14"))

ext_modules = [
    Pybind11Extension("sam_reader_cpp", ["src/sam_reader.cpp"], cxx_std=cxx_std),
    Pybind11Extension("fasta_reader_cpp", ["src/fasta_reader.cpp"], cxx_std=cxx_std),
]

setup(
        name='sequence_analysis',
        version=sequence_analysis.__version__,
        description='simple tools for analyzing biological sequences',
        url='https://github.com/sodiumnitrate/sequence_analysis.git',
        author='Irem Altan',
        author_email='irem.altan@yale.edu',
        license='',
        packages=find_packages(),
        install_requires=['numpy','matplotlib','biopython>=1.80','scikit-learn','colorama','networkx', 'pybind11'],
        python_requires='>=3.6',
        ext_modules=ext_modules
    )
