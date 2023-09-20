import sys
import os
from glob import glob
#from pybind11 import get_cmake_dir
from pybind11.setup_helpers import Pybind11Extension
from setuptools import setup, find_packages
import sequence_analysis

cxx_std = int(os.environ.get("CMAKE_CXX_STANDARD", "20"))

ext_modules = [
    Pybind11Extension("sequence_analysis_cpp", sorted(glob("src/*.cpp")), cxx_std=cxx_std)
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
        install_requires=['numpy','matplotlib','colorama','networkx', 'pybind11', 'biopython', 'pydot', 'pytest'],
        python_requires='>=3.10',
        ext_modules=ext_modules
    )
