import sys
import os
from setuptools import setup, find_packages
import sequence_analysis

setup(
        name='sequence_analysis',
        version=sequence_analysis.__version__,
        description='simple tools for analyzing biological sequences',
        url='https://github.com/sodiumnitrate/sequence_analysis.git',
        author='Irem Altan',
        author_email='irem.altan@yale.edu',
        license='',
        packages=find_packages(),
        install_requires=['numpy','matplotlib','biopython>=1.80','sklearn','colorama','networkx'],
        python_requires='>=3.6'
        #ext_modules=ext_modules
    )
