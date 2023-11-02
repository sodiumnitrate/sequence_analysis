"""
Tree builder that runs IQTREE.

Assumes iqtree2 is available as a command line tool.

TODO: add more of the possible options
"""

from shutil import which
import os
import subprocess

from sequence_analysis import SeqSet

class TreeBuilder:
    """
    TreeBuilder object holds info about aligned sequences to build trees from.
    """
    def __init__(self, sset, output_folder="", name=None):
        self.sset = sset
        self.output_folder = output_folder
        self.name = name

        self.iqtree = None
        self.check_for_iqtree()

        self.bootstrap_replicates = 1000

        self.run_info = None
        self.run_stderr = None
        self.run_result = None

    def check_for_iqtree(self):
        """
        Check that iqtree is installed and is accessible.
        """
        if which('iqtree2') is None:
            self.iqtree = False
        else:
            self.iqtree = True

    def preprocess(self):
        """
        Trim gaps, remove duplicates.
        """
        self.sset.trim_gaps()
        self.sset.remove_duplicates()

        if len(self.sset) <= 5:
            raise ValueError

    def build(self, preprocess=True):
        """
        Run tree builder.
        """
        if not self.iqtree:
            raise OSError

        if preprocess:
            self.preprocess()

        if not os.path.exists(self.output_folder) and self.output_folder != "":
            # TODO: correct exception?
            raise FileExistsError

        cwd = os.getcwd()
        if self.output_folder != "":
            os.chdir(self.output_folder)
        if self.name is None:
            name = "sequences.fasta"
        else:
            name = self.name
        self.sset.write_fasta(name)

        run_str = f"iqtree2 -s {name} -alrt {self.bootstrap_replicates} -B {self.bootstrap_replicates} -T AUTO"
        result = subprocess.run([run_str], shell=True, capture_output=True, text=True, check=True)

        self.run_str = run_str
        self.run_info = result.stdout
        self.run_stderr = result.stderr
        self.run_result = result

        os.chdir(cwd)