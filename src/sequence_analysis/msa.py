"""
Hold the multiple sequence alignment class. Runs using the mafft executable.

MAFFT: https://mafft.cbrc.jp/alignment/software/

Requires mafft to be installed such that mafft can be run in the command line.

TODO:
- implement options
- custom matrix?
- make compatible on windows
- is there any other alignment software I can use that has an API of sorts?
"""

from shutil import which
import subprocess
import os
import time

from sequence_analysis import SeqSet

class MSA:
    """
    MSA object holds info about sequences to be aligned, alignment options, and alignment results.
    """
    def __init__(self, sequences, scratch_folder_name="tmp"):
        self.sequences = sequences
        self.check_input_data()

        self.mafft = None
        self.check_for_mafft()

        self.input_name = None
        self.output_name = None

        self.scratch_folder_name = scratch_folder_name

        # gap opening penalty
        self.op = 1.53
        # gap extension penalty
        self.ep = 0.123

        # to be populated after alignment
        self.alignment_info = None
        self.aligned_sequences = None
        self.elapsed = None

    def check_input_data(self):
        """
        Check that sequences are input as a SeqSet object.
        TODO: make it possible for input to be a fasta filename as well?
        """
        if not isinstance(self.sequences, SeqSet):
            print("ERROR: sequences must be a SeqSet object.")
            raise TypeError

    def check_for_mafft(self):
        """
        Check that mafft is installed and is accessible.
        """
        if which('mafft') is None:
            self.mafft = False
        else:
            self.mafft = True

    def default_output_name(self, aligned=False, overwrite=False):
        """
        Generate .fasta file names to be used to align using mafft.
        """
        if aligned:
            file_name = "tmp_aligned.fasta"
        else:
            file_name = "tmp.fasta"

        if not overwrite:
            ct = 0
            while file_name in os.listdir(self.scratch_folder_name):
                if aligned:
                    file_name = f"tmp_{ct}_aligned.fasta"
                else:
                    file_name = f"tmp_{ct}.fasta"
                ct += 1

        if aligned:
            self.output_name = file_name
        else:
            self.input_name = file_name

    def align(self, overwrite=False, clean=False):
        """
        Align sequences. Save info about alignment, etc.

        If clean=True, remove tmp files so that the aligned sequences live within
        the MSA object, in memory only.
        """
        if not self.mafft:
            print("ERROR: mafft either not installed, or executable not in path.")
            # TODO: is this the correct way?
            raise OSError

        if self.scratch_folder_name[-1] == '/':
            self.scratch_folder_name = self.scratch_folder_name[:-1]

        created = False
        if not os.path.exists(self.scratch_folder_name):
            created = True
            os.mkdir(self.scratch_folder_name)

        begin = time.time()
        if self.output_name is None:
            self.default_output_name(overwrite=overwrite, aligned=True)

        if self.input_name is None:
            self.default_output_name(overwrite=overwrite, aligned=False)

        # write sequences as a fasta file
        self.sequences.write_fasta(self.scratch_folder_name + '/' + self.input_name)

        run_str = f"mafft --op {self.op} --ep {self.ep} {self.scratch_folder_name + '/' + self.input_name} >  {self.scratch_folder_name + '/' + self.output_name}"

        result = subprocess.run([run_str], shell=True, capture_output=True, text=True, check=True)

        # not sure why mafft outputs to stderr instead of stdout, but here we are
        self.alignment_info = result.stderr
        self.aligned_sequences = SeqSet()
        self.aligned_sequences.read_fasta(self.scratch_folder_name + '/' + self.output_name)
        
        if clean:
            os.remove(self.scratch_folder_name + '/' + self.output_name)
            os.remove(self.scratch_folder_name + '/' + self.input_name)
            if created:
                os.rmdir(self.scratch_folder_name)

        end = time.time()

        self.elapsed = end - begin
