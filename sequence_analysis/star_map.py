"""
Class for running STAR for illumina reads under the hood.

IMPORTANT: use of this module requires having MMseqs2 installed and
in your system path as a command line tool.

Also, it assumes we are trying to align short reads, so STAR itself
is sufficient.

TODO: python bindings?
"""
import subprocess
import os
import shutil

# STAR --runMode alignReads --readMapNumber -1 -readNameSeparator _ --runThreadN 18 --genomeDir ../STAR --readFilesIn ../../squid_transcriptome_from_genome_paper/SRR18071790_1.fastq

# STAR --runMode genomeGenerate --genomeDir . --genomeFastaFiles ../../squid_genome/GCA_023376005.1_UCB_Dpea_1_genomic.fna --runThreadN 18

class StarMap:
    def __init__(self, n_threads=1):
        self.n_threads = n_threads
        self.genome_dir = None
        self.genome_indexed = False
        self.query_file = None

        self.exec_str = None

        self.check_for_STAR()

    def check_for_STAR(self):
        """
        We have to make sure that STAR is installed and is in the
        system path.
        """
        path = shutil.which("STAR")
        if path is None:
            print("ERROR: STAR not found.")
            raise OSError

    def index_genome(self, genome_fasta, genome_dir):
        """
        Given a .fasta file for the genome, generate indices.
        """
        exec_string = f"STAR --runMode genomeGenerate --genomeDir {genome_dir} --genomeFastaFiles {genome_fasta} --runThreadN {self.n_threads}"

        self.exec_str = exec_string

        result = subprocess.run([exec_string], shell=True, capture_output=True, text=True)
        self.genome_dir = genome_dir

    def set_genome_dir(self, genome_dir):
        """
        If the genome is already indexed, set its directory.
        """
        self.genome_dir = genome_dir

        # TODO: some checks to see if indices are actually there

    def align_reads(self):
        """
        Align reads given a query fasta file and genome dir with indices.
        """
        if self.genome_dir is None:
            print("ERROR: genomeDir not found. Please set it.")
            raise ValueError

        exec_str = f"STAR --runMode alignReads --readMapNumber -1 --runThreadN {self.n_threads} --genomeDir {self.genome_dir} --readFilesIn {self.query_file}"

        result = subprocess.run([exec_str], shell=True, capture_output=True, text=True)

        self.exec_str = exec_str