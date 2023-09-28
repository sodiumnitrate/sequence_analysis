"""
This class runs mmseqs easy-search in the bg to search
for a sequence (or list of sequences) within another list of sequences.

IMPORTANT: use of this module requires having MMseqs2 installed
and in your system path as a command line tool.

see: https://github.com/soedinglab/MMseqs2

TODO: python bindings?
TODO: wrap tmp in a python tempdir and delete after exec
"""
import shutil
import subprocess
import os

class EasySearch:
    def __init__(self, query_file, reference_file):
        self.query_file = query_file
        self.reference_file = reference_file
        if not os.path.isfile(self.query_file):
            print(f"ERROR: {self.query_file} can't be found.")
            raise ValueError

        if not os.path.isfile(self.reference_file):
            print(f"ERROR: {self.reference_file} can't be found.")
            raise ValueError

        # check that we have mmseqs2 installed
        self.mmseqs2 = None
        self.check_for_mmseqs2()

        self.search_type = None
        self.tmp_folder_name = "tmp"
        self.max_e_value = None
        self.cov_mode = None
        # sam files by default
        self.format_mode = 1
        self.output_name = None

        self.exec_string = None

    def check_for_mmseqs2(self):
        """
        We have to make sure that mmseqs2 is installed and is in the
        system path.
        """
        path = shutil.which("mmseqs")
        if path is None:
            self.mmseqs2 = False
            return

        result = subprocess.run(["mmseqs"], shell=True, capture_output=True, text=True)
        
        if 'MMseqs2' not in result.stdout:
            self.mmseqs2 = False
            return

        self.mmseqs2 = True


    def set_search_parameters(self,
                              search_type=3,
                              output_name="out.sam",
                              max_e_value=1e-16,
                              cov_mode=2,
                              tmp_folder_name="tmp"):
        """
        Set search parameters.

        Only allow for SAM output for now (format_mode=1).
        """

        self.search_type = search_type

        # check if output_name exists already
        if os.path.isfile(output_name):
            print(f"ERROR: file {output_name} already exists.")
            raise ValueError

        self.output_name = output_name
        self.max_e_value = max_e_value
        self.cov_mode = cov_mode
        self.tmp_folder_name = tmp_folder_name

    def run_search(self):
        """
        Given the options, run mmseqs2.
        """
        if not self.mmseqs2:
            print("ERROR: mmseqs not found.")
            # is this the right exception to raise?
            raise OSError

        if any([self.search_type is None, self.max_e_value is None, self.cov_mode is None, self.output_name is None]):
            print("WARNING: setting default search parameters.")
            self.set_search_parameters()

        # construct the command to execute
        self.exec_string = f"mmseqs easy-search --search-type {self.search_type}"
        self.exec_string += f" {self.query_file} {self.reference_file}"
        self.exec_string += f" {self.output_name} {self.tmp_folder_name} -e {self.max_e_value}"
        self.exec_string += f" --cov-mode {self.cov_mode} --format-mode {self.format_mode}"

        result = subprocess.run([self.exec_string], shell=True)


        # TODO: why doesn't the following work?
        """
        popen = subprocess.Popen(self.exec_string, stdout=subprocess.PIPE, universal_newlines=True)
        for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line
        popen.stdout.close()
        return_code = popen.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, self.exec_string)
        """