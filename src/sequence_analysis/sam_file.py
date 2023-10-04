"""
Python bindings for the SamFile object.
"""
from .sequence_analysis_cpp import SamFile as SamFile_cpp

import os
from concurrent.futures import ThreadPoolExecutor

class SamFile(SamFile_cpp):
    """
    Python bindings for the SamFile class.
    """

    def read_multiple_files(self, file_names):
        """
        Run the read function in parallel to read multiple files.
        """
        n_files = len(file_names)

        # TODO: a better way??
        if n_files > os.cpu_count():
            print("ERROR: you provided more files than the number of possible threads.")
            raise ValueError

        with ThreadPoolExecutor(n_files) as executor:
            results = executor.map(self.read_parallel, file_names)

        for result in results:
            self.add_sam_file(result)
            
    def read_parallel(self, file_name):
        sam = SamFile()
        sam.copy_filters_from_another(self)
        sam.read(file_name)
        return sam
