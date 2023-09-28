from .sequence_analysis_cpp import FastaIterator as FastaIterator_cpp

class FastaIterator(FastaIterator_cpp):
    def read(self):
        """
        Read one by one. Return None if we reach the end of file.

        // TODO: a more robust way to signal EOF from C++ side?
        // TODO: how to make this a true blue python iterator?
        """
        seq = self.get_next()
        if seq.seq_str == "":
            return None

        return seq
