# sequence_analysis

This python module puts together some tools to do very basic (protein/dna/rna) sequence analysis. Some of these were created for research and some for my own education. Please use at your own risk. 

## Installation

First, make sure your python version is `>=3.10`. Then, make sure `wheel` is installed: `pip install wheel`. Then:

    git clone https://github.com/sodiumnitrate/sequence_analysis.git
    pip install ./sequence_analysis

## Basic usage

### Sequence manipulation

Creating and manipulating Sequence objects:

```python
from sequence_analysis import Sequence

# create seq
seq = Sequence("ACGT")

# change its name
seq.name = "test"

# the type (in this case 'dna') should be set automatically, and be printed with
print(seq.type)

# override type by (note that this can break future functionality, and should only be used if automatic detection doesn't work.)
seq.type = "protein"

# find index to a specific codon (assumes the current reading frame is valid)
idx = seq.find_codon("ACG")

# reverse complement (also assumes the current reading frame is valid)
seq_rc = seq.reverse_complement()

# translate 
prot = seq.translate()

# write single sequence to fasta file
prot.write_fasta("test.fasta")

# get open reading frames (assumes a minimum orf length of 87)
dna_seq = Sequence("AGTGCC...")
orfs = dna_seq.get_open_reading_frames()
# print protein sequence of one of the orfs
print(orfs[0].protein_sequence)
```
### Sequence sets

There is a SeqSet object that contains a vector of Sequence objects. More manipulations specific to more than one sequence can be done. (There are also future plans to allow for parallel manipulation of contained sequences through this object, such as tranlating all sequences in parallel.)

Create sequence sets:

```python
from sequence_analysis import SeqSet
from sequence_analysis import Sequence

sset = SeqSet()
s1 = Sequence("ACGT")
s2 = Sequence("GGCCTT")
# setting records also assigns type automatically
sset.records = [s1, s2]
```
Note that setting records does so by copying the sequences in memory, so further manipulation of `s1` or `s2` won't affect `sset`. 

Reading from and writing to `.fasta` files is possible:

```python
sset = SeqSet()
sset.read_fasta("example.fasta")
sset.write_fasta("new_file.fasta")
```
Alphabetizing and/or removing duplicates:

```python
# sort sequences such that they are in alphabetical order
sset.alphabetize()
# remove duplicates (does not require prior alphabetization, but results in alphabetized sequences)
sset.remove_duplicates()
```
Removing duplicates sorts the sequences in alphabetical order, and extends the name of each sequence such that the names of the removed, identical sequences are appended to the remaining sequence.

## Features to implement

This is still a work in progress. The goal is to implement the following features (in no particular order):

- ~~Move from `setuptools` to `CMake` to compile `C++` code~~
- Implement unit tests on the `C++` side within `CMake`
- Docs
- `OpenMP` under the hood for some functions (to work both on macOS & linux -- the macOS requirement is hard.)
- Parallelize some functionality on `python` side, with careful removal of GIL on `C++` side?
- Finish `PairwiseAligner`
- Clustering with pairwise distance calculation using `edlib`?
- The ability to stop execution with Ctrl+C
