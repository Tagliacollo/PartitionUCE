from functions import *
from Bio import AlignIO, SeqIO
import numpy as np

# open alignment
alignment = AlignIO.read(open('matrix_gen1.data'), "nexus")


# test the base frequencies function
p = bp_freqs_calc(alignment) # the alignment should includes 0.25 of each bases 
print(p)

e = entropy_calc(p)
print(e)


# test the entropy function
p1 = np.array([0.25, 0.25, 0.25, 0.25]) # should be an entropy 2
entropy_calc(p1)

p2 = np.array([0.0, 0.0, 0.5, 0.5]) # should be an entropy of 1
entropy_calc(p2)
