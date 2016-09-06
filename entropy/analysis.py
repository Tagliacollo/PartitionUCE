from functions import *
from Bio import AlignIO, SeqIO
import numpy as np
import pandas as pd

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


####### DATA FRAME OF ENTROPIES

uce_dat = split_charsets_to_list('matrix.nex')

aln = uce_dat['gen1'] # TODO loop through each key in the dictionary uce_dat / problem: how to get the index instead of a key name. 

i = 0 
site_list = []
while i < aln.get_alignment_length():
	site = aln[:,i]
	bp_freqs = np.array([site.count(base)/len(site) for base in ['A', 'C', 'G', 'T']]) # my function bp_freqs_calc is giving me error messagens. I'll work on that. 
	value = entropy_calc(bp_freqs)
	site_list.append(value)
	i = i + 1

df = pd.DataFrame(site_list)
df.to_csv('uce_dat.csv')# test large datsets
dat = split_charsets_to_list('Crawford_2012.nex')

uce_dict = {}
for key in dat:

	aln = dat[key]
	site_list =[]
	uce_dict[key] = site_list
	
	i = 0 
		while i < aln.get_alignment_length():
		site = aln[:,i].replace('-','')
		bp = np.array([site.count(base)/len(site) for base in ['A', 'C', 'G', 'T']]) # my function bp_freqs_calc isn't working. The problem is the concatenation of each sequence in a long string
		val = entropy_calc(bp)
		site_list.append(val)
		i = i + 1

print(uce_dict)

#df = pd.DataFrame(uce_dict)
#df.to_csv('uce_dat.csv')
