from functions import *
from Bio import AlignIO, SeqIO
import numpy as np
import pandas as pd


# split alignment into dict of charsets
uce_alns = split_charsets_to_list('Crawford_2012.nex')

# get entropies for each charset
uce_entropies = {}
for key in uce_alns:
	uce_entropies[key] = sitewise_entropies(uce_alns[key])
	

# write csv
outfile = open('Crawford_2012.csv', 'w')
outfile.write("name,site,entropy\n")
for key in uce_entropies:
	ent = uce_entropies[key]
	print(key)
	# define middle site as zero
	middle = int(float(len(ent)) / 2.0)
	sites = np.array(range(len(ent))) - middle

	names = [key] * len(ent)

	# write file
	for i in range(len(ent)):
		outfile.write("%s,%d,%f\n" %(names[i], sites[i], ent[i])

	outfile.close()

