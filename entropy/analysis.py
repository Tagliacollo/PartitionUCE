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
	

df = pd.DataFrame(uce_entropies['chr20_5193'], )
df.to_csv('uce_dat.csv')


