from Bio.Nexus import Nexus
from Bio import AlignIO
import numpy as np
from tqdm import tqdm
from collections import defaultdict

dataset_path = '/Users/Tagliacollo/Desktop/scripts_UCE_test/matrix.nex'
#dataset_path = '/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/test.nex'

def charset_uce_aln(aln):

    dat = Nexus.Nexus()
    dat.read(aln)
    aln = AlignIO.read(open(aln), "nexus")

    uce_aln  = []

    for name in tqdm(dat.charsets):
        sites = dat.charsets[name]
        start = min(sites)
        stop = max(sites) + 1
        # slice the alignment to get the UCE
        uce_aln.append(aln[:, start : stop])

    return(uce_aln)
