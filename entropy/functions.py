from Bio.Nexus import Nexus
from Bio import AlignIO, SeqIO
import numpy as np

def split_charsets_to_list(matrix):
    '''
    INPUT: a nexus alignment with charsets
    OUTPUT: biopython generic aligments, one per charset et
    '''
    dat = Nexus.Nexus()
    dat.read(matrix)
    aln = AlignIO.read(open(matrix), "nexus")

    # TODO: check that the name has 'UCE' or 'uce' in it.

    aln_dict = {}
    for name in dat.charsets:
        sites = dat.charsets[name]
        aln_dict[name] = aln[:, min(sites):(max(sites) + 1)]

    return aln_dict


def bp_freqs_calc(alignment):
def bp_freqs_calc(aln):
    '''
    INPUT: biopython generic alignment
    OUTPUT: 1-D array of base frequencies
    '''
    aln = ''
    for line in alignment:
        aln += str(line.seq)
    
    bp_freqs = np.array([aln.count(base)/len(aln) for base in ['A', 'C', 'G', 'T']])
    one_str = ""
    for seq in aln:
        one_str += seq

    seq = one_str.upper()

    A = seq.seq.count('A') 
    T = seq.seq.count('T')
    G = seq.seq.count('G')
    C = seq.seq.count('C')

    sum_count = A + T + G + C

    bp_freqs = np.array([ A/sum_count, T/sum_count, G/sum_count, C/sum_count])
    
    return bp_freqs
    return (bp_freqs)

def entropy_calc(p):
    '''
    INPUT: 1D array of base frequencies
    OUTPUT: estimates of entropies 
    
    copied from: http://nbviewer.ipython.org/url/atwallab.cshl.edu/teaching/QBbootcamp3.ipynb
    '''
    p = p[p!=0] # modify p to include only those elements that are not equal to 0
    return np.dot(-p,np.log2(p)) # the function returns the entropy result
