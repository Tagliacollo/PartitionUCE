from Bio.Nexus import Nexus
from Bio import AlignIO, SeqIO
import numpy as np
import os

def write_csv(uce_dict, outfilename, parameter_name):
    '''
    write a csv file of a uce dictionsary
    where the keys are the names
    and the values are lists of something like entropies
    '''
    outfile = open(outfilename, 'w')
    outfile.write("name,site,%s\n" %(parameter_name))
    for key in uce_dict:
        ent = uce_dict[key]

        # define middle site as zero
        middle = int(float(len(ent)) / 2.0)
        sites = np.array(range(len(ent))) - middle

        names = [key] * len(ent)

        # write that UCE to file
        for i in range(len(ent)):
            outfile.write("%s,%d,%f\n" %(names[i], sites[i], ent[i]))

    outfile.close()


def alignment_entropy(aln):
    '''
    input: biopython generic alignment
    output: entropy

    '''
    bp_freqs = bp_freqs_calc(aln)
    entropy = entropy_calc(bp_freqs)
    return (entropy)


def sitewise_TIGER(aln, tigger_path):
    '''
    input: biopython generic alignment, and path to tigger binary
    output: list of tigger rates
    '''

    # write alignment as phylip file

    # get path of phylip file


    os.system("%s %s" %(tigger_path, aln_path))
    tigger_output = ''.join(aln_path.rstrip("phy"), "tigger")

    with open(tigger_output) as f:
        lines = f.read().splitlines()

    tiggers = [float(l) for l in lines]

    os.remove(tigger_output)
    os.remove(alignment_file)

    return(tiggers)


def sitewise_entropies(aln):

    entropies = []
    for i in range(aln.get_alignment_length()):

        site_i = aln[:,i:i+1]
        ent_i = alignment_entropy(site_i)
        entropies.append(ent_i)

    return (entropies)


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


def bp_freqs_calc(aln):
    '''
    INPUT: biopython generic alignment
    OUTPUT: 1-D array of base frequencies
    '''
    one_str = ""
    for seq in aln:
        one_str += seq

    seq = one_str.upper()

    A = seq.seq.count('A') 
    T = seq.seq.count('T')
    G = seq.seq.count('G')
    C = seq.seq.count('C')

    sum_count = A + T + G + C

    bp_freqs = np.array([ A, C, G, T])/float(sum_count)
    
    return (bp_freqs)


def entropy_calc(p):
    '''
    INPUT: 1D array of base frequencies
    OUTPUT: estimates of entropies 
    
    copied from: http://nbviewer.ipython.org/url/atwallab.cshl.edu/teaching/QBbootcamp3.ipynb
    '''
    p = p[p!=0] # modify p to include only those elements that are not equal to 0
    return np.dot(-p,np.log2(p)) # the function returns the entropy result
