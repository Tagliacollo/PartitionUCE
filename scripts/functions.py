from Bio.Nexus import Nexus
from Bio import AlignIO, SeqIO, SeqUtils
import Bio
import numpy as np
import os, subprocess
from glob import glob
from itertools import combinations


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

def write_phylip(aln, aln_path):
    '''
    I look forward to the day I don't have to write
    custom alignment export functions
    '''
    aln_handle = open(aln_path, "w")
    header = "%d\t%d\n" %(len(aln), aln.get_alignment_length())
    aln_handle.write(header)

    for s in aln:
        aln_handle.write("%s    %s\n" %(s.name, str(s.seq.upper())))        
    aln_handle.close()

def sitewise_TIGER(aln, tigger_path):
    '''
    input: biopython generic alignment, and path to tigger binary
    output: list of tigger rates
    '''

    # we'll use this as the temp dir...
    tigger_dir = os.path.dirname(tigger_path)
    aln_path = os.path.join(tigger_dir, "aln.phy")

    write_phylip(aln, aln_path)    

    subprocess.call([tigger_path, aln_path], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    tigger_output = ''.join([aln_path.rstrip("phy"), "tigger"])

    with open(tigger_output) as f:
        lines = f.read().splitlines()

    tiggers = [float(l) for l in lines]

    # clean up
    os.remove(tigger_output)
    os.remove(aln_path)

    return(tiggers)


def sitewise_entropies(aln):

    entropies = []
    for i in range(aln.get_alignment_length()):

        site_i = aln[:,i:i+1]
        ent_i = alignment_entropy(site_i)
        entropies.append(ent_i)

    return (entropies)

def sitewise_gc(aln):

    gc = []
    for i in range(aln.get_alignment_length()):
        site_i = aln[:,i]
        gc_i = SeqUtils.GC(site_i)
        gc.append(gc_i)

    return (gc)



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


def nexus_concat(dataset_path):
    '''
    INPUT: The path where all nexus files are allocated.
    ex. os.chdir('/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/Faircloth_2013')

    OUTPUT: nexus supermatrix already including charsets 
    '''

    nexus_list = glob('*.nex')
    mtx = [(nex, Nexus.Nexus(nex)) for nex in nexus_list]
    supermtx = Nexus.combine(mtx)

    return supermtx


def get_sse(metric):
    '''
    input: list of values
    output: sum of squared errors
    '''
    sse = np.sum((metric - np.mean(metric))**2)

    return sse


def get_all_wd(aln, k):
    '''
    Input: MultipleSeqAlignment
    Output: list of tuples [ (start : end) ]
    k = minimum wd size 
    '''
    aln_size = aln.get_alignment_length()
    
    if type(aln) != Bio.Align.MultipleSeqAlignment:
        print ('alignment has to be of the class Bio.Align.MultipleSeqAlignment')

    elif aln_size <= 1:
        print ('alignment is too small: %i bp' %(aln_size))

    else: 
        res = list(combinations(range(aln_size), 2))

        # exclude wd that gives only two partition
        num = [ val for val in res if val[0] != 0 and val[1] != (aln_size) ]
    
        # set a minimum wd size (k)
        out = [ i for i in num if i[1] - i[0] > k ]

        return out

def get_sum_sse_uce_partition(tuple_list, metric):
    '''
    input: 
        1) a tuple list from get_all_wd
        2) a list including values of a metric (eg. entropies) 
    output: list of sse values for each partitioning scheme
    '''
    wd_values = []
    for dat in tuple_list:

        wd_left = get_sse(metric[  : dat[0]])
        wd_core = get_sse(metric [ dat[0] : dat[1] ])
        wd_right = get_sse(metric [ dat[1] : ])

        res = float(wd_left + wd_core + wd_right)
    
        wd_values.append(res)

    return wd_values



