from Bio.Nexus import Nexus
from Bio import AlignIO, SeqIO, SeqUtils
import Bio
import numpy as np
import os, subprocess
from glob import glob
from itertools import combinations
from itertools import islice


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


def write_csv_windows(uce_dict, outfilename):
    '''
    write a csv file of a uce dictionsary
    where the keys are the names
    and the values are lists of something like entropies
    '''
    outfile = open(outfilename, 'w')
    outfile.write("name,start,stop\n")
    for key in uce_dict:
        window_tuple = uce_dict[key]

        # write that UCE to file
        outfile.write("%s,%d,%d\n" %(key, window_tuple[0], window_tuple[1]))

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

def sum_entropies_gc(list_entropies, list_gc):
def sitewise_entropies_and_gc(aln):

    norm_entropies = np.array([list_entropies])/max(list_entropies)
    norm_gc = np.array([list_gc])/max(list_gc)
    entropies = np.array([sitewise_entropies(aln)])
    gc = np.array([sitewise_gc(aln)])

    norm_entropies = entropies/np.max(entropies)
    norm_gc = gc/np.max(gc)

    sum_norm = norm_entropies + norm_gc

    return (sum_norm)


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


def get_all_wd(list_of_values):
    '''
    Input: MultipleSeqAlignment
    Output: list of tuples [ (start : end) ]
    k = minimum wd size 
    '''

    minimum_window_size = 50
    length = len(list_of_values)

    keep_windows = []

    # some lists of values are too small to split
    if length < 2*minimum_window_size:
        return ([(0, length)])

    for window in combinations(range(length), 2):
        start = window[0]
        stop = window[1]

        # left flank size
        if start < minimum_window_size:
            continue
        # right flank size
        if (length - stop) < minimum_window_size:
            continue
        # central window
        if (stop - start) < minimum_window_size:
            continue

        keep_windows.append(window)

    return (keep_windows)

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
        wd_core = get_sse(metric[ dat[0] : dat[1] ])
        wd_right = get_sse(metric[ dat[1] : ])

        res = float(wd_left + wd_core + wd_right)
    
        wd_values.append(res)

    return wd_values

def get_best_window(all_windows, values):

    all_sses = get_sum_sse_uce_partition(all_windows, values)

    min_sse = min(all_sses)
    
    # we could be smart and get all occurrences, then choose wisely
    # like this
    #min_indices = [i for i, x in enumerate(all_sses) if x == min_sse]

    # This just gets the first occurrence of the minimum value.
    best_window = all_windows[all_sses.index(min_sse)]

    return best_window


def get_best_windows(uce_dict):
    '''
    input: a dict of UCEs, where the keys are the names and the 
           values are lists of some metric
    output: a dict with UCE names as keys, and the best 
            tuple as the value
    '''

    best_windows = {}
    for key in uce_dict:
        values = uce_dict[key]
        all_windows = get_all_wd(values)
        best_window = get_best_window(all_windows, values)
        best_windows[key] = best_window
        print(key, len(values), best_window)

    return(best_windows)


def take(n, iterable):
    # taken from: http://stackoverflow.com/questions/7971618/python-return-first-n-keyvalue-pairs-from-dict
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


