from Bio.Nexus import Nexus
from Bio import AlignIO, SeqIO, SeqUtils
import Bio
import numpy as np
import os, subprocess
from glob import glob
from itertools import combinations
from itertools import islice


def process_dataset(dataset_path, metrics, weights, outfilename):
    '''
    input: dataset_path: path to a nexus alignment with UCE charsets
           metrics: a list of 'gc', 'entropy' or both
           weights: a 1D np.array with weights for entropy, gc, multi *in that order*
    output: csv files written to disk
    '''
    outfile = open(outfilename, 'w')
    outfile.write("name,uce_site,aln_site,window_start,window_stop,%s\n" %(','.join(metrics)))
    outfile.close()

    dat = Nexus.Nexus()
    dat.read(dataset_path)
    aln = AlignIO.read(open(dataset_path), "nexus")

    for name in dat.charsets:
        print name
        sites = dat.charsets[name]
        start = min(sites)
        stop = max(sites) + 1
        # slice the alignment to get the UCE
        uce_aln = aln[:, start:stop]

        best_window, metric_array = process_uce(uce_aln, metrics, weights)

        write_csvs(best_window, metric_array, sites, name, outfilename)

def write_csvs(best_window, metrics, aln_sites, name, outfilename):

    N = len(aln_sites)
    middle = int(float(N) / 2.0)
    uce_sites = np.array(range(N)) - middle

    outfile = open(outfilename, 'a')

    names = [name]*N
    window_start = [best_window[0]]*N
    window_stop = [best_window[1]]*N
    metrics = metrics.tolist()

    all_info = [names, uce_sites, aln_sites, window_start, window_stop]

    for m in metrics:
        all_info.append(m)

    all_info = zip(*all_info)

    for i in all_info:
        line = [str(thing) for thing in i]
        line = ','.join(line)
        outfile.write(line)
        outfile.write("\n")
    outfile.close()


def process_uce(aln, metrics):
        
    windows = get_all_windows(aln)
    
    entropy = sitewise_entropies(aln)
    gc      = sitewise_gc(aln)
    multi   = sitewise_multi(aln)

    metrics = np.array([entropy, gc, multi])

    best_window = get_best_window(metrics, windows)

    return(best_window, metrics)

def get_best_window(metrics, windows):
    '''
    values: an a n-dimensional numpy array, 
            each column is a site in the alignment
            each row is some metric appropriately normalised
    '''

    #1. Mak an empty array:
    #   columns = number of things in windows
    #   rows = number of rows in metrics
    all_sses = np.empty((metrics.shape[0], len(windows) ))
    all_sses[:] = np.NAN # safety first

    # 2. Get SSE for each cell in array
    for i, window in enumerate(windows):
        # get SSE's for a given window
        all_sses[:,i] = get_sses(metrics, window)

    # 3. Sum each column of array -> 1D array
    all_sses = np.sum(all_sses, 0)

    # 4. get index of minimum value in 1D array
    best_index = np.argmin(all_sses)

    # 5. best window is windows[index]
    best_window = windows[i]

    return best_window

def get_sses(metrics, window, weights):
    '''
    input: metrics is an array where each row is a metric
            and each column is a site
           window gives slice instructions for the array
    output: an array with one col and the same number of 
            rows as metrics, where each entry is the SSE
    '''

    sses = np.apply_along_axis(get_sse, 1, metrics, window)

    # apply weights
    mean_sses = np.apply_along_axis(np.mean, 1, sses)
    scaling_factors = weights * means
    sses = sses / scaling_factors[:,None]


    return(sses)


def get_sse(metric, window):
    '''
    slice the 1D array metric, add up the SSES
    '''

    left  = sse(metric[ : window[0]])
    core  = sse(metric[window[0] : window[1]])
    right = sse(metric[window[1] : ])

    res = np.sum(left + core + right)

    return(res)


def sse(metric):
    '''
    input: list of values
    output: sum of squared errors
    '''
    sse = np.sum((metric - np.mean(metric))**2)

    return sse

def get_all_windows(aln, minimum_window_size=50):
    '''
    Input: aln: multiple sequence alignment
        minimum_window_size: smallest allowable window 
    Output: list of tuples [ (start : end) ]
    '''

    minimum_window_size = 50
    length = aln.get_alignment_length()

    keep_windows = []

    if length < 3*minimum_window_size:
        # some things can't be split
        return ([(0, length)])

    for window in combinations(range(length), 2):
        start = window[0]
        stop = window[1]

        if start < minimum_window_size:
            continue
        if (length - stop) < minimum_window_size:
            continue
        if (stop - start) < minimum_window_size:
            continue

        keep_windows.append(window)

    return (keep_windows)


def alignment_entropy(aln):
    '''
    input: biopython generic alignment
    output: entropy

    '''
    bp_freqs = bp_freqs_calc(aln)
    entropy = entropy_calc(bp_freqs)
    return (entropy)

def sitewise_entropies(aln):

    entropies = []
    for i in range(aln.get_alignment_length()):

        site_i = aln[:,i:i+1]
        ent_i = alignment_entropy(site_i)
        entropies.append(ent_i)

    # normalise and weight
    entropies = np.array(entropies)

    return (entropies)

def sitewise_gc(aln):

    gc = []
    for i in range(aln.get_alignment_length()):
        site_i = aln[:,i]
        gc_i = SeqUtils.GC(site_i)
        gc.append(gc_i)

    gc = np.array(gc)

    return (gc)

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


