from Bio.Nexus import Nexus
from Bio import AlignIO, SeqIO, SeqUtils
import Bio
import numpy as np
import os, subprocess
from glob import glob
from itertools import combinations
from itertools import islice
from pathlib2 import Path
from math import factorial

def output_paths(dataset_path):
    '''
    Input: 
        dataset_path: path to a nexus alignment with UCE charsets 
    Ouput: 
        folder path with the name of the nexus UCE dataset  
    '''
    
    dataset_name = os.path.basename(dataset_path).rstrip(".nex")

    repository_dir      = Path(dataset_path).parents[1]
    processed_data_dir  = os.path.join(str(repository_dir), "processed_data")

    output_path = os.path.join(processed_data_dir, dataset_name)
    
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    return (output_path)


def process_dataset(dataset_path, metrics, outfilename):
    '''
    Input: 
        dataset_path: path to a nexus alignment with UCE charsets
        metrics: a list of 'gc', 'entropy' or 'multi'
        outfilename: name for the csv file 
    Output: 
        csv files written to disk
    '''
    outfile = open(outfilename, 'w')
    outfile.write("name,uce_site,aln_site,window_start,window_stop,%s\n" %(','.join(metrics)))
    outfile.close()

    dat = Nexus.Nexus()
    dat.read(dataset_path)
    aln = AlignIO.read(open(dataset_path), "nexus")

    for name in dat.charsets:
        print (name)
        sites = dat.charsets[name]
        start = min(sites)
        stop = max(sites) + 1
        # slice the alignment to get the UCE
        uce_aln = aln[:, start:stop]

        best_windows, metric_array = process_uce(uce_aln, metrics)

        #write_csvs(best_windows, metric_array, sites, name, outfilename)
        print(best_windows)

def write_csvs(best_windows, metrics, aln_sites, name, outfilename):
    '''
    Input: 
        best_windows: the best window for each metric
        metrics: a list with values of 'gc', 'entropy' and 'multi'
        aln_sites: site number in the alignment 
        name: UCE name
        outfilename: name for the csv file
    Output: 
        csv files written to disk
    '''

    N = len(aln_sites)
    middle = int(float(N) / 2.0)
    uce_sites = np.array(range(N)) - middle

    outfile = open(outfilename, 'a')

    # TODO: account for different windows from different metrics

    names = [name]*N
    
    window_start = []
    window_stop = []
    for w in best_windows:
        window_start.append([w[0]]*N)
        window_stop.append([w[1]]*N)

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
    '''
    Input:
        aln: biopython generic alignment
        metrics: a list with values of 'gc', 'entropy' or 'multi'
    Output:
        best_window: the best window for each metric
        metrics: a list with values of 'gc', 'entropy' or 'multi'
    '''
        
    windows = get_all_windows(aln)
    
    entropy = sitewise_entropies(aln)
    gc      = sitewise_gc(aln)
    multi   = sitewise_multi(aln)

    metrics = np.array([entropy, gc, multi])

    best_window = get_best_windows(metrics, windows)

    return (best_window, metrics)

def get_best_windows(metrics, windows):
    '''
    Input: 
        an a n-dimensional numpy array, 
        each column is a site in the alignment
        each row is some metric appropriately normalised
    Output: 
        the best window for each metric
    '''

    #1. Mak an empty array:
    #   columns = number of things in windows
    #   rows = number of rows in metrics
    all_sses = np.empty((metrics.shape[0], len(windows) ))
    all_sses[:] = np.NAN # safety first

    #print(metrics, windows)

    # 2. Get SSE for each cell in array
    for i, window in enumerate(windows):
        # get SSE's for a given window
        all_sses[:,i] = get_sses(metrics, window)

    # 3. get index of minimum value in 1D array FOR each metric
    best_indices = np.apply_along_axis(np.argmin, 1, all_sses)
    
    # 4. best window is windows[index]
    # get the best window for each metric
    # i.e. best_indices gives the index for each row

    entropy = windows[best_indices[0]]
    gc      = windows[best_indices[1]]
    multi   = windows[best_indices[2]]

    best_windows = [entropy, gc, multi]

    return (best_windows)

def get_sses(metrics, window):
    '''
    Input: 
        metrics is an array where each row is a metric
        and each column is a site window gives slice 
        instructions for the array
    Output: 
        an array with one col and the same number of 
        rows as metrics, where each entry is the SSE
    '''

    sses = np.apply_along_axis(get_sse, 1, metrics, window)

    return (sses)


def get_sse(metric, window):
    '''
    slice the 1D array metric, add up the SSES
    '''

    left  = sse(metric[ : window[0]])
    core  = sse(metric[window[0] : window[1]])
    right = sse(metric[window[1] : ])

    res = np.sum(left + core + right)

    return (res)


def sse(metric):
    '''
    input: 
        list with values of 'gc', 'entropy' and 'multi'
    output: 
        sum of squared errors
    '''
    sse = np.sum((metric - np.mean(metric))**2)

    return (sse)

def get_all_windows(aln, minimum_window_size=50):
    '''
    Input: 
        aln: multiple sequence alignment
        minimum_window_size: smallest allowable window 
    Output: 
        list of all possible tuples [ (start : end) ]
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
    Input: 
        aln: biopython generic alignment
    Output: 
        array with values of entropies
    '''

    bp_freqs = bp_freqs_calc(aln)
    entropy = entropy_calc(bp_freqs)
    
    return (entropy)

def sitewise_entropies(aln):
    '''
    Input: 
        aln: biopython generic alignment
    Output: 
        array with values of entropies per site
    '''

    entropies = []
    for i in range(aln.get_alignment_length()):

        site_i = aln[:,i:i+1]
        ent_i = alignment_entropy(site_i)
        entropies.append(ent_i)

    entropies = np.array(entropies)

    return (entropies)

def sitewise_gc(aln):
    '''
    Input: 
        aln: biopython generic alignment
    Output: 
        array with values of gc per site
    '''

    gc = []
    for i in range(aln.get_alignment_length()):
        site_i = aln[:,i]
        gc_i = SeqUtils.GC(site_i)
        gc.append(gc_i)

    gc = np.array(gc)

    return (gc)

def sitewise_multi(aln):
    '''
    Input: 
        aln: biopython generic alignment
    Output: 
        1D numpy array with multinomial values for each site
    '''

    number_ssp = len(aln)
    bp_freqs = bp_freqs_calc(aln)

    prop_A = bp_freqs[0]
    prop_C = bp_freqs[1]
    prop_G = bp_freqs[2]
    prop_T = bp_freqs[3]

    multinomial_results = []
    for i in range(aln.get_alignment_length()):
        site = aln[:,i]
    
        count_A = site.count('A')
        count_C = site.count('C')
        count_G = site.count('G')
        count_T = site.count('T')

        # Function to calculate multinomial - OBS: numpy has no function for factorial calculations
        # for factorial calculations: from math import factorial
        N = factorial(number_ssp)
        K = factorial(count_A) * factorial(count_C) * factorial(count_G) * factorial(count_T)
        J = (prop_A**count_A * prop_C**count_C * prop_G**count_G * prop_T**count_T)

        multinomial_cal = (N/K)*J
        
        multinomial_results.append(multinomial_cal)

    return (np.array(multinomial_results))  

def bp_freqs_calc(aln):
    '''
    Input: 
        aln: biopython generic alignment
    Output: 
        1-D array of base frequencies
    '''
    one_str = ""
    for seq in aln:
        one_str += seq

    seq = one_str.upper()

    A = seq.seq.count('A') 
    C = seq.seq.count('C')
    G = seq.seq.count('G')
    T = seq.seq.count('T')

    sum_count = A + C + G + T

    bp_freqs = np.array([ A, C, G, T])/float(sum_count)
    
    return (bp_freqs)


def entropy_calc(p):
    '''
    Input: 
       p: 1D array of base frequencies
    Output: 
        estimates of entropies 
    
    copied from: http://nbviewer.ipython.org/url/atwallab.cshl.edu/teaching/QBbootcamp3.ipynb
    '''
    p = p[p!=0] # modify p to include only those elements that are not equal to 0
    
    return np.dot(-p,np.log2(p)) # the function returns the entropy result


