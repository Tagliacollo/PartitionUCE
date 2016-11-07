from utilities import *
from Bio.Nexus import Nexus
from Bio import AlignIO, SeqIO, SeqUtils
import Bio
import numpy as np
import os, subprocess
from math import factorial
import itertools
from tqdm import tqdm
from functions_full_multi import *


def process_dataset_metrics(dataset_path, metrics, minimum_window_size, outfilename):
    '''
    Input: 
        dataset_path: path to a nexus alignment with UCE charsets
        metrics: a list of 'gc', 'entropy' or 'multi'
        outfilename: name for the csv file 
    Output: 
        csv files written to disk
    '''

    print("Sitewise metrics analysis")

    
    dataset_name = os.path.basename(dataset_path).rstrip(".nex")

    outfile = open(outfilename, 'w')
    outfile.write("name,uce_site,aln_site,window_start,window_stop,type,value,plot_mtx\n")
    outfile.close()

    # write the start blocks of the partitionfinder files
    for m in metrics:
        pfinder_config_file = open('%s_%s_partition_finder.cfg' % (dataset_name, m), 'w')
        pfinder_config_file.write(p_finder_start_block(dataset_name))
        pfinder_config_file.close()

    dat = Nexus.Nexus()
    dat.read(dataset_path)
    aln = AlignIO.read(open(dataset_path), "nexus")

    for name in tqdm(dat.charsets):

        sites = dat.charsets[name]
        start = min(sites)
        stop = max(sites) + 1
        # slice the alignment to get the UCE
        uce_aln = aln[:, start:stop]

        best_windows, metric_array = process_uce(uce_aln, metrics, minimum_window_size)

        for i, best_window in enumerate(best_windows):
            pfinder_config_file = open('%s_%s_partition_finder.cfg' % (dataset_name, metrics[i]), 'a')
            pfinder_config_file.write(blocks_pfinder_config(best_window, name, start, stop, uce_aln)) 
            break

        write_csvs(best_windows, metric_array, sites, name, outfilename)

    # write the end blocks of the partitionfinder files
    for m in metrics:
        pfinder_config_file = open('%s_%s_partition_finder.cfg' % (dataset_name, m), 'a')
        pfinder_config_file.write(p_finder_end_block(dataset_name))
        pfinder_config_file.close()


def write_csvs(best_windows, metric_array, aln_sites, name, outfilename):
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

    # make our long lists of things common to all metrics
    uce_sites = [uce_sites] * 3
    uce_sites_single_matrix = list(itertools.chain.from_iterable(uce_sites))
    outfile = open(outfilename, 'a')
    names = [name] * N * 3
    aln_sites = [np.array(aln_sites)] * 3 
    aln_sites_single_matrix = list(itertools.chain.from_iterable(aln_sites))

    # have to separate this in three lists
    window_start = []
    window_stop = []
    matrix_value = []
    
    for w in best_windows:
        window_start.append([w[0]]*N)
        window_stop.append([w[1]]*N)
        
        matrix_value.append(csv_col_to_plot_matrix(w, N))

    window_start_entropy = window_start[0]
    window_start_gc      = window_start[1]
    window_start_multi   = window_start[2]

    window_stop_entropy = window_stop[0]
    window_stop_gc      = window_stop[1]
    window_stop_multi   = window_stop[2]

    matrix_value_single_matrix = list(itertools.chain.from_iterable(matrix_value))

    window_starts = np.concatenate([window_start_entropy, window_start_gc, window_start_multi])
    window_stops = np.concatenate([window_stop_entropy, window_stop_gc, window_stop_multi])

    # build our types
    entropy_type = ['entropy'] * N
    gc_type = ['gc'] * N
    multi_type = ['multi'] * N

    metric_type = np.concatenate([entropy_type, gc_type, multi_type])

    metric_value = np.concatenate([metric_array[0], metric_array[1], metric_array[2]])

    all_info = [names, uce_sites_single_matrix, aln_sites_single_matrix, 
                window_starts, window_stops, metric_type, metric_value, 
                matrix_value_single_matrix]

    all_info = zip(*all_info)

    for i in all_info:
        line = [str(thing) for thing in i]
        line = ','.join(line)
        outfile.write(line)
        outfile.write("\n")
    outfile.close()


def process_uce(aln, metrics, minimum_window_size):
    '''
    Input:
        aln: biopython generic alignment
        metrics: a list with values of 'gc', 'entropy' or 'multi'
    Output:
        best_window: the best window for each metric
        metrics: a list with values of 'gc', 'entropy' or 'multi'
    '''
        
    windows = get_all_windows(aln, minimum_window_size)
    
    entropy = sitewise_entropies(aln)
    gc      = sitewise_gc(aln)
    multi   = sitewise_multi(aln)

    metrics = np.array([entropy, gc, multi])

    # sometimes we can't split a UCE, in which case there's one
    # window and it's the whole UCE
    


    if(len(windows)>1):
        best_window = get_best_windows(metrics, windows, aln.get_alignment_length())
    else:
        best_window = [windows[0], windows[0], windows[0]]
    
    return (best_window, metrics)

def get_best_windows(metrics, windows, aln_length):
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

    # 3. get index of minimum value for each metric
    entropy_sses = all_sses[0]
    gc_sses = all_sses[1]
    multi_sses = all_sses[2]

    entropy_mins = np.where(entropy_sses == entropy_sses.min())[0]  
    gc_mins = np.where(gc_sses == gc_sses.min())[0]
    multi_mins = np.where(multi_sses == multi_sses.min())[0]    


    # choose the windows with the minimum variance in length of
    #l-flank, core, r-flank
    entropy_wins = [windows[i] for i in entropy_mins]
    gc_wins = [windows[i] for i in gc_mins]
    multi_wins = [windows[i] for i in multi_mins]


    entropy = get_min_var_window(entropy_wins, aln_length)
    gc = get_min_var_window(gc_wins, aln_length)
    multi = get_min_var_window(multi_wins, aln_length)

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

def sitewise_multi(uce_aln):
    '''
    Input: 
        aln: biopython generic alignment
    Output: 
        1D numpy array with multinomial values for each site
    '''

    # preprocess stuff we need a lot
    uce_counts  = sitewise_base_counts(uce_aln)
    uce_sums = uce_counts.sum(axis = 0)

    uce_n_obs_factorial = []
    for s in uce_sums:
        f = factorial(s)
        uce_n_obs_factorial.append(f)


    uce_factorials = factorial_matrix(uce_counts)

    # sanity checks
    m1 = np.min(uce_factorials)
    m2 = np.min(uce_n_obs_factorial)
    m3 = min([m1, m2])
    if(m3<0):
        raise ValueError('One of your factorials is <0, this is bad, quitting)')

    sitewise_likelihoods = sitewise_multi_count(uce_counts, uce_factorials, uce_n_obs_factorial)

    log_likelihoods = np.log(sitewise_likelihoods)

    return(log_likelihoods)


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



