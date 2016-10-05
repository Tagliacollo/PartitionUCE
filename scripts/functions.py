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

def blocks_pfinder_config(best_window, name, start, stop, outfinename):

    # left UCE
    left_start  = start
    left_end = left_start + best_window[0]
    left_UCE = '%s_left = %s-%s;\n' % (name, left_start, left_end)

    # core UCE
    core_start = left_end + 1
    core_end = core_start + (best_window[1] - best_window[0] - 1)
    core_UCE = '%s_core = %s-%s;\n' % (name, core_start, core_end)

    #right UCE
    right_start = core_end + 1
    right_end = stop
    right_UCE = '%s_right = %s-%s;\n' % (name, right_start, right_end)

    return (left_UCE + core_UCE + right_UCE)


def p_finder_start_block(dataset_name, branchlengths = 'linked', models = 'all', model_selection = 'aicc'):
    begin_block = str('## ALIGNMENT FILE ##\n' + 
                      'alignment = %s.phy;\n\n' % (dataset_name) +  


                      '## BRANCHLENGTHS: linked | unlinked ##\n' +
                      'branchlengths = %s;\n\n' % (branchlengths) +

                       '## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai <list> ##\n' +
                       'models = %s;\n\n' % (models) + 

                       '# MODEL SELECCTION: AIC | AICc | BIC #' +
                       'model_selection = %s;\n\n' % (model_selection) +

                       '## DATA BLOCKS: see manual for how to define ##\n' +
                       '[data_blocks]\n')

    return (begin_block)

def p_finder_end_block(dataset_name, search = 'rcluster'):
    end_block = str('\n' +
                    '## SCHEMES, search: all | user | greedy | rcluster | hcluster | kmeans ##\n' +
                    '[schemes]\n' +
                    'search = %s;\n\n' % (search) +

                    '#user schemes go here if search=user. See manual for how to define.#')

    return (end_block)

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
    dataset_name = os.path.basename(dataset_path).rstrip(".nex")

    outfile = open(outfilename, 'w')
    outfile.write("name,uce_site,aln_site,entropy_wd_start,entropy_wd_stop,gc_wd_start,gc_wd_stop,multi_wd_start,multi_wd_stop,%s\n" %(','.join(metrics)))
    outfile.close()

    # write the start blocks of the partitionfinder files
    for m in metrics:
        pfinder_config_file = open('%s_%s_partition_finder.cfg' % (dataset_name, m), 'w')
        pfinder_config_file.write(p_finder_start_block(dataset_name))
        pfinder_config_file.close()

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

        for i, best_window in enumerate(best_windows):
            pfinder_config_file = open('%s_%s_partition_finder.cfg' % (dataset_name, metrics[i]), 'a')
            pfinder_config_file.write(blocks_pfinder_config(best_window, name, start, stop, 
                                           outfinename = pfinder_config_file))


        write_csvs(best_windows, metric_array, sites, name, outfilename)

    # write the end blocks of the partitionfinder files
    for m in metrics:
        pfinder_config_file = open('%s_%s_partition_finder.cfg' % (dataset_name, m), 'a')
        pfinder_config_file.write(p_finder_end_block(dataset_name))
        pfinder_config_file.close()


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

    names = [name]*N
    

    # have to separate this in three lists
    window_start = []
    window_stop = []
    for w in best_windows:
        window_start.append([w[0]]*N)
        window_stop.append([w[1]]*N)

    window_start_entropy = window_start[0]
    window_start_gc      = window_start[1]
    window_start_multi   = window_start[2]

    window_stop_entropy = window_stop[0]
    window_stop_gc      = window_stop[1]
    window_stop_multi   = window_stop[2]


    metrics = metrics.tolist()

    all_info = [names, uce_sites, aln_sites, window_start_entropy, window_stop_entropy, 
                window_start_gc, window_stop_gc, window_start_multi, window_stop_multi]

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

def full_sitewise_likelihood_multi(aln, window):
    # aln[species :, aln_start : aln_end]
    uce_left  = sitewise_multi(aln[ :, : window[0][0]])
    print(uce_left)
    uce_core  = sitewise_multi(aln[ :, window[0][0] : window[0][1]])
    uce_right = sitewise_multi(aln[ :, window[0][1] : ])
    print(uce_right)

    return (np.concatenate([uce_left,uce_core,uce_right]))


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


