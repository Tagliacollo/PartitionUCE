from functions_metrics import *
import numpy as np
from Bio import AlignIO
from Bio.Nexus import Nexus
from collections import defaultdict
from tqdm import tqdm
from utilities import *
import os
fv = np.vectorize(factorial)

def process_dataset_full_multi(dataset_path, minimum_window_size, outfilename):

    print("Full multinomial likelihood analysis")
    dataset_name = os.path.basename(dataset_path).rstrip(".nex")

    dat = Nexus.Nexus()
    dat.read(dataset_path)
    aln = AlignIO.read(open(dataset_path), "nexus")


    pfinder_config_file = open('%s_full_multi_partition_finder.cfg' % (dataset_name), 'w')
    pfinder_config_file.write(p_finder_start_block(dataset_name))
    pfinder_config_file.close()

    lik_windows = defaultdict(list)
    for name in tqdm(dat.charsets):
        sites = dat.charsets[name]
        start = min(sites)
        stop = max(sites) + 1
        # slice the alignment to get the UCE
        uce_aln = aln[:, start:stop]

        # preprocess stuff we need a lot
        uce_counts  = sitewise_base_counts(uce_aln)
        uce_n_obs_factorial = fv(uce_counts.sum(axis = 0))
        uce_factorials = factorial_matrix(uce_counts)

        windows = get_all_windows(uce_aln, minimum_window_size)
 
        best_likelihood = np.inf * -1  # start at the bottom

        for window in windows:
            log_likelihood = multinomial_likelihood(uce_counts, uce_factorials, uce_n_obs_factorial, window)

            if(log_likelihood > best_likelihood):
                best_window = window
                best_likelihood = log_likelihood

        pfinder_config_file = open('%s_full_multi_partition_finder.cfg' % (dataset_name), 'a')
        pfinder_config_file.write(blocks_pfinder_config(best_window, name, start, stop, uce_aln))
    pfinder_config_file = open('%s_full_multi_partition_finder.cfg' % (dataset_name), 'a')
    pfinder_config_file.write(p_finder_end_block(dataset_name))
    pfinder_config_file.close()

    return (best_window)

def multinomial_likelihood(counts, factorials, Ns, window):

    sitewise_likelihoods = sitewise_full_multi(counts, factorials, Ns, window)
    log_likelihoods = np.log(sitewise_likelihoods)
    log_likelihood = np.sum(log_likelihoods)

    return(log_likelihood)
    
def sitewise_full_multi(counts, factorials, Ns, window):
    # aln[species :, aln_start : aln_end]
    uce_left  = sitewise_multi_count(counts[ :, : window[0]], factorials[ : window[0]], Ns[ : window[0]])
    uce_core  = sitewise_multi_count(counts[ :, window[0] : window[1]], factorials[window[0] : window[1]], Ns[window[0] : window[1]])
    uce_right = sitewise_multi_count(counts[ :, window[1] : ], factorials[window[1] : ], Ns[window[1] : ])

    return (np.concatenate([uce_left,uce_core,uce_right]))


def sitewise_multi_count(counts, factorials, Ns):
    '''
    Input: 
        counts: counts of A,C,G,T, in 4xN matrix, where N is number of sites
        factorials: products of the factorials of the counts (i.e. 1xN array)
        Ns: factorials of the sums of the counts (i.e. 1xN array)
    Output: 
        1D numpy array with multinomial values for each site
    '''

    n_sites = counts.shape[1]

    background_base_counts = counts.sum(axis = 1)
    background_base_freqs = background_base_counts/sum(background_base_counts)

    multinomial_likelihoods = np.zeros(counts.shape[1])

    for i in range(n_sites):

        K = factorials[i]
        N = Ns[i]
        counts_i = counts[:,i]
        J = np.product(background_base_freqs**counts_i)
        
        L = (N/K)*J
        
        multinomial_likelihoods[i] = L

    return (multinomial_likelihoods)  


