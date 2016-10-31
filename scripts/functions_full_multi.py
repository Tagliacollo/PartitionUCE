from functions_metrics import *
import numpy as np
from Bio import AlignIO
from Bio.Nexus import Nexus
from collections import defaultdict
from tqdm import tqdm
from utilities import *
import os

def process_dataset_full_multi(dataset_path, minimum_window_size, outfilename):

    print("Full multinomial likelihood analysis")
    dataset_name = os.path.basename(dataset_path).rstrip(".nex")

    outfile = open(outfilename, 'w')
    outfile.write("name, uce_site, aln_site, window_start, window_stop, type, value\n")
    outfile.close()

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

        windows = get_all_windows(uce_aln, minimum_window_size)
 
        best_likelihood = np.inf * -1  # start at the bottom

        log_likelihoods = []
        for window in windows:
            lnL = multinomial_likelihood(uce_counts, uce_factorials, uce_n_obs_factorial, window)
            log_likelihoods.append(lnL)

        log_likelihoods = np.array(log_likelihoods)

        # indices of ML solutions
        ML_indices = np.where(log_likelihoods == log_likelihoods.max())[0]  

        ML_wins = [windows[i] for i in ML_indices]

        # if >1 alignment with equal likelihod, choose the one with minimum variance in lenghts
        best_window = get_min_var_window(ML_wins, aln.get_alignment_length())

        ## At this point, we need to go and re-calculate the likelihoods of
        # the best window, so we can output them to the csv file
        best_sitewise_likelihoods = sitewise_full_multi(uce_counts, uce_factorials, uce_n_obs_factorial, best_window)
        log_likelihoods = np.log(best_sitewise_likelihoods)

        #TODO: now save those log_likeilhoods to the csv file
        write_csvs_full(best_window, log_likelihoods, sites, name, outfilename)

        
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



def write_csvs_full(best_window, metric_full, aln_sites, name, outfilename):
    '''
    Input: 
        best_windows: the best window for full_multi
        metric: a list with values of log likelihoods
        aln_sites: site number in the alignment 
        name: UCE name
        outfilename: name for the csv file
    Output: 
        csv files written to disk
    '''

    N = len(aln_sites)
    middle = int(float(N) / 2.0)
    uce_sites = np.array(range(N)) - middle

    wd_start = [best_window[0]]*N
    wd_stop  = [best_window[1]]*N 

    name = [name]*N
    metric_type = ['full_multi']*N
    
    all_info = [name, uce_sites, aln_sites, wd_start, wd_stop,
                metric_type, metric_full]

    all_info = zip(*all_info)

    outfile = open(outfilename, 'a')
    
    for i in all_info:
        line = [str(thing) for thing in i]
        line = ','.join(line)
        outfile.write(line)
        outfile.write("\n")
    outfile.close()

