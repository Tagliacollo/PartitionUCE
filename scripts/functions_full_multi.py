from functions_metrics import sitewise_multi, get_all_windows
import numpy as np
from Bio import AlignIO
from Bio.Nexus import Nexus
from collections import defaultdict
from tqdm import tqdm
from utilities import *
import os

def process_dataset_full_multi(dataset_path, minimum_window_size, outfilename):

    dataset_name = os.path.basename(dataset_path).rstrip(".nex")

    dat = Nexus.Nexus()
    dat.read(dataset_path)
    aln = AlignIO.read(open(dataset_path), "nexus")


    pfinder_config_file = open('%s_full_multi_partition_finder.cfg' % (dataset_name), 'w')
    pfinder_config_file.write(p_finder_start_block(dataset_name))
    pfinder_config_file.close()

    lik_windows = defaultdict(list)
    for name in dat.charsets:
        print(name)
        sites = dat.charsets[name]
        start = min(sites)
        stop = max(sites) + 1
        # slice the alignment to get the UCE
        uce_aln = aln[:, start:stop]

        windows = get_all_windows(uce_aln, minimum_window_size)
    
        lik_windows = defaultdict(list)
        for window in tqdm(windows):
            lik_windows[window] = best_window_full_log_multi(uce_aln, window)

        best_window = max(lik_windows.items(), key=lambda a: a[1]) # should we take max or min? 
        pfinder_config_file = open('%s_full_multi_partition_finder.cfg' % (dataset_name), 'a')
        pfinder_config_file.write(export_charset_full_multi(name, best_window, start, stop, 
                                                            outfilename = pfinder_config_file))

    pfinder_config_file = open('%s_full_multi_partition_finder.cfg' % (dataset_name), 'a')
    pfinder_config_file.write(p_finder_end_block(dataset_name))
    pfinder_config_file.close()

    return (best_window)

def best_window_full_log_multi(aln, window):

   all_lik_multi = np.sum(np.log(sitewise_full_multi(aln, window)))

   return(all_lik_multi)
    
def sitewise_full_multi(aln, window):
    # aln[species :, aln_start : aln_end]
    uce_left  = sitewise_multi(aln[ :, : window[0]])
    uce_core  = sitewise_multi(aln[ :, window[0] : window[1]])
    uce_right = sitewise_multi(aln[ :, window[1] : ])

    return (np.concatenate([uce_left,uce_core,uce_right]))

def export_charset_full_multi(charset_name, best_window, start, stop, outfilename):
 
        # left UCE
        left_start  = start + 1
        left_end = left_start + best_window[0][0]
        left_UCE = '%s_left = %s-%s;\n' % (charset_name, left_start, left_end)
    
        # core UCE
        core_start = left_end + 1
        core_end = core_start + (best_window[0][1] - best_window[0][0] - 1)
        core_UCE = '%s_core = %s-%s;\n' % (charset_name, core_start, core_end)

        #right UCE
        right_start = core_end + 1
        right_end = stop
        right_UCE = '%s_right = %s-%s;\n' % (charset_name, right_start, right_end)

        # sometimes we couldn't split the window so it's all core
        if(best_window[0][1]-best_window[0][0] == stop-start):
            return(core_UCE)
        else:
            return(left_UCE + core_UCE + right_UCE)
    
