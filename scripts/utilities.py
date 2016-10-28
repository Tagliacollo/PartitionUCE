import os
from pathlib2 import Path
from Bio import AlignIO, SeqIO, SeqUtils
from itertools import combinations
import numpy as np
from math import factorial
from Bio.Nexus import Nexus

def check_taxa(matrices):
    '''Checks that nexus instances in a list [(name, instance)...] have
        the same taxa, provides useful error if not and returns None if
        everything matches
        From: http://biopython.org/wiki/Concatenate_nexus
    '''
    first_taxa = matrices[0][1].taxlabels
    for name, matrix in matrices[1:]:
        first_only = [t for t in first_taxa if t not in matrix.taxlabels]
        new_only = [t for t in matrix.taxlabels if t not in first_taxa]

        if first_only:
            missing = ', '.join(first_only)
            msg = '%s taxa %s not in martix %s' % (matrices[0][0], missing, name)
            raise Nexus.NexusError(msg)

        elif new_only:
            missing = ', '.join(new_only)
            msg = '%s taxa %s not in all matrices'  % (name, missing)
            raise Nexus.NexusError(msg)

    return None # will only get here if it hasn't thrown an exception


def concat(mypath, same_taxa):
    ''' Combine multiple nexus data matrices in one partitioned file.
        By default this will only work if the same taxa are present in each file
        use  same_taxa=False if you are not concerned by this
        From: http://biopython.org/wiki/Concatenate_nexus
        small change: added onlyfiles block to remove hidden files
    '''

    onlyfiles = []
    for item in os.listdir(mypath):
        if not item.startswith('.') and os.path.isfile(os.path.join(mypath, item)):
            onlyfiles.append(item)
    
    nexi = []
    for nex in onlyfiles:
        nex_open = open(nex, 'r')
        nex_save = Nexus.Nexus(nex_open)
        nexi.append((nex, nex_save))
         
    if same_taxa:
        if not check_taxa(nexi):
            return Nexus.combine(nexi)
    else:
        return Nexus.combine(nexi)


def output_conc_nex(mypath, outfilename, same_taxa=False):

    os.chdir(mypath)
    combined = concat(mypath, same_taxa)
    combined.write_nexus_data(filename=open('%s.nex' % (outfilename), 'w'))

    return None


def blocks_pfinder_config(best_window, name, start, stop, uce_aln):

    # sometimes we couldn't split the window so it's all together
    if(best_window[1]-best_window[0] == stop-start):
        whole_UCE = '%s_all = %s-%s;\n' % (name, start+1, stop)
        return (whole_UCE)

    else:
        # left UCE
        left_start = start + 1
        left_end = left_start + best_window[0]

        # core UCE
        core_start = left_end + 1
        core_end = start + best_window[1]

        #right UCE
        right_start = core_end + 1
        right_end = stop

    # do not output any undetermined blocks - if this happens, just output the whole UCE
    if(any_undetermined_blocks(best_window, uce_aln)==True):
        whole_UCE = '%s_all = %s-%s;\n' % (name, start+1, stop)
        return (whole_UCE)
    else:
        core_UCE = '%s_core = %s-%s;\n' % (name, core_start, core_end)
        left_UCE = '%s_left = %s-%s;\n' % (name, left_start, left_end)
        right_UCE = '%s_right = %s-%s;\n' % (name, right_start, right_end)

        return (left_UCE + core_UCE + right_UCE)


def any_undetermined_blocks(best_window, uce_aln):
    # Return TRUE if there are any blocks with only undeteremined characters
    # Defined as anything other than ACGT

    left_aln = uce_aln[:, 0 : best_window[0]]
    core_aln = uce_aln[:, best_window[0] : best_window[1]]
    right_aln = uce_aln[:, best_window[1] : uce_aln.get_alignment_length()]

    l_freq = bp_freqs_calc(left_aln)
    c_freq = bp_freqs_calc(core_aln)
    r_freq = bp_freqs_calc(right_aln)

    if(np.isnan(l_freq.max()) or np.isnan(c_freq.max()) or np.isnan(r_freq.max())):
        return(True)
    else:
        return(False)

def get_all_windows(aln, minimum_window_size):
    '''
    Input: 
        aln: multiple sequence alignment
        minimum_window_size: smallest allowable window 
    Output: 
        list of all possible tuples [ (start : end) ]
    '''

    length = aln.get_alignment_length()

    keep_windows = []

    if length < 3*minimum_window_size:
        # some things can't be split
        return ([(0, length)])

    for window in combinations(range(length), 2):
        start = window[0]
        stop = window[1]

        if start <= minimum_window_size:
            continue
        if (length - stop) <= minimum_window_size:
            continue
        if (stop - start) <= minimum_window_size:
            continue

        keep_windows.append(window)

    return (keep_windows)


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

def p_finder_start_block(dataset_name, branchlengths = 'linked', models = 'GTR+G', model_selection = 'aicc'):
    begin_block = str('## ALIGNMENT FILE ##\n' + 
                      'alignment = %s.phy;\n\n' % (dataset_name) +  


                      '## BRANCHLENGTHS: linked | unlinked ##\n' +
                      'branchlengths = %s;\n\n' % (branchlengths) +

                       '## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai <list> ##\n' +
                       'models = %s;\n\n' % (models) + 

                       '# MODEL SELECCTION: AIC | AICc | BIC #\n' +
                       'model_selection = %s;\n\n' % (model_selection) +

                       '## DATA BLOCKS: see manual for how to define ##\n' +
                       '[data_blocks]\n')

    return (begin_block)


def p_finder_end_block(dataset_name, search = 'rclusterf'):
    '''
    args: 
        dataset_name: name of the dataset
        search: pFinder input arguments
    return:
        str with information about the end block pFinder config block
    '''
    
    end_block = str('\n' +
                    '## SCHEMES, search: all | user | greedy | rcluster | hcluster | kmeans ##\n' +
                    '[schemes]\n' +
                    'search = %s;\n\n' % (search)
                    )

    return (end_block)

def sitewise_base_counts(aln):
    '''
    Input: 
        aln: biopython generic alignment
    Output: 
        4xN array of base counts (A,C,G,T) by site
    '''
    n_sites = aln.get_alignment_length()

    base_counts = np.zeros((4, n_sites))

    for i in range(n_sites):

        site_i = aln[:,i:i+1]
        counts_i = count_bases(site_i)
        base_counts[:,i] = counts_i

    return(base_counts)

def factorial_matrix(counts):
    '''
    input an array of integers
    output an array of the colum-wise products of the factorials of those integers
    '''
    # NB: I used to do this like this, which seemed more sensible and almost certainly quicker,
    #fv = np.vectorize(factorial)
    #factorials = fv(counts)
    #f_product = factorials.prod(axis = 0)
    # however, I kept getting overflow errors that I couldn't figure out, so I gave up.

    f_product = []

    for c in counts.T:
        cf = []
        for i in c:
            cf.append(factorial(i))
        f_product.append(np.product(cf))


    return(f_product)


def count_bases(aln):

    one_str = ""
    for seq in aln:
        one_str += seq

    seq = one_str.upper()

    A = seq.seq.count('A') 
    C = seq.seq.count('C')
    G = seq.seq.count('G')
    T = seq.seq.count('T')

    return(np.array([A, C, G, T]))

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
