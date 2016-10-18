import os
from pathlib2 import Path
from Bio import AlignIO
from itertools import combinations
from functions_metrics import bp_freqs_calc

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

def get_all_windows(aln, minimum_window_size=50):
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

        if start < minimum_window_size:
            continue
        if (length - stop) < minimum_window_size:
            continue
        if (stop - start) < minimum_window_size:
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

def p_finder_start_block(dataset_name, branchlengths = 'linked', models = 'all', model_selection = 'aicc'):
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


def p_finder_end_block(dataset_name, search = 'rcluster'):
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
                    'search = %s;\n\n' % (search) +

                    '#user schemes go here if search=user. See manual for how to define.#')

    return (end_block)