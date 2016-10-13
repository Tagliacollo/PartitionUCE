from Bio.Nexus import Nexus
from Bio import AlignIO, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, UnknownSeq
import numpy as np
from tqdm import tqdm
from collections import defaultdict, OrderedDict
import sys, os
from pathlib2 import Path

def process_dataset(aln, outfilename):
    read_aln = charset_uce_aln(aln)

    uce_dicts = []
    for uce in tqdm(read_aln):
        uce_dicts.append(dict_uce_sites(uce))

    conc_dicts = OrderedDict(conc_dicts_by_key(uce_dicts))

    merge_small_aln = conc_end_flanks(conc_dicts, min_aln_size = 100)

    return (merge_small_aln)

    
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

def charset_uce_aln(aln):

    dat = Nexus.Nexus()
    dat.read(aln)
    aln = AlignIO.read(open(aln), "nexus")

    uce_aln  = []

    for name in dat.charsets:
        sites = dat.charsets[name]
        start = min(sites)
        stop = max(sites) + 1
        # slice the alignment to get the UCE
        uce_aln.append(aln[:, start : stop])

    return(uce_aln)

def dict_uce_sites(aln):
    N        = aln.get_alignment_length()
    middle   = int(float(N) / 2.0)
    site_num = np.array(range(N)) - middle

    sites = {}
    for i, site in enumerate(range(aln.get_alignment_length())):
        sites[site_num[i]] = (aln[: , site : site+1])

    return(sites)

def conc_dicts_by_key(tuple_dicts):
# inspired by http://stackoverflow.com/questions/5946236/how-to-merge-multiple-dicts-with-same-key
    uce_dicts = defaultdict(list)

    for dicts in tuple_dicts:
        for key, value in dicts.items():
            uce_dicts[key].append(value)

    return(uce_dicts)

def conc_end_flanks(conc_dicts, min_aln_size):
    for i in range(min(conc_dicts),0):
        key = i
        if len(conc_dicts[key]) < min_aln_size:
            conc_dicts[key + 1] = conc_dicts[key] + conc_dicts[key + 1]
            del conc_dicts[key]

    for i in range(0, max(conc_dicts)):
        key = i
        if len(conc_dicts[key]) < min_aln_size:
            conc_dicts[key + 1] = conc_dicts[key] + conc_dicts[key + 1]
            del conc_dicts[key]

    if len(conc_dicts[max(conc_dicts)]) < min_aln_size:

        index_max = list(conc_dicts.keys()).index(max(conc_dicts))
        key_max   = list(conc_dicts.keys())[index_max]
    
        index_bef_max = index_max - 1
        key_bef_max   = list(conc_dicts.keys())[index_bef_max]

        conc_dicts[key_bef_max] = conc_dicts[key_bef_max] + conc_dicts[key_max]
        del conc_dicts[key_max]

    if len(conc_dicts[min(conc_dicts)]) < min_aln_size:

        index_min = list(conc_dicts.keys()).index(min(conc_dicts))
        key_min  = list(conc_dicts.keys())[index_min]
    
        index_bef_min = index_min + 1
        key_bef_min   = list(conc_dicts.keys())[index_bef_min]

        conc_dicts[key_bef_min] = conc_dicts[key_bef_min] + conc_dicts[key_min]
        del conc_dicts[key_min]

    return(conc_dicts)

def conc_site_in_aln(list_aln):
    # from https://gist.github.com/kgori/f0532cff6500e839cb29
    #       function name changed from concatenate to conc_site_in_aln
    #       argument name changed from alignment to list_aln
    
    """
    Concatenates a list of Bio.Align.MultipleSeqAlignment objects.
    If any sequences are missing the are padded with unknown data
    (Bio.Seq.UnknownSeq).
    Returns a single Bio.Align.MultipleSeqAlignment.
    Limitations: any annotations in the sub-alignments are lost in
    the concatenated alignment.
    """
    
    # Get the full set of labels (i.e. sequence ids) for all the alignments
    all_labels = set(seq.id for aln in list_aln for seq in aln)
    
    # Make a dictionary to store info as we go along
    # (defaultdict is convenient -- asking for a missing key gives back an empty list)
    tmp = defaultdict(list)
    
    # Assume all alignments have same alphabet
    alphabet = list_aln[0]._alphabet
    
    for aln in list_aln:
        length = aln.get_alignment_length()
        
        # check if any labels are missing in the current alignment
        these_labels = set(rec.id for rec in aln)
        missing = all_labels - these_labels
        
        # if any are missing, create unknown data of the right length,
        # stuff the string representation into the tmp dict
        for label in missing:
            new_seq = UnknownSeq(length, alphabet=alphabet)
            tmp[label].append(str(new_seq))
        
        # else stuff the string representation into the tmp dict
        for rec in aln:
            tmp[rec.id].append(str(rec.seq))
            
    # Stitch all the substrings together using join (most efficient way),
    # and build the Biopython data structures Seq, SeqRecord and MultipleSeqAlignment

    return MultipleSeqAlignment(SeqRecord(Seq(''.join(v), alphabet=alphabet), id=k) 
               for (k,v) in tmp.items())



