from Bio.Nexus import Nexus
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq, UnknownSeq
import numpy as np
from tqdm import tqdm
from collections import defaultdict

dataset_path = '/Users/Tagliacollo/Desktop/scripts_UCE_test/matrix.nex'
#dataset_path = '/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/test.nex'

def charset_uce_aln(aln):

    dat = Nexus.Nexus()
    dat.read(aln)
    aln = AlignIO.read(open(aln), "nexus")

    uce_aln  = []

    for name in tqdm(dat.charsets):
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
    dd = defaultdict(list)

    for dicts in tuple_dicts:
        for key, value in dicts.items():
            dd[key].append(value)

    return(dd)

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




