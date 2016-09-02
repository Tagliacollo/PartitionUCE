import os
os.chdir("/Users/Tagliacollo/Desktop/ANU_Australia/scripts/")

######
from Bio.Nexus import Nexus
from Bio import AlignIO, SeqIO
from math import log

def entropy_calc(p):
    p = p[p!=0] # modify p to include only those elements that are not equal to 0
    return np.dot(-p,np.log2(p)) # the function returns the entropy result



###### read and split partition
def get_UCEs(matrix):
    aln = Nexus.Nexus()
    aln.read(matrix)
    uce_loci = aln.write_nexus_data_partitions(charpartition = aln.charsets)
    return uce_loci

x = get_UCEs('matrix.nex')
print (x)

###### loci in fasta format
def get_UCE_locus(partition):
    value = AlignIO.read(open(partition), "nexus")
    i = 0
    out =[]
    while i < value.get_alignment_length():
        loci = value[:,i]
        fas = '>' + 'locus_' + str(i) + '\n' + loci.replace('-', '') + '\n'
        out.append(fas)
        i = i + 1
    return out

y = get_UCE_locus('matrix_gen1.data')
print (y)

###### base frequencies
def get_bp_freqs(fasta):
    out = ['A', 'C', 'G', 'T']
    print (out)
    for line in SeqIO.parse(fasta, "fasta"):
        freqs = [line.seq.count(base)*1.0/len(line.seq) for base in ['A', 'C', 'G', 'T']]
        print(freqs)
        out.append(freqs)
    return out

z = get_bp_freqs("uce_loci.fasta")
print (z)

###### estimates of entropy per locus
def logent(x):
    if x <= 0:
        return 0
    else:
        return -x*log(x)

def UCE_entropy(fasta):
    out = []
    for line in SeqIO.parse(fasta, "fasta"):
        freqs = [line.seq.count(base)*1.0/len(line.seq) for base in ['A', 'C', 'G', 'T']]
        out.append(freqs)
        res = sum([logent(elem) for elem in freqs])
        out.append(res)
    return out

w = UCE_entropy("uce_loci.fasta")
print (w)


