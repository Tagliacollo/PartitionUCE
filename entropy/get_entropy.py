# working directory 
import os
os.chdir("/Users/Tagliacollo/Desktop/ANU_Australia/scripts/")

### packages
from Bio.Nexus import Nexus
from Bio import AlignIO, SeqIO

### read alignmet
# matrix.nex is a hypothetical dataset. 
aln = Nexus.Nexus()
aln.read("matrix.nex") 

### isolate partitions
uce = aln.write_nexus_data_partitions(charpartition = aln.charsets)

### Get loci in fasta format
# This part of the code works, but the output is been overwrited (i.e. I get only one output instead of two).  
# There is certainly a easier way of doing this loop. Do you have any suggestion to help me out?
# One more question: should we save outputs in fasta format to run later?
for key in uce:
    value = AlignIO.read(open(uce[key]), "nexus")
    print (uce[key])
    output = open("uce_loci.fasta", "w")
    i = 0
    while i < value.get_alignment_length():
        loci = value[:,i]
        print ('>' + 'locus_' + str(i) + '\n' + loci.replace('-', '') + '\n')
        output.write ('>' + 'locus_' + str(i) + '\n' + loci.replace('-', '') + '\n')
        i = i + 1

### calculate entropy from loci
# I'm still trying to understand how this piece of code works. I copy it from http://stackoverflow.com/questions/37909873/how-to-calculate-the-entropy-of-a-dna-sequence-in-a-fasta-file
# It isn't working because 'SeqRecord' object (from SeqIO) has no attribute 'count'. What I've understood so far is that "*.count" is a built in function in python. Why doesn't it work in my code?
from math import log  

def logent(x):  
    if x <= 0:     
        return 0  
    else:  
        return - x*log(x)  

def entropy(lis):   
    return sum([logent(elem) for elem in lis])

for seq in SeqIO.parse("uce_loci.fasta", "fasta"):
    lisfreq1 = [seq.count(base)*1.0/len(seq) for base in ["A", "C","G","T"]]
    x = entropy(lisfreq1)
    print (x)



