from functions import *
from Bio import AlignIO, SeqIO
import numpy as np
from pathlib2 import Path
import sys


datasets = ["/Users/roblanfear/Documents/github/PartitionUCE/raw_data/test.nex"]


#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Smith_2014.nex"
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Moyle_2016.nex"
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Crawford_2012.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Faircloth_2013.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/McCormack_2013.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Meiklejohn_2016.nex"]

#### Calculate statistics on all datasets #####
for dataset_path in datasets:

	print ("\n")
	print (dataset_path)

	process_dataset(dataset_path, ['entropy'])
	process_dataset(dataset_path, ['gc'])
	process_dataset(dataset_path, ['entropy', 'gc'])


