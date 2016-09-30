from functions import *
from Bio import AlignIO, SeqIO
import numpy as np
import sys, os

datasets = ['/Users/roblanfear/Documents/github/PartitionUCE/raw_data/test.nex']
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/Crawford_2012.nex',
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/Faircloth_2013.nex',
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/McCormack_2013.nex',
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/Meiklejohn_2016.nex',
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/Moyle_2016.nex',
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/Smith_2014.nex']

#		   ["/Users/roblanfear/Documents/github/PartitionUCE/raw_data/test.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Smith_2014.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Moyle_2016.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Crawford_2012.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Faircloth_2013.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/McCormack_2013.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Meiklejohn_2016.nex"]

#### Calculate statistics on all datasets #####
for dataset_path in datasets:

	print ("\n")
	print (dataset_path)

	# the order of weights has to be 'entropy', 'gc', 'multi' 
	weights = np.array([[1, 0, 0],
					    [0, 1, 0],
	                    [0, 0, 1],
	                    [1, 1, 0],
	                    [1, 0, 1],
	                    [0, 1, 1],
	                    [1, 1, 1]])

	weights = np.array([[1, 1, 1]])

	for weight in weights:
		os.chdir(output_paths(dataset_path, weight))
		process_dataset(dataset_path, ['entropy', 'gc', 'multi'], weights = weight, outfilename = 'VAT_test')


