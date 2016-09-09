from functions import *
from Bio import AlignIO, SeqIO
import numpy as np
from pathlib2 import Path
import sys


datasets = ["/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Moyle_2016.nex"]

#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Crawford_2012.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Faircloth_2013.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/McCormack_2013.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Meiklejohn_2016.nex"]

tigger_path = "/Users/roblanfear/Documents/github/PartitionUCE/programs/tigger"

for dataset_path in datasets:

	print ("\n")
	print (dataset_path)

	# 1. Split alignment into dict of charsets
	uce_alns = split_charsets_to_list(dataset_path)
	dataset_name = os.path.basename(dataset_path).rstrip(".nex")

	# 2. Set up empty dicts for results
	uce_entropies   = {}
	uce_gc 			= {}
	uce_tiger 		= {}
	for i, key in enumerate(uce_alns):

		# progress
		done = (float(i)*100.0/float(len(uce_alns.keys())))
		sys.stdout.write("\r%.2f%% done (now analysing %s)\t\t\t\t" %(done, key))
		sys.stdout.flush()

		uce_entropies[key] 	= sitewise_entropies(uce_alns[key])	
		uce_gc[key] 		= sitewise_gc(uce_alns[key])	
		uce_tiger[key] 		= sitewise_TIGER(uce_alns[key], tigger_path)	

	# write csv files
	repo_directory = Path(dataset_path).parents[1]
	processed_data_dir = os.path.join(str(repo_directory), "processed_data")
	output_file_base = os.path.join(processed_data_dir, dataset_name)

	write_csv(uce_entropies, "%s_entropy.csv" %(output_file_base), 'value')
	write_csv(uce_gc, "%s_gc.csv" %(output_file_base), 'value')
	write_csv(uce_tiger, "%s_tiger.csv" %(output_file_base), 'value')

