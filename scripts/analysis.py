from functions import *
from Bio import AlignIO, SeqIO
import numpy as np
from pathlib2 import Path

datasets = ["/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Crawford_2012.nex"]

for dataset_path in datasets:

	# 1. Split alignment into dict of charsets
	uce_alns = split_charsets_to_list(dataset_path)
	dataset_name = os.path.basename(dataset_path).rstrip(".nex")

	# 2. Set up empty dicts for results
	uce_entropies = {}
	uce_gc = {}

	for key in uce_alns:

		uce_entropies[key] 	= sitewise_entropies(uce_alns[key])	
		uce_gc[key] 		= sitewise_gc(uce_alns[key])	

	# write csv files
	repo_directory = Path(dataset_path).parents[1]
	processed_data_dir = os.path.join(str(repo_directory), "processed_data")
	output_file_base = os.path.join(processed_data_dir, dataset_name)
	write_csv(uce_entropies, "%s_entropy.csv" %(output_file_base), 'entropy')
	write_csv(uce_gc, "%s_gc.csv" %(output_file_base), 'gc')

