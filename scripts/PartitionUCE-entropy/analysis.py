from utilities import *
from functions_metrics import *
import os
from time import asctime


datasets = ['/Users/Tagliacollo/Google-Drive/00-research/research-projects/in-progress/VAT-PartitionUCE/raw-data/test.nex']


for dataset_path in datasets:

	print ("\n")
	print (dataset_path)

	os.chdir(output_paths(dataset_path))

	name = os.path.basename(dataset_path).rstrip(".nex")

	print(asctime())
	process_dataset_metrics(dataset_path, ['entropy'], minimum_window_size = 50, outfilename = '%s.csv' % (name))
