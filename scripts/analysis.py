from utilities import *
from functions_metrics import *
from functions_site_position import *
from functions_full_multi import *
import os
from time import asctime


datasets = ['/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Moyle_2016.nex']
			#'/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Meiklejohn_2016.nex']
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/Crawford_2012.nex',
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/Faircloth_2013.nex',
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/McCormack_2013.nex',
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/Meiklejohn_2016.nex',
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/Moyle_2016.nex',
			#'/Users/Tagliacollo/Desktop/ANU_Australia/PartitionUCE/raw_data/Smith_2014.nex']

#			['/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Crawford_2012_bugtest.nex']
#		   ["/Users/roblanfear/Documents/github/PartitionUCE/raw_data/test.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Smith_2014.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Moyle_2016.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Crawford_2012.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Faircloth_2013.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/McCormack_2013.nex",
#			"/Users/roblanfear/Documents/github/PartitionUCE/raw_data/Meiklejohn_2016.nex"]

for dataset_path in datasets:

	print ("\n")
	print (dataset_path)

	os.chdir(output_paths(dataset_path))

	name = os.path.basename(dataset_path).rstrip(".nex")

	# WARNING: metrics must ALWAYS be ['entropy', 'gc', 'multi'] IN THAT ORDER!!!

	print(asctime())
	process_dataset_metrics(dataset_path, ['entropy', 'gc', 'multi'], minimum_window_size = 50, outfilename = '%s.csv' % (name))
	print(asctime())
	process_dataset_site_position(dataset_path, minimum_window_size = 50, outfilename = name)
	print(asctime())
	process_dataset_full_multi(dataset_path, minimum_window_size = 50, outfilename = name)
	print(asctime())



