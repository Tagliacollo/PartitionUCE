from utilities import *
from functions_metrics import *
from functions_site_position import *
import os
from time import asctime


datasets = ['https://raw.githubusercontent.com/Tagliacollo/PartitionUCE/master/raw_data/Crawford_2012.nex',
			'https://raw.githubusercontent.com/Tagliacollo/PartitionUCE/master/raw_data/Harrington_2016.nex',
			'https://raw.githubusercontent.com/Tagliacollo/PartitionUCE/master/raw_data/McCormack_2013.nex',
			'https://raw.githubusercontent.com/Tagliacollo/PartitionUCE/master/raw_data/Meiklejohn_2016.nex',
			'https://raw.githubusercontent.com/Tagliacollo/PartitionUCE/master/raw_data/Moyle_2016.nex']

for dataset_path in datasets:

	print ("\n")
	print (dataset_path)

	os.chdir(output_paths(dataset_path))

	name = os.path.basename(dataset_path).rstrip(".nex")

	# WARNING: metrics must ALWAYS be ['entropy', 'gc', 'multi'] IN THAT ORDER!!!

	print(asctime())
	process_dataset_metrics(dataset_path, ['entropy', 'gc', 'multi'], minimum_window_size = 50, outfilename = '%s.csv' % (name))
	print(asctime())
	process_dataset_site_position(dataset_path, minimum_window_size = 5, outfilename = name)
	print(asctime())


