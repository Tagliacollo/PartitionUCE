from permutation import *
import os

n_permut = 25

datasets = ['/Volumes/VATagliacollo/GitHub/PartitionUCE/raw_data/Harrington_2016.nex',
			'/Volumes/VATagliacollo/GitHub/PartitionUCE/raw_data/Crawford_2012.nex',
			'/Volumes/VATagliacollo/GitHub/PartitionUCE/raw_data/McCormack_2013.nex',
			'/Volumes/VATagliacollo/GitHub/PartitionUCE/raw_data/Meiklejohn_2016.nex',
			'/Volumes/VATagliacollo/GitHub/PartitionUCE/raw_data/Moyle_2016.nex']

#### Generate random alns
for dataset_path in datasets:
	print ("\n")
	print (dataset_path)

	dataset_name   = os.path.basename(dataset_path).rstrip(".nex")
	repository_dir = Path(dataset_path).parents[1]
	processed_dir  = os.path.join(str(repository_dir), "processed_data/permutations/")
	output_path    = os.path.join(processed_dir, dataset_name)

	if not os.path.exists(output_path):
		os.makedirs(output_path)

	os.chdir(output_path)

	process_permutations(dataset_path, n_permut)

#### Generate random alns