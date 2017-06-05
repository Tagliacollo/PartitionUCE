
for folder in $(find $PWD -type d); do
	python /data/programs/github/partitionfinder/PartitionFinder.py --raxml --quick $folder
done
