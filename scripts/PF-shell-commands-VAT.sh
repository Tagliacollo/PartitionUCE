#fi="Crawford_2012"
#fi="Harrington_2016"
#fi="McCormack_2013"
#fi="Meiklejohn_2016"
fi="Moyle_2016"

#word="subset_uce"
word="subset_aln"

ty="PF-mp-"$word
pf="*"$word"_partition_finder.cfg"
pd="/home/vtagliacollo/GitHub/PartitionUCE/"
rd="/home/vtagliacollo/GitHub/PartitionUCE/raw_data/"
hm=$pd$fi"/"

ad=$hm$ty
mkdir $ad
cp $hm$pf  $ad"/partition_finder.cfg"
cp $rd$fi".phy" $ad"/"$fi".phy"

python /disks/dacelo/data/programs/github/partitionfinder/PartitionFinder.py $ad --raxml -n -p 10
