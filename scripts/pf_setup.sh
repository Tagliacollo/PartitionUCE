fi="Harrington_2016"
word="entropy"


ty="PF"$word
pf="*"$word"_partition_finder.cfg"
pd="/home/rlanfear/github/PartitionUCE/processed_data/"
rd="../raw_data/"
hm=$pd$fi"/"

unzip $pd$fi".zip" -d $pd
ad=$hm$ty
mkdir $ad
cp $hm$pf  $ad"/partition_finder.cfg"
cp $rd$fi".phy" $ad"/"$fi".phy"
python ~/github/partitionfinder/PartitionFinder.py $ad -q --force --raxml -p 10
