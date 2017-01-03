##### PFinder
#fi="Crawford_2012"
#fi="Harrington_2016"
#fi="McCormack_2013"
fi="Meiklejohn_2016"
#fi="Moyle_2016"

word="-permut-11"


ty="PF-mp"$word

pf=$fi$word"_entropy_partition_finder.cfg"

pd="/home/vtagliacollo/GitHub/PartitionUCE/processed_data/Permutations/PartitionFinder/"
rd="/home/vtagliacollo/GitHub/PartitionUCE/processed_data/Permutations/UCE-subsets/"
nd="/home/vtagliacollo/GitHub/PartitionUCE/processed_data/Permutations/alignments/"

hm=$pd$fi"/"

ad=$hm$ty
mkdir $ad



cp $rd$fi"/"$fi$word"/"$pf  $ad"/partition_finder.cfg" # problem is here
cp $nd$fi"/"$fi$word".phy" $ad"/"$fi$word".phy"



#python /disks/dacelo/data/programs/github/partitionfinder/PartitionFinder.py $ad --raxml -n -p 10
python /disks/dacelo/data/programs/github/partitionfinder/PartitionFinder.py $ad"/partition_finder.cfg" $ad"/"$fi$word".phy" --raxml -n -p 10

##### IQ-Tree
name="Crawford_2012"
#name="Harrington_2016"
#name="McCormack_2013"
#name="Meiklejohn_2016"
#name="Moyle_2016"

phy=$name".phy"
nex=$name".nex"

cwd="/home/vtagliacollo/GitHub/PartitionUCE/processed_data/IQ-Tree/"
swd=$cwd$name
cd $swd

$swd"/bin/"iqtree-omp -s $phy -spp $nex -b 100 -nt 10 -wba

# https://groups.google.com/forum/#!searchin/iqtree/bootstraps/iqtree/EokGXmxYB1M/323etcJiBwAJ
# https://groups.google.com/forum/#!topic/iqtree/4MqPgFnSIDQ





