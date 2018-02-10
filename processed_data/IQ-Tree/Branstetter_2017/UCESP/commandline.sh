inputa="*.phy"
inputp="*.nex"
threads=20

./iqtree-omp -s $inputa -sp $inputp -nt $threads -b 100
