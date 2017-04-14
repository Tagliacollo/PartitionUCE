library(treespace)
library(ape)
library(phangorn)

Single="/Users/roblanfear/Documents/github/PartitionUCE/processed_data/IQ-Tree/Crawford_2012/Single/Crawford_2012.nex.treefile"
SWSC_EN="/Users/roblanfear/Documents/github/PartitionUCE/processed_data/IQ-Tree/Crawford_2012/SWSC-EN/Crawford_2012.nex.treefile"

Single = read.tree(Single)
SWSC_EN = read.tree(SWSC_EN)

cophyloplot(SWSC_EN, Single)
sum(Single$edge.length)
sum(SWSC_EN$edge.length)

# if we root them first
plotTreeDiff(SWSC_EN,Single,use.edge.length=FALSE)
