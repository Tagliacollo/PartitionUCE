library(ape)
library(phangorn)
library(phytools)
dataset="Harrington_2016"

t1=read.tree(paste("/Users/roblanfear/Documents/github/PartitionUCE/processed_data/IQ-Tree/", dataset, "/SWSC-EN/", dataset, ".nex.treefile", sep=""))
t2=read.tree(paste("/Users/roblanfear/Documents/github/PartitionUCE/processed_data/IQ-Tree/", dataset, "/SWSC-GC/", dataset, ".nex.treefile", sep=""))
t3=read.tree(paste("/Users/roblanfear/Documents/github/PartitionUCE/processed_data/IQ-Tree/", dataset, "/SWSC-MU/", dataset, ".nex.treefile", sep=""))
t4=read.tree(paste("/Users/roblanfear/Documents/github/PartitionUCE/processed_data/IQ-Tree/", dataset, "/UCESP/", dataset, ".nex.treefile", sep=""))
t5=read.tree(paste("/Users/roblanfear/Documents/github/PartitionUCE/processed_data/IQ-Tree/", dataset, "/PF-UCE/", dataset, ".nex.treefile", sep=""))
t6=read.tree(paste("/Users/roblanfear/Documents/github/PartitionUCE/processed_data/IQ-Tree/", dataset, "/UCE/", dataset, ".nex.treefile", sep=""))
t7=read.tree(paste("/Users/roblanfear/Documents/github/PartitionUCE/processed_data/IQ-Tree/", dataset, "/Single/", dataset, ".nex.treefile", sep=""))

assoc<-cbind(t1$tip.label,t1$tip.label)
plot(cophylo(t1, t2, assoc))
plot(cophylo(t1, t3, assoc))
plot(cophylo(t1, t4, assoc))
plot(cophylo(t1, t5, assoc))
plot(cophylo(t1, t6, assoc))
plot(cophylo(t1, t7, assoc))

sum(t1$edge.length) / sum(t2$edge.length)
sum(t1$edge.length) / sum(t3$edge.length)
sum(t1$edge.length) / sum(t4$edge.length)
sum(t1$edge.length) / sum(t5$edge.length)
sum(t1$edge.length) / sum(t6$edge.length)
sum(t1$edge.length) / sum(t7$edge.length)

path.dist(t1, t2)
path.dist(t1, t3)
path.dist(t1, t4)
path.dist(t1, t5)
path.dist(t1, t6)
path.dist(t1, t7)


