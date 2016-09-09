library(ggplot2)

ent = read.csv("~/Documents/github/PartitionUCE/processed_data/Crawford_2012_entropy.csv")
gc  = read.csv("~/Documents/github/PartitionUCE/processed_data/Crawford_2012_gc.csv")
tig = read.csv("~/Documents/github/PartitionUCE/processed_data/Crawford_2012_tiger.csv")

ent$type = "entropy"
gc$type = "gc"
tig$type = "tiger"

crawford = rbind(ent, gc, tig)

# the picture acros all UCEs together
p1 = ggplot(crawford, aes(x = site, y = value))
p1 + geom_smooth() + facet_wrap(~type, scales = "free_y", ncol = 1)


# the picture across all UCEs individually
p2 = ggplot(crawford, aes(x = site, y = value, group = name))
p2 + geom_line(stat="smooth", method = "loess", size = 0.1, alpha = 0.3) + facet_wrap(~type, scales = "free_y", ncol = 1)

# look at individual UCEs with point data (shows there's not much info on a site-by-site basis)
for(i in sample(1:length(levels(ent$name)), 10)){
  p3 = ggplot(subset(ent, name == levels(ent$name)[i]), aes(x = site, y = entropy))
  plot(p3 + geom_smooth() + geom_point())
  cat ("Press [enter] to continue")
  line <- readline()
}

g3 = ggplot(subset(gc, name == "chrz_14545"), aes(x = site, y = gc))
g3 + geom_smooth() + geom_point()
