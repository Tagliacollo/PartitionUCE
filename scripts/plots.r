library(ggplot2)

ent = read.csv("~/Documents/github/PartitionUCE/processed_data/Crawford_2012_entropy.csv")
gc  = read.csv("~/Documents/github/PartitionUCE/processed_data/Crawford_2012_gc.csv")

names(ent) = c("name", "site", "value")
ent$type = "entropy"

names(gc) = c("name", "site", "value")
gc$type = "gc"

crawford = rbind(ent, gc)

# the picture acros all UCEs together
p1 = ggplot(ent, aes(x = site, y = entropy))
p1 + geom_smooth()

g1 = ggplot(gc, aes(x = site, y = gc))
g1 + geom_smooth()

# the picture across all UCEs individually
p2 = ggplot(ent, aes(x = site, y = entropy, group = name))
p2 + geom_line(stat="smooth", method = "loess", size = 0.1, alpha = 0.5)

g2 = ggplot(gc, aes(x = site, y = gc, group = name))
g2 + geom_line(stat="smooth", method = "loess", size = 0.1, alpha = 0.1)


# look at individual UCEs with point data (shows there's not much info on a site-by-site basis)
for(i in sample(1:length(levels(ent$name)), 10)){
  p3 = ggplot(subset(ent, name == levels(ent$name)[i]), aes(x = site, y = entropy))
  plot(p3 + geom_smooth() + geom_point())
  cat ("Press [enter] to continue")
  line <- readline()
}

  g3 = ggplot(subset(gc, name == "chrz_14545"), aes(x = site, y = gc))
  g3 + geom_smooth() + geom_point()
