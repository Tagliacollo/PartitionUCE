library(ggplot2)

ent = read.csv("~/Documents/github/PartitionUCE/entropy/Crawford_2012.csv")

# the picture acros all UCEs together
p1 = ggplot(ent, aes(x = site, y = entropy))
p1 + geom_smooth()

# the picture across all UCEs individually
p2 = ggplot(ent, aes(x = site, y = entropy, group = name))
p2 + geom_line(stat="smooth", method = "loess", size = 0.1, alpha = 0.5)

# look at individual UCEs with point data (shows there's not much info on a site-by-site basis)
for(i in sample(1:length(levels(ent$name)), 10)){
    p3 = ggplot(subset(ent, name == levels(ent$name)[i]), aes(x = site, y = entropy))
    plot(p3 + geom_smooth() + geom_point())
    cat ("Press [enter] to continue")
    line <- readline()
}