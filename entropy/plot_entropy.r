library(ggplot2)

ent = read.csv("~/Documents/github/PartitionUCE/entropy/Crawford_2012.csv")

p1 = ggplot(ent, aes(x = site, y = entropy))
p1 + geom_smooth()

# look at individual
p2 = ggplot(subset(ent, name == levels(ent$name)[80]), aes(x = site, y = entropy))
p2 + geom_smooth() + geom_point()
