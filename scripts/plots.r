library(ggplot2)


load.dataset = function(dataset_name){

  base_file = "~/Documents/github/PartitionUCE/processed_data/"
  
  ent = read.csv(paste(base_file, dataset_name, "_entropy.csv", sep = ""))
  tig = read.csv(paste(base_file, dataset_name, "_tiger.csv", sep = ""))
  gc  = read.csv(paste(base_file, dataset_name, "_gc.csv", sep = ""))
                                
  ent$type = "entropy"
  gc$type = "gc"
  tig$type = "tiger"
  
  dataset = rbind(ent, gc, tig)

  return(dataset)
}

Faircloth = load.dataset("Faircloth_2013")
Faircloth$dataset = "Faircloth_2013"
McCormack = load.dataset("McCormack_2013")
McCormack$dataset = "McCormack_2013"


# bind it all up and reorder the factor
data = rbind(Faircloth, McCormack)

# the picture acros all UCEs together
p1 = ggplot(data, aes(x = site, y = value))
p1 + geom_smooth() + facet_grid(type~dataset, scales = "free")

# the picture across all UCEs individually
p2 = ggplot(data, aes(x = site, y = value, group = name))
p2 + geom_line(stat="smooth", method = "loess", size = 0.1, alpha = 0.1) + facet_grid(type~dataset, scales = 'free')

# look at individual UCEs with point data (shows there's not much info on a site-by-site basis)
for(i in sample(1:length(levels(ent$name)), 10)){
  p3 = ggplot(subset(ent, name == levels(ent$name)[i]), aes(x = site, y = entropy))
  plot(p3 + geom_smooth() + geom_point())
  cat ("Press [enter] to continue")
  line <- readline()
}

g3 = ggplot(subset(gc, name == "chrz_14545"), aes(x = site, y = gc))
g3 + geom_smooth() + geom_point()
