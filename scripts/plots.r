library(ggplot2)


load.dataset = function(dataset_name){

  base_file = "~/Documents/github/PartitionUCE/processed_data/"
  
  ent = read.csv(paste(base_file, dataset_name, "_entropy.csv", sep = ""))
  tig = read.csv(paste(base_file, dataset_name, "_tiger.csv", sep = ""))
  gc  = read.csv(paste(base_file, dataset_name, "_gc.csv", sep = ""))
                                
  ent$type = "Entropy"
  gc$type = "%GC"
  tig$type = "TIGER"
  
  dataset = rbind(ent, gc, tig)

  return(dataset)
}

Faircloth = load.dataset("Faircloth_2013")
Faircloth$dataset = "Faircloth_2013"
McCormack = load.dataset("McCormack_2013")
McCormack$dataset = "McCormack_2013"
Meiklejohn = load.dataset("Meiklejohn_2016")
Meiklejohn$dataset = "Meiklejohn_2016"
Crawford = load.dataset("Crawford_2012")
Crawford$dataset = "Crawford_2012"
Moyle = load.dataset("Moyle_2016")
Moyle$dataset = "Moyle_2016"



# bind it all up and reorder the factor
data = rbind(Faircloth, McCormack, Meiklejohn, Crawford, Moyle)

# the picture acros all UCEs together
p1 = ggplot(data, aes(x = site, y = value))
p1 + geom_smooth() + facet_grid(type~dataset, scales = "free")

# the picture across all UCEs individually
p2 = ggplot(data, aes(x = site, y = value, group = name))
p2 + geom_line(stat="smooth", method = "loess", size = 0.1, alpha = 0.2) + facet_grid(type~dataset, scales = 'free')

