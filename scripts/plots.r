library('reshape2')
library(plotly)


load.dataset = function(dataset_name){

  base_file = "~/Documents/github/PartitionUCE/processed_data/"
  
  ent = read.csv(paste(base_file, dataset_name, "_entropy.csv", sep = ""))
  gc  = read.csv(paste(base_file, dataset_name, "_gc.csv", sep = ""))

  ent_min = ddply(ent, 'name', summarise, min = min(site))
  gc_min = ddply(gc, 'name', summarise, min = min(site))
  
  ent_windows = read.csv(paste(base_file, dataset_name, "_entropy_windows.csv", sep = ""))
  gc_windows  = read.csv(paste(base_file, dataset_name, "_gc_windows.csv", sep = ""))
  
  ent_w = merge(ent_min, ent_windows, by = 'name')
  ent_w$start = ent_w$min + ent_w$start
  ent_w$stop = ent_w$min + ent_w$stop
  ent_w = ent_w[,-2]
  gc_w = merge(gc_min, gc_windows, by = 'name')
  gc_w$start = gc_w$min + gc_w$start
  gc_w$stop = gc_w$min + gc_w$stop
  gc_w = gc_w[,-2]

  ent$type = "Entropy"
  gc$type = "%GC"

  ent = merge(ent, ent_w, by = 'name')
  gc = merge(gc, gc_w, by = 'name')
  
  dataset = rbind(ent, gc)

  dataset$region = 'core'
  dataset$region[which(dataset$site < dataset$start)] = 'left_flank'
  dataset$region[which(dataset$site > dataset$stop)] = 'right_flank'
  
  dataset$dataset = dataset_name
  
  return(dataset)
}

Faircloth = load.dataset("Faircloth_2013")
McCormack = load.dataset("McCormack_2013")
Meiklejohn = load.dataset("Meiklejohn_2016")
Crawford = load.dataset("Crawford_2012")
Moyle = load.dataset("Moyle_2016")
Smith = load.dataset("Smith_2014")



# bind it all up and reorder the factor
data = rbind(Faircloth, McCormack, Meiklejohn, Crawford, Moyle, Smith)

# the picture acros all UCEs together
p1 = ggplot(data, aes(x = site, y = value))
p1 + geom_smooth() + facet_grid(type~dataset, scales = "free")

# the picture across all UCEs individually
p2 = ggplot(data, aes(x = site, y = value, group = name, color = region))
p2 +  geom_line(stat="smooth", method = "loess", size = 0.1, alpha = 0.2) + 
      facet_grid(type~dataset, scales = 'free') +
      scale_color_identity()


# we can look at individual UCEs like this:
dataset = Crawford
names = levels(dataset$name)
one_uce = subset(dataset, name == sample(names, 1))
p3 = ggplot(one_uce, aes(x = site, y = value, color = region))
p3 + geom_smooth() + geom_point() + facet_wrap(~type, ncol = 1, scales = "free")


# we can look at the lengths of the windows like this
data$length = data$stop - data$start
lengths = ddply(data, .(name, dataset, type), summarise, length = max(length))
p4 = ggplot(lengths, aes(x = length))
p4 + geom_histogram() + facet_grid(type~dataset, scales = "free")


################### VICTOR
library(ggplot2)

load.data = function(dataset_name){
  
  base_path = '/Volumes/VATagliacollo/GitHub/PartitionUCE/processed_data/'
  data_path = paste0(base_path, dataset_name, '/', dataset_name)
  
  df    = read.csv(paste0(data_path, ".csv", sep = ""))
  
  return(df)
  
}

df = load.data('test')


p = ggplot(df, aes(uce_site, name)) + 
  geom_tile(aes(fill = plot_mtx), colour = "white") + 
  scale_fill_gradient(low = "black", high = "steelblue") + 
  #facet_grid(type ~ .) # version 1  
  facet_wrap(~ type, ncol = 2) # vesion 2
  
base_size = 9
p + theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
        legend.position = "none")




