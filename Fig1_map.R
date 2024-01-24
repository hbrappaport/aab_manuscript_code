library(ggplot2)
library(tidyverse)
library(maps)

####code adapted from https://sarahleejane.github.io/learning/r/2014/09/20/plotting-beautiful-clear-maps-with-r.html####

# make world map using map_data() function
world_map <- map_data("world")

#import data set with all 500 sourdough starter coordinates
map_shape_500 <- read.csv("raw_data/aab_shape_500.csv", header = TRUE, stringsAsFactors = FALSE)

#import data set with all AAB amplicon and genome coordinates
map_color_500 <- read.csv("raw_data/aab_color_500.csv", header = TRUE, stringsAsFactors = FALSE)

#Create a base plot with gpplot2
p <- ggplot() + coord_fixed() +
  xlab("") + ylab("")

#Add map to base plot
base_world_messy <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), 
                                     colour="#CACFD2", fill="#CACFD2")

base_world_messy

cleanup <- 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill = 'white', colour = 'white'), 
        axis.line = element_line(colour = "white"), legend.position="right",
        axis.ticks=element_blank(), axis.text.x=element_blank(),
        axis.text.y=element_blank())

base_world <- base_world_messy + cleanup

base_world


map_plot <-
  base_world +
  geom_point(data = map_shape_500, aes(x = CurrentLong, y = CurrentLat, pch = Shape),
             size = 2.5, stroke = 0.2) +
  scale_shape_manual(values = c(1, 21, 1, 21)) +
  geom_point(data = map_color_500, aes(x = CurrentLong, y = CurrentLat, fill = Color),
             size = 2.5, stroke = 0.2, pch = 21, color = "black") +
  scale_fill_manual(values = c("#a2a9c8", "#d5936c","#4363ef" ))

map_plot