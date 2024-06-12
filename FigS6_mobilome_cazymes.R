library(ggplot2)
library(readr)
library(ggpattern)

#results of DRAMv and DRAM on plasmids
cazy_all <- read.table(file = 'raw_data/cazy_all.tsv', sep = '\t', header = TRUE)

#group by & count number of cazy best hits
cazy_all_group <- group_by(cazy_all, cazy_best_hit) %>%
  mutate(cazy_best_hit, n()) %>%
  rename("best_hit_count" = "n()")

#plot best hits colored by source (plasmid or prophage)
cazy_plot = ggplot(data=cazy_all_group, aes(y=cazy_best_hit, fill = source)) +
  geom_bar() +
  scale_fill_manual(values = c("#c6c9ea", "#4363ef")) +
  theme_classic() 

cazy_plot

#Plot Fig. S6, cazymes colored by species recovered from and patterened by source
ggplot(data = cazy_all_group, aes(y = cazy_best_hit, fill = species, pattern = source)) +
  geom_bar_pattern(
    color = "black", 
    #pattern_fill = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_spacing = 0.01,
    pattern_key_scale_factor = 0.6) + 
  scale_fill_manual(values = c(c("#4363EF","#909EDD", "#48558E", "#fceab4", "#CEAA3A", "#EEC33F", "#5a8c7b","#D1536B", "#B2374F", "#D5936C", "#9F5C35"))) +
  scale_pattern_manual(values = c(prophage = "stripe", plasmid = "none")) +
  #labs(x = "Class", y = "Number of Students", pattern = "Nerd?") + 
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none"))) +
  theme_classic()
