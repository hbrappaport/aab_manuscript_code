library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(scales)
library(ggpol)
library(forcats)
library(stringr)
library(tidyr)
library(broom)

#############################
#all 3 species clusters together
#############################
within_species_functions = read.csv("raw_data/strain_diversity_pan_summary.csv")

#remove unknown functions
withinspecies_filtered = filter(within_species_functions, Shortened_Category != "Function unknown")

#removes blank
withinspecies_filtered = filter(withinspecies_filtered, Shortened_Category != "0")

#removes general function prediction
withinspecies_filtered = filter(withinspecies_filtered, Shortened_Category != "General function prediction only")

grouped_bycat_all <- withinspecies_filtered %>%                                   
  group_by(Shortened_Category, genome_name, Species) %>%                   
  count(Shortened_Category, core_accessory) %>%          # counting each categorical variable
  mutate(percent = n/sum(n))                           # creating the "%" variable and values

grouped_bycat_all_accessory <- grouped_bycat_all %>%
  filter(core_accessory == "accessory")

grouped_bycat_all_shortened = grouped_bycat_all_accessory %>%
  filter(Shortened_Category == "Mobilome: prophages, transposons" | Shortened_Category 
         == "Carbohydrate transport and metabolism" | Shortened_Category ==
           "Translation, ribosomal structure and biogenesis" | Shortened_Category == "Inorganic ion transport and metabolism"|
         Shortened_Category == "Energy production and conversion")


###Fig 4B
core_acc_by_func_all_shortened <- 
  grouped_bycat_all_accessory_shortened %>%
  ggplot(data = grouped_bycat_all_shortened, mapping = aes(x = fct_reorder(Shortened_Category, -percent, .fun = mean, .desc =FALSE), y = percent, fill = Species)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, color = "#363636", lwd = 0.3) +
  geom_jitter(set.seed(666), show.legend = FALSE, pch = 21, size = 3, color = "#000000", stroke = 0.3)+
  scale_fill_manual(values = c("#5a738c", "#6c5a8c", "#7b5a8c")) +                           # coloring the plot
  labs(x = "Cog20 Category",                                              # labelling x axis
       y = "Percentage of Accessory Genes", # labeling y axis
       fill = "Core vs. Accessory") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 35,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  ) +
  theme(axis.text.x=element_text(size=8, angle=35,hjust=0.95,vjust=0.2)) +
  scale_x_discrete(labels = c("Mobilome", "Defense mechanisms", "Carbohydrate transport and metabolism", "Energy production and conversion", "Translation, ribosomal structure, and biogenesis"))+
  facet_grid("Species") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust =1, color = "black", family = "Arial")) +
  theme(axis.text.y = element_text(size = 8, color = "black", family = "Arial")) +
  theme(axis.title.y = element_text(family = "Arial")) +
  theme(axis.title.x = element_text(family = "Arial"))
core_acc_by_func_all_shortened


###Fig SX
core_acc_by_func_all <- 
  grouped_bycat_all_accessory %>%
  ggplot(data = grouped_bycat_all_accessory, mapping = aes(x = fct_reorder(Shortened_Category, -percent, .fun = mean, .desc =FALSE), y = percent, fill = Species)) +
  geom_boxplot(outlier.shape = NA, show.legend = FALSE, color = "#363636", lwd = 0.3) +
  geom_jitter(set.seed(666), show.legend = FALSE, pch = 21, size = 3, color = "#000000", stroke = 0.3)+
 scale_fill_manual(values = c("#5a738c", "#6c5a8c", "#7b5a8c")) +                           # coloring the plot
  labs(x = "Cog20 Category",                                              # labelling x axis
       y = "Percentage of Accessory Genes", # labeling y axis
       fill = "Core vs. Accessory") +                               # legend
  scale_y_continuous(labels = scales::percent_format()) +         # changing the y axis nber format
  theme(
    axis.text.x = element_text(angle = 35,                        # rotating the x axis text
                               vjust = 0.5),                      # adjusting the position
    axis.title.x = element_text(face = "bold"),                   # face the x axit title/label
    axis.title.y = element_text(face = "bold"),                   # face the y axis title/label
    plot.title = element_text(hjust = 0.5),                       # positioning the plot title
    legend.title = element_text(face = "bold")                    # face the legend title
  ) +
  theme(axis.text.x=element_text(size=8, angle=35,hjust=0.95,vjust=0.2)) +
  facet_grid("Species") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 8, angle = 30, hjust =1, color = "black", family = "Arial")) +
  theme(axis.text.y = element_text(size = 8, color = "black", family = "Arial")) +
  theme(axis.title.y = element_text(family = "Arial")) +
  theme(axis.title.x = element_text(family = "Arial"))
core_acc_by_func_all
