
# libraries
library(vegan)
library(tidyverse)

# read in data
ani_matrix = read.csv("raw_data/ani_matrix.csv", sep = ",", row.names = 1)
pH_matrix = read.csv("raw_data/pH_matrix.csv", sep = ",", row.names = 1)
cfu_matrix = read.csv("raw_data/cfu_matrix.csv", row.names = 1)
ph_meta = read.csv("raw_data/ph_metadata.csv")
cfus_meta = read.csv("raw_data/cfus_metadata.csv")

# convert to dist matrices
ani_matrix = as.dist(1-ani_matrix)
pH_matrix = as.dist(pH_matrix)
cfu_matrix = as.dist(cfu_matrix)

# mantel tests
set.seed(1)
mantel(xdis = ani_matrix, ydis = pH_matrix, method = "spearman", permutations = 9999)
mantel(xdis = ani_matrix, ydis = cfu_matrix, method = "spearman", permutations = 9999)

# corr tests
cor.test(cfus_meta$avg_pH, cfus_meta$CFU_A, method = "spearman")
cor.test(cfus_meta$avg_pH, cfus_meta$CFU_Y, method = "spearman")
cor.test(cfus_meta$avg_pH, cfus_meta$CFU_L, method = "spearman")


# plot CFUs vs. pH
ggplot(cfus_meta, aes(avg_pH, CFU_A)) +
  geom_point(fill="#3A62EC", pch=21, size=3,
             stroke=0.2) + theme_classic() +
  geom_text(hjust=-0.1, vjust=0, aes(label=strain_treatment), size =3) + 
  geom_smooth(span=1.2, se=F, color="black", lwd = 0.5)

# ANI vs. pH
dm_ani_ph = as.data.frame(cbind(ani_matrix, pH_matrix))

ggplot(dm_ani_ph, aes(ani_matrix, pH_matrix)) +
  geom_point(fill="#3A62EC", pch=21, size=3,
             stroke=0.2) +
  geom_smooth(span=1, se=F, color="black", lwd = 0.5) +
  theme_classic() + labs(x="ANI distance",
                         y="pH distance")

