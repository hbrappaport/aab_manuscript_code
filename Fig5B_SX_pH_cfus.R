library(dplyr)
library(ggplot2)

pH_YL2_day4 = read.csv("raw_data/pH_data_yl2_day4.csv", sep = ",")

cfus_data = read.csv("raw_data/pH_by_cfus_avg.csv", sep = ",")

no_yeast_sparse = read.csv("raw_data/deleted_yeast_sparse.csv", sep = ",")


####CFUs plot Fig 5B
#pivot longer
cfus_data_long = cfus_data %>%
  select("CFU_Y", "CFU_L", "CFU_A", "comm_type") %>%
  pivot_longer(-comm_type, names_to = "CFUs", values_to = "value")

ggplot(cfus_data_long, aes(CFUs, log(value), fill = CFUs)) + 
  geom_point(size = 5, pch = 21, color = "black") +
  scale_fill_manual(values = c("#4363ef","#d5936c", "#c6c9ea"))+
  facet_grid(~factor(comm_type, levels = c("control", "YL", "A8", "A1", "A6", "A9",
                                      "A3", "A5", "A4", "A10", "A2", "A7"))) +
  theme_bw()

#####pH plot Fig 5B
ggplot(pH_YL2_day4, aes(y = pH, x = factor(comm_type, levels = c("control", "YL", "A8", "A1", "A6", "A9",
                                                                 "A3", "A5", "A4", "A10", "A2", "A7"))))+
  geom_boxplot(color = "dark grey") +
  geom_jitter() +
  theme_bw()

#####All pH Fig SX
ggplot(no_yeast_sparse, aes(y = pH, x = factor(comm_day, levels = c("C_Day2", "C_Day4", "YL_Day2", "YL_Day4", "A8_Day2", "A8_Day4", "A1_Day2", "A1_Day4",
   "A6_Day2", "A6_Day4", "A9_Day2", "A9_Day4", "A3_Day2", "A3_Day4", "A5_Day2", "A5_Day4", "A4_Day2", "A4_Day4",
        "A10_Day2", "A10_Day4", "A2_Day2", "A2_Day4", "A7_Day2", "A7_Day4")), color = day, shape = community)) +
  scale_color_manual(values = c("#d5936c", "#4363ef"))+
  #scale_shape_manual(values = c("triangle", "circle"))
  geom_boxplot(outlier.shape = NA, color = "dark grey") +
  geom_jitter() +
  #facet_grid(~factor(comm_type, levels = c("control", "YL", "A8", "A1", "A6", "A9",
                                      #     "A3", "A5", "A4", "A10", "A2", "A7")))+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1), panel.background = element_rect(fill = "white", color = "grey50"), 
        panel.grid.minor.y = element_line(color = "light grey", linewidth = 0.1), panel.grid.major.y = element_line(color = "light grey", linewidth =  0.1))


##ANOVA 
##pH YL2 day 4
aov_pH_yl2 = aov(pH~comm_type,data=pH_YL2_day4)
summary(aov_pH_yl2) 
boxplot(pH~comm_type,data=pH_YL2_day4)
TukeyHSD(aov_pH_yl2)
plot(TukeyHSD(aov_pH_yl2))

##all pH exclusing sparse yeast
aov_pH_all = aov(pH~comm_type_full, data=no_yeast_sparse)
summary(aov_pH_all)
boxplot(pH~comm_type_full, data=no_yeast_sparse)
pH_tukey_all = TukeyHSD(aov_pH_all)
plot(TukeyHSD(aov_pH_all))
pH_tukey_all <-as.data.frame(pH_tukey_all[1])
