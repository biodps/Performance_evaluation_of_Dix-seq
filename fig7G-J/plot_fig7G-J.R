################################################################################ 
# This script has been used to visualize the statistical results of Figures 7G-J
# in the manuscript "Dix-seq: An integrated pipeline for fast amplicon data analysis".  
# Pensgheng Dong @ Zhengzhou, China
# dpsh@henau.edu.cn 
################################################################################

#setting the pathway of the Working Directory
setwd("E:/YOUR_PATHWAY/fig7G-J/")

#loading packages and functions
tryCatch({
  library(librarian)
}, error=function(e) {
  install.packages("librarian")
  library(librarian)
})

librarian::shelf(ggplot2, tidyverse, RColorBrewer, ggh4x, ggpubr)

# reading input data
for (i in dir()[grep("rds$",dir())]) {
  t1<- gsub("\\.rds", "", i)
  assign(t1, readRDS(i))
}

# calculating the average of F values
IQR_PF_depth<- aggregate(value ~ method1 + index1, 
                         data = subset(PF_depth_df, depth1=="30k"), median)

# ploting Fig7 by using ggplot2
p1<- ggplot(subset(PF_depth_df,method1=="dix-seq"), aes(x=depth1, y=value,colour = index1)) +
  geom_point(alpha = 0.25) +
  geom_boxplot(outlier.shape =  NA) +
  geom_hline(data = subset(IQR_PF_depth,method1=="dix-seq"), aes(yintercept = value),
             linetype="dashed", linewidth=1.2, color="#ffa567") +
  facet_wrap( ~ index1, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("#000080", "#c80000")) +
  theme_bw() +
  # rotate_x_text(angle = 45) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(color=NA, fill=NA),
        panel.border = element_rect(colour = "black"),
        legend.position = 'none',
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black")) +
  labs(x="Sequencing depth", y=" Pseudo F")

p2<- ggplot(subset(Pearson_depth_df, method1=="dix-seq"), 
            aes(x=depth1, y=value, colour = index1)) +
  geom_point(alpha = 0.25) +
  geom_boxplot(outlier.shape =  NA) +
  facet_wrap( ~ index1 , ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("#2ca02c", "#36cabb")) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(color=NA, fill=NA),
        panel.border = element_rect(colour = "black"),
        legend.position = 'none',
        axis.text = element_text(colour = "black"),
        axis.title = element_text(colour = "black")) +
  labs(x="Sequencing depth", y="Pearson r")

p_dixsq<- ggarrange(p2,p1, nrow = 2, align = "v", common.legend = F)
ggsave(filename = "fig7G-J_dixsq_depth_box.pdf", p_dixsq, width = 8, height = 4.7)



