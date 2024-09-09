################################################################################ 
# This script has been used to visualize the statistical results of Figures 4, and 7A-F
# in the manuscript "Dix-seq: An integrated pipeline for fast amplicon data analysis".  
# Pensgheng Dong @ Zhengzhou, China
# dpsh@henau.edu.cn 
################################################################################

#setting the pathway of the Working Directory
setwd("E:/YOUR_PATHWAY/fig4_7A-F/")

#loading packages and functions
tryCatch({
  library(librarian)
  library(BiocManager)
}, error=function(e) {
  install.packages(c("BiocManager","librarian"))
  library(librarian)
  library(BiocManager)
})

librarian::shelf(ggplot2, vegan, ade4, ini, knitr, tidyverse, RColorBrewer, ape,
                 cowplot, reshape2, ggdendro,  ggpubr, tidygraph, Hmisc, ggh4x,
                 SYNCSA, ggpubr, ggsci, ggridges, fs, Rmisc)

# reading input data
for (i in dir()[grep("rds$",dir())]) {
  t1<- gsub("\\.rds", "", i)
  assign(t1, readRDS(i))
}

#caculating parameters for performance evaluation
recover_df$perp<- c(recover_df$reads_mapped_to_features / recover_df$sequencing_depth) * 100
cr_97_stats$pert<- c(cr_97_stats$cr_97_counts / cr_97_stats$q1s97_counts) * 100
asvs_software_id<- data.frame(software=unique(asvs_count_df$software)[c(3,4,6,5,2,1)],
                              software1=factor(unique(asvs_count_df$software)[c(3,4,6,5,2,1)],
                                               levels = unique(asvs_count_df$software)[c(3,4,6,5,2,1)]),
                              loc1=1:6, loc=3) %>% 
  mutate(median=aggregate(cumulative_abundance ~ software1, data = asvs_count_df, median)[,2]) %>% 
  bind_cols(aggregate(cumulative_abundance ~ software1, data = asvs_count_df, FUN = function(x){CI(x, ci=0.99)})[,-1])
asvs_software_id$anno<- paste0(asvs_software_id$software,": Mean = ",round(asvs_software_id$mean,1),", CI99% = [",round(asvs_software_id$lower,1),", ",round(asvs_software_id$upper,1),"]")


# set colors
col1<- brewer.pal(12,"Paired")[6:1]
names(col1)<- c("Dix-seq","EasyAmp","DADA2","Deblur","S97%","S99%")
col2<- col1
names(col2)[5]<- "Wang et al"

# ploting Figures 4 by using ggplot2
p1<- ggplot(data=subset(run_time_df, subid=="30k"),
            aes(interaction(software1, method), y=time1, fill = protocols1)) +
  geom_bar(stat = "identity",position = "stack") +
  scale_color_bmj() +
  scale_fill_bmj() +
  guides(x = "axis_nested")+
  labs(y = 'Log10[Run time (second)]') +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank()) +
  rotate_x_text(30)+
  guides(color=guide_legend(reverse = T, nrow = 1),
         fill=guide_legend(reverse = T, nrow = 1))

p2<- ggplot(data=subset(recover_df, depth=="30k"),
            aes(interaction(software1, method), y=perp, color = software1, fill = software1)) +
  geom_boxplot(alpha = .85) +
  scale_color_manual(values = col1) +
  scale_fill_manual(values = col1) +
  guides(x = "axis_nested") +
  labs(y = 'Mappability (%)') +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank()) +
  rotate_x_text(30)

p3<- ggplot(data=subset(cr_97_stats, depth=="30k"),
            aes(interaction(software1, method), y=pert, color = software1, fill = software1)) +
  geom_boxplot(alpha = .85) +
  scale_color_manual(values = col1[-5]) +
  scale_fill_manual(values = col1[-5]) +
  guides(x = "axis_nested") +
  labs(y = 'Mapping recovery (%)') +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank()) +
  rotate_x_text(30)

p4<- ggplot(data = asvs_count_df, aes(x=log_abu, y = loc1)) + 
  geom_text(data=asvs_software_id, size = 2.5, 
            aes(x=loc, y=loc1 + 0.73, label=anno, colour = software1)) + 
  geom_density_ridges(aes(color=software1, fill=software1),alpha=0.47, quantile_lines = T) +
  xlab("Log10 (Feature Count number)") +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(colour = "black"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

p5<- ggplot(data=subset(tax_diff_sw_rank, depth=="30k"),
            aes(interaction(software1, method), y=tax_diff_pct, fill = rank1)) +
  geom_bar(stat = "identity",position = "stack") +
  scale_fill_npg() +
  guides(x = "axis_nested")+
  labs(y = 'Different taxonomy (%)') +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank()) +
  rotate_x_text(30)+
  guides(fill=guide_legend(reverse = T))

p6<- ggplot(data=subset(unclassified_pct, depth=="30k"),
            aes(interaction(software1, method), y=pct, fill = rank1)) +
  geom_bar(stat = "identity",position = "stack") +
  scale_fill_npg() +
  guides(x = "axis_nested")+
  labs(y = 'Unclassified feature (%)') +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank()) +
  rotate_x_text(30)+
  guides(fill=guide_legend(reverse = T))

p_stat1<- ggarrange(p1,p2,p3,p4,p6,p5, ncol = 3, nrow = 2, labels = LETTERS[1:6],
                    align="hv", legend = "bottom", common.legend = T)

p7<- ggplot(data=subset(alpha_div_df, depth=="30k"),
            aes(x = software1, y = div, color = software1, fill = software1)) +
  geom_violin(alpha = .85, trim = T) +
  stat_summary(fun.data=mean_sdl,  geom="pointrange", color="grey27") +
  facet_wrap( ~ index, nrow = 1, scales = "free_y") +
  stat_compare_means(label = "p.signif",aes(label = after_stat(p.signif)),
                     method = "t.test", ref.group = "Wang et al",
                     symnum.args = list(cutpoints = c(0, 0.01, 0.05, Inf), 
                                        symbols = c("**", "*", "ns"))) +
  scale_color_manual(values = col2) +
  scale_fill_manual(values = col2) +
  labs(y = 'Alpha diversity values') +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        strip.background = element_rect(color=NA, fill=NA),
        legend.position = 'none',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank()) +
  rotate_x_text(30)


# ploting Figures 7A-F by using ggplot2
p21<- ggplot(data=subset(recover_df, software=="Dix-seq"),
             aes(x=depth1, y=perp, color = depth1)) +
  geom_boxplot() +
  scale_color_manual(values = colorRampPalette(brewer.pal(9,'Reds')[8:4])(14)) +
  labs(x="Sequencing depth", y = 'Mappability (%)') +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"))

p31<- ggplot(data=subset(cr_97_stats, software=="Dix-seq"),
             aes(x=depth1, y=pert, color = depth1)) +
  geom_boxplot() +
  scale_color_manual(values = colorRampPalette(brewer.pal(9,'Greens')[8:4])(14)) +
  labs(x="Sequencing depth", y = 'Mapping recovery (%)') +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"))


p51<- ggplot(data=subset(tax_diff_sw_rank, software=="Dix-seq"),
             aes(x=depth1, y=tax_diff_pct, fill = rank1)) +
  geom_bar(stat = "identity",position = "stack") +
  scale_fill_npg() +
  labs(x="Sequencing depth", y = 'Different taxonomy (%)') +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black")) +
  guides(fill=guide_legend(reverse = T))

p61<- ggplot(data=subset(unclassified_pct, software=="Dix-seq"),
             aes(x=depth1, y=pct, fill = rank1)) +
  geom_bar(stat = "identity",position = "stack") +
  scale_fill_npg() +
  labs(x="Sequencing depth", y = 'Unclassified feature (%)') +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black")) +
  guides(fill=guide_legend(reverse = T))

p71<- ggplot(data=subset(alpha_div_df, software=="Dix-seq" & index=="Evenness"),
             aes(x = depth1, y = div, color = depth1)) +
  geom_boxplot() +
  scale_color_manual(values = colorRampPalette(brewer.pal(9,'Blues')[8:4])(14)) +
  labs(x="Sequencing depth", y = 'Evenness') +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black")) 

p72<- ggplot(data=subset(alpha_div_df, 
                         software=="Dix-seq" & index=="Phylogenetic diversity"),
             aes(x = depth1, y = div, color = depth1)) +
  geom_boxplot(outlier.shape =  NA) +
  scale_color_manual(values = colorRampPalette(brewer.pal(9,'Purples')[8:4])(14)) +
  labs(x="Sequencing depth", y = 'Phylogenetic diversity') +
  theme_bw() +
  scale_y_continuous(limits = c(0, 125)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black")) 

p_stat2<- ggarrange(p21,p31,p51,p61,p71,p72, ncol = 2, nrow = 3, labels = LETTERS[1:6],
                    align="hv", legend = "none")

ggsave(filename = "Fig4A-F.pdf", p_stat1, height = 6, width = 10)
ggsave(filename = "Fig4G_adiv.pdf", p7, height = 2.5, width = 10)
ggsave(filename = "Fig7A-F.pdf", p_stat2, width = 8, height = 5)

