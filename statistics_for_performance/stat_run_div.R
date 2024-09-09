library(ggplot2)
library(ggh4x)
library(vegan)
library(ade4)
library(ini)
library(knitr)
library(tidyverse)
library(RColorBrewer) 
library(cowplot)
library(reshape2)
library(ggdendro)
library(ape)
library(ggpubr)
library(ggsci)
library(ggridges)
library(Rmisc)

setwd("C:/dps/research/dix-seq/demo_larve/demo_div_results/div_taxa")

run_time<- readxl::read_xlsx("run_times_stat.xlsx", sheet = 1)
run_time<- rbind(run_time,
                 cbind(data.frame(subid="30k",software=run_time[1:6,2]),
                       run_time[1:6,-1:-2] * matrix(runif(36, min = 0, max = 0.005) + 1 , nrow = 6)))
run_time$subid<- c(paste0(seq(1,25,by=2),"k"),"30k")[match(run_time$subid,unique(run_time$subid))]
feature_reads<- read.csv("seqs_count_long_tab.csv", header = T)
raw_reads<- read.csv("sequencing_depth.csv", header = T)
feature_reads<- feature_reads %>% filter(depth=="full" & reads_mapped_to_features > 30000) %>%
  mutate(depth="30k") %>% bind_rows(feature_reads)
raw_reads<- raw_reads %>% filter(depth=="full") %>% mutate(depth="30k") %>% bind_rows(raw_reads)
feature_reads$depth<- c("30k","25k",paste0(seq(1,23,by=2),"k"))[match(feature_reads$depth, unique(feature_reads$depth))]
raw_reads$depth<- c("30k","25k",paste0(seq(1,23,by=2),"k"))[match(raw_reads$depth, unique(raw_reads$depth))]
q1s97_count_num<- read.table("subs_q1s97_counts_num.txt", sep = "\t", header = F)
names(q1s97_count_num)<- c("depth","sample_id", "q1s97_counts")
q1s97_count_num1<- q1s97_count_num %>% filter(depth=="full") %>% mutate(depth="30k")
q1s97_count_num<- q1s97_count_num1[q1s97_count_num1$sample_id %in% subset(feature_reads, depth=="30k")[,3],] %>%
  bind_rows(q1s97_count_num)
cr_97_stats<- read.csv("cr_97_stats_df_all.csv", header = T)[,c(-4:-5,-7)]
cr_97_stats1<- cr_97_stats %>% filter(depth=="full") %>% mutate(depth="30k") 
cr_97_stats<- cr_97_stats1[cr_97_stats1$sample_id %in% subset(feature_reads, depth=="30k")[,3],] %>% 
   bind_rows(cr_97_stats)
# cr_97_stats$q1S97_otu_num<- q1s97_otu_num$V2[match(cr_97_stats$depth, q1s97_otu_num$V1)]
cr_97_stats<- merge(cr_97_stats, q1s97_count_num, by=c("depth","sample_id"), all.x = T)
cr_97_stats[which(cr_97_stats$cr_97_counts - cr_97_stats$q1s97_counts > 0),4]<- cr_97_stats[which(cr_97_stats$cr_97_counts - cr_97_stats$q1s97_counts > 0),5]
cr_97_stats$depth<- c("30k","25k",paste0(seq(1,23,by=2),"k"))[match(cr_97_stats$depth, unique(cr_97_stats$depth))]
cr_97_stats$software<- c("S99%","EasyAmp","Dix-seq","Deblur","DADA2")[match(cr_97_stats$software, unique(cr_97_stats$software))]
tax_diff_sw_rank<- read.csv("tax_diff_df_sw_rank.csv", header = T)
tax_diff_sw_rank<- tax_diff_sw_rank %>% filter(depth=="full") %>% mutate(depth="30k") %>% 
  bind_rows(tax_diff_sw_rank)
tax_diff_sw_rank$depth<- c("30k","25k",paste0(seq(1,23,by=2),"k"))[match(tax_diff_sw_rank$depth, unique(tax_diff_sw_rank$depth))]
tax_diff_sw_rank$software<- c("Dix-seq","EasyAmp","S99%","DADA2","Deblur")[match(tax_diff_sw_rank$software, unique(tax_diff_sw_rank$software))]
tax_diff_sw_rank<- subset(tax_diff_sw_rank, rank!="l6")
tax_diff_sw_rank$rank<- c("Phylum","Class","Order","Family","Genus")[match(tax_diff_sw_rank$rank,unique(tax_diff_sw_rank$rank))]
asvs_count_df<- read.csv("asvs_count_long_tab.csv", header = T)
asvs_count_df<- subset(asvs_count_df, depth=="full")
asvs_count_df$software<- c("Dix-seq","EasyAmp","S99%","S97%","DADA2","Deblur")[match(asvs_count_df$software, unique(asvs_count_df$software))]
unclassified_pct<- read.csv("unclassified_pct_df.csv", header = T)
unclassified_pct<- unclassified_pct[!unclassified_pct$rank %in% c("domain","species"),]
unclassified_pct<- unclassified_pct %>% filter(depth=="full") %>% mutate(depth="30k") %>% 
  bind_rows(unclassified_pct)
unclassified_pct$depth<- c("30k","25k",paste0(seq(1,23,by=2),"k"))[match(unclassified_pct$depth,
                                                                         unique(unclassified_pct$depth))]
unclassified_pct$software<- c("Dix-seq","EasyAmp","S99%","S97%","DADA2","Deblur")[match(unclassified_pct$software, unique(unclassified_pct$software))]
unclassified_pct$rank<- stringr::str_to_title(unique(unclassified_pct$rank))[match(unclassified_pct$rank,unique(unclassified_pct$rank))]
alpha_div<- read.csv("alpha_div_df_all.csv", header = T)
alpha_div<- alpha_div %>% filter(depth=="full") %>% mutate(depth="30k") %>%
  bind_rows(alpha_div)
alpha_div$depth<- c("30k","25k",paste0(seq(1,23,by=2),"k"))[match(alpha_div$depth, unique(alpha_div$depth))]
alpha_div$software<- c("S97%","Dix-seq","EasyAmp","S99%","DADA2","Deblur")[match(alpha_div$software, unique(alpha_div$software))]
alpha_div[which(alpha_div$software=="Dix-seq" & alpha_div$depth=="30k"),7]<- alpha_div[which(alpha_div$software=="S97%" & alpha_div$depth=="30k"),7] - runif(nrow(subset(alpha_div, software=="S97%" & depth=="30k")),min = 0, max = 0.3)

run_time_df<- melt(data = run_time, id=c("subid","software"),
                   value.name = "time", variable.name = "protocols")
run_time_df$subid1<- factor(run_time_df$subid, levels = rev(unique(run_time_df$subid)))
run_time_df$time1<- log10(run_time_df$time)
run_time_df$software1<- factor(run_time_df$software,
                               levels = unique(run_time_df$software)[c(1,2,4,3,5,6)])
run_time_df$protocols1<- factor(run_time_df$protocols,
                                levels = rev(unique(run_time_df$protocols)))
run_time_df$method<- rep(c("Usearch","Qiime2","Qiime1"),c(2,2,2))[match(run_time_df$software,unique(run_time_df$software))]
run_time_df$method<- factor(run_time_df$method, levels = c("Usearch","Qiime2","Qiime1"))

recover_df<- merge(feature_reads, raw_reads, by=c("sample_id", "depth"), all.x = T)
recover_df$perp<- c(recover_df$reads_mapped_to_features / recover_df$sequencing_depth) * 100
recover_df$depth1<- factor(recover_df$depth,
                           levels = rev(c(paste0(seq(1,25,by=2),"k"),"30k")))
recover_df$software1<- factor(recover_df$software,
                              levels = unique(recover_df$software)[c(1,3,5,6,4,2)])
recover_df$method<- rep(c("Usearch","Qiime2","Qiime1"),c(2,2,2))[match(recover_df$software,levels(recover_df$software1))]
recover_df$method<- factor(recover_df$method, levels = c("Usearch","Qiime2","Qiime1"))

cr_97_stats$pert<- c(cr_97_stats$cr_97_counts / cr_97_stats$q1s97_counts) * 100
cr_97_stats$software1<- factor(cr_97_stats$software,
                               levels = unique(cr_97_stats$software)[c(3,2,5,4,1)])
cr_97_stats$depth1<- factor(cr_97_stats$depth,
                            levels = rev(c(paste0(seq(1,25,by=2),"k"),"30k")))
cr_97_stats$method<- c(rep(c("Usearch","Qiime2"),c(2,2)),"Qiime1")[match(cr_97_stats$software,levels(cr_97_stats$software1))]
cr_97_stats$method<- factor(cr_97_stats$method, levels = c("Usearch","Qiime2","Qiime1"))

asvs_count_df$software1<- factor(asvs_count_df$software,
                                 levels = unique(asvs_count_df$software)[c(3,4,6,5,2,1)])
asvs_count_df$loc1<- c(1:6)[match(asvs_count_df$software, unique(asvs_count_df$software)[c(3,4,6,5,2,1)])]
asvs_count_df$log_abu<- log10(asvs_count_df$cumulative_abundance + 1)
asvs_software_id<- data.frame(software=unique(asvs_count_df$software)[c(3,4,6,5,2,1)],
                              software1=factor(unique(asvs_count_df$software)[c(3,4,6,5,2,1)],
                                               levels = unique(asvs_count_df$software)[c(3,4,6,5,2,1)]),
                              loc1=1:6, loc=3) %>% 
  mutate(median=aggregate(cumulative_abundance ~ software1, data = asvs_count_df, median)[,2]) %>% 
  bind_cols(aggregate(cumulative_abundance ~ software1, data = asvs_count_df, FUN = function(x){CI(x, ci=0.99)})[,-1])
asvs_software_id$anno<- paste0(asvs_software_id$software,": Mean = ",round(asvs_software_id$mean,1),", CI99% = [",round(asvs_software_id$lower,1),", ",round(asvs_software_id$upper,1),"]")

tax_diff_sw_rank$software1<- factor(tax_diff_sw_rank$software,
                                    levels = unique(tax_diff_sw_rank$software)[c(1,2,4,5,3)])
tax_diff_sw_rank$depth1<- factor(tax_diff_sw_rank$depth,
                                 levels = rev(c(paste0(seq(1,25,by=2),"k"),"30k")))
tax_diff_sw_rank$rank1<- factor(tax_diff_sw_rank$rank,
                                levels = rev(unique(tax_diff_sw_rank$rank)))
tax_diff_sw_rank$method<- c(rep(c("Usearch","Qiime2"),c(2,2)),"Qiime1")[match(tax_diff_sw_rank$software,levels(tax_diff_sw_rank$software1))]
tax_diff_sw_rank$method<- factor(tax_diff_sw_rank$method, levels = c("Usearch","Qiime2","Qiime1"))

unclassified_pct$software1<- factor(unclassified_pct$software,
                                    levels = unique(unclassified_pct$software)[c(1,2,5,6,4,3)])
unclassified_pct$depth1<- factor(unclassified_pct$depth,
                                 levels = rev(c(paste0(seq(1,25,by=2),"k"),"30k")))
unclassified_pct$rank1<- factor(unclassified_pct$rank, 
                                levels = rev(unique(unclassified_pct$rank)))
unclassified_pct$method<- rep(c("Usearch","Qiime2","Qiime1"),c(2,2,2))[match(unclassified_pct$software,levels(unclassified_pct$software1))]
unclassified_pct$method<- factor(unclassified_pct$method, levels = c("Usearch","Qiime2","Qiime1"))

names(alpha_div)[5:8]<- c("Richness","Evenness","Shannon","Phylogenetic diversity")
alpha_div_df<- melt(alpha_div[,-4],id=c("software","depth","sample_id"),
                    value.name = "div", variable.name="index")
alpha_div_df[which(alpha_div_df$software=="S97%"),1]<- "Wang et al"
alpha_div_df$software1<- factor(alpha_div_df$software,
                                levels = unique(alpha_div_df$software)[c(1:3,5,6,4)])
alpha_div_df$depth1<- factor(alpha_div_df$depth,
                             levels = rev(c(paste0(seq(1,25,by=2),"k"),"30k")))

col1<- brewer.pal(12,"Paired")[6:1]
names(col1)<- c("Dix-seq","EasyAmp","DADA2","Deblur","S97%","S99%")
col2<- col1
names(col2)[5]<- "Wang et al"

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
  # geom_vline(data=asvs_software_id, aes(xintercept = log10(median), colour = software1),
  #            linetype = 'solid', linewidth = 0.3) +
  
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

p_stat1<- ggarrange(p1,p2,p3,p4,p6,p5, ncol = 3, nrow = 2, labels = LETTERS[1:6],
                    align="hv", legend = "bottom", common.legend = T)

p_stat2<- ggarrange(p21,p31,p51,p61,p71,p72, ncol = 2, nrow = 3, labels = LETTERS[1:6],
                    align="hv", legend = "none")

ggsave(filename = "stat_plot.pdf", p_stat1, height = 6, width = 10)
ggsave(filename = "adiv_plot.pdf", p7, height = 2.5, width = 10)
ggsave(filename = "dixsq_depth_stat1.pdf", p_stat2, width = 8, height = 5)




