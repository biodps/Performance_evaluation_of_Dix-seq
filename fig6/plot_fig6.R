################################################################################ 
# This script has been used to visualize the statistical results of Figures 6
# in the manuscript "Dix-seq: An integrated pipeline for fast amplicon data analysis".  
# Pensgheng Dong @ Zhengzhou, China
# dpsh@henau.edu.cn 
################################################################################

#setting the pathway of the Working Directory
setwd("E:/YOUR_PATHWAY/fig6/")

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
                 ggpubr, ggsci, ggridges, fs, Rmisc, parallel)

# reading input data
for (i in dir()[grep("rds$",dir())]) {
  t1<- gsub("\\.rds", "", i)
  assign(t1, readRDS(i))
}
design1<- rbind(design, design %>% mutate(type="all"))

# PCoa variance explained for PC1, PC2, and PC3
all_pcoa_eig<- NA
for (a in ls()[grep("_[b,w][c,u]$", ls())]) {
  for (aa in unique(design1$type)) {
    a1<- subset(design1, type==aa)
    a2<- cmdscale(as.matrix(get(a)[a1$Run,a1$Run]), k=2, eig=T)
    a3<- data.frame(eig=round(a2$eig/sum(a2$eig)*100, 3),
                    dis=gsub("^.*_", "", a),
                    type=aa,
                    method=gsub("_.*$", "", a))[1:3,] %>% mutate(axis=paste0("PCoA",1:3))
    all_pcoa_eig<- rbind(all_pcoa_eig, a3)
    
  }
  
}

all_pcoa_eig<- all_pcoa_eig[complete.cases(all_pcoa_eig),]
all_pcoa_eig$dis1<- c("Bray-Curtis", "Weighted UniFrac")[match(all_pcoa_eig$dis, c("bc","wu"))]
all_pcoa_eig$dis1<- factor(all_pcoa_eig$dis1, levels = c("Bray-Curtis", "Weighted UniFrac"))
all_pcoa_eig$method1<- c("DADA2","Deblur","Dix-seq","EasyAmp","S97%","S99%")[match(all_pcoa_eig$method, unique(all_pcoa_eig$method))]
all_pcoa_eig$method1<- factor(all_pcoa_eig$method1,
                              levels = c("Dix-seq","EasyAmp","DADA2","Deblur","S97%","S99%"))
all_pcoa_eig$method2<- rep(c("Qiime2","Usearch","Qiime1"),c(2,2,2))[match(all_pcoa_eig$method, unique(all_pcoa_eig$method))]
all_pcoa_eig$method2<- factor(all_pcoa_eig$method2, levels = c("Usearch","Qiime2","Qiime1"))
all_pcoa_eig$pc<- paste0(all_pcoa_eig$dis, all_pcoa_eig$axis)
all_pcoa_eig$pc1<- factor(all_pcoa_eig$pc, levels = c(paste0("bcPCoA",3:1),paste0("wuPCoA",3:1)))
all_pcoa_eig$axis1<- factor(all_pcoa_eig$axis, levels = paste0("PCoA",3:1))
all_pcoa_eig$type1<- c("Larve","Water","All")[match(all_pcoa_eig$type,unique(design1$type))]
all_pcoa_eig$type1<- factor(all_pcoa_eig$type1, levels = c("Larve","Water","All"))
all_pcoa_eig$x<- c(1:6)[match(all_pcoa_eig$method1, levels(all_pcoa_eig$method1))]

# PERMANOVA pseudo-F statistics and R2 for the separation between samples from stage and habitats.

all_adonis_dis<- NA
for (b in ls()[grep("_[b,w][c,u]$", ls())]) {
  b1<- adonis2(as.dist(as.matrix(get(b))) ~ type + stage,
               data = design, permutations = 999, parallel=22)
  b2<- adonis2(as.dist(as.matrix(get(b))) ~ type * stage,
               data = design, permutations = 999, parallel=22)
  all_adonis_dis<- rbind(all_adonis_dis,
                         data.frame(method=rep(gsub("_.*", "", b),3),
                                    dis=rep(gsub(".*_", "", b),3),
                                    group=c("Habitat", "Stage", "Interaction"),
                                    R2=c(b1$R2[1:2],b2$R2[3]),
                                    Fvalue=c(b1$F[1:2],b2$F[3]),
                                    pvalue=c(b1[[5]][1:2],b2[[5]][3])))
  
}
all_adonis_dis<- all_adonis_dis[complete.cases(all_adonis_dis),]
all_adonis_dis$dis1<- c("Bray-Curtis", "Weighted UniFrac")[match(all_adonis_dis$dis, c("bc","wu"))]
all_adonis_dis$dis1<- factor(all_adonis_dis$dis1, levels = c("Bray-Curtis", "Weighted UniFrac"))
all_adonis_dis$method1<- c("DADA2","Deblur","Dix-seq","EasyAmp","S97%","S99%")[match(all_adonis_dis$method, unique(all_adonis_dis$method))]
all_adonis_dis$method1<- factor(all_adonis_dis$method1,
                                levels = c("Dix-seq","EasyAmp","DADA2","Deblur","S97%","S99%"))
all_adonis_dis$method2<- rep(c("Qiime2","Usearch","Qiime1"),c(2,2,2))[match(all_adonis_dis$method, unique(all_adonis_dis$method))]
all_adonis_dis$method2<- factor(all_adonis_dis$method2, levels = c("Usearch","Qiime2","Qiime1"))
all_adonis_dis$group1<- factor(all_adonis_dis$group, levels = c("Habitat", "Stage", "Interaction"))

# Ratio of Nauplius stage vs. others.

nauplius_id<- design[grep("N",design$stage),1]
all_mratios_dis<- NA
for (c in ls()[grep("_[b,w][c,u]$", ls())]) {
  c1<- get(c)
  for (c2 in nauplius_id){
    #c2 sample type: c3
    c3<- design[which(design$Run %in% c2), 3]
    #non-nauplius samples in the same habitat type of c2
    c4<- intersect(design[which(!design$Run %in% nauplius_id), 1],
                   design[which(design$type %in% c3), 1])
    c5<- c1[c4,c2]
    #non-nauplius samples in the different habitat type of c2
    c6<- intersect(design[which(!design$Run %in% nauplius_id), 1],
                   design[which(!design$type %in% c3), 1])
    c7<- c1[c6,c2]
    #non-nauplius samples in all two habitat types of c2
    c8<- c1[design[which(!design$Run %in% nauplius_id), 1],c2]
    all_mratios_dis<- rbind(all_mratios_dis,
                            data.frame(sample=c2,ratio=1-c5,habtat=c3,
                                       method=gsub("_.*", "", c),
                                       dis=gsub(".*_", "", c),
                                       treat="Proximity in same habtat"),
                            data.frame(sample=c2,ratio=1-c7,habtat=c3,
                                       method=gsub("_.*", "", c),
                                       dis=gsub(".*_", "", c),
                                       treat="Proximity in different habtat"),
                            data.frame(sample=c2,ratio=1-c8,habtat=c3,
                                       method=gsub("_.*", "", c),
                                       dis=gsub(".*_", "", c),
                                       treat="Proximity in all habtats"))
  }
  
}

all_mratios_dis<- all_mratios_dis[complete.cases(all_mratios_dis),]
all_mratios_dis$dis1<- c("Bray-Curtis", "Weighted UniFrac")[match(all_mratios_dis$dis, c("bc","wu"))]
all_mratios_dis$dis1<- factor(all_mratios_dis$dis1, levels = c("Bray-Curtis", "Weighted UniFrac"))
all_mratios_dis$method1<- c("DADA2","Deblur","Dix-seq","EasyAmp","S97%","S99%")[match(all_mratios_dis$method, unique(all_mratios_dis$method))]
all_mratios_dis$method1<- factor(all_mratios_dis$method1,
                                 levels = c("Dix-seq","EasyAmp","DADA2","Deblur","S97%","S99%"))
all_mratios_dis$method2<- rep(c("Qiime2","Usearch","Qiime1"),c(2,2,2))[match(all_mratios_dis$method, unique(all_mratios_dis$method))]
all_mratios_dis$method2<- factor(all_mratios_dis$method2, levels = c("Usearch","Qiime2","Qiime1"))
all_mratios_dis$treat1<- factor(all_mratios_dis$treat,
                                levels = c("Proximity in same habtat", "Proximity in different habtat","Proximity in all habtats"))


#set colors
col1<- c(brewer.pal(9,"Blues")[c(9,7,5)], brewer.pal(9,"Oranges")[c(6,5,3)])
names(col1)<- unique(all_pcoa_eig$pc)
col2<- col1[1:3]
names(col2)<- paste0("PCoA",1:3)
col3<- col1[c(1,4)]
names(col3)<- c("Bray-Curtis", "Weighted UniFrac")

# plot Fig6 by using ggplt2
p1<- ggplot() +
  geom_bar(data=subset(all_pcoa_eig, dis=="bc"), 
           aes(x=x, y=eig, fill = pc1), 
           stat = "identity",position = "stack",width=0.43) +
  geom_bar(data=subset(all_pcoa_eig, dis=="wu"), 
           aes(x=x+0.3+0.1, y=eig, fill = pc1), 
           stat = "identity",position = "stack",width=0.43) +
  facet_grid( ~ type1, space = "free_x", scales = "free_x") +
  scale_y_continuous(expand=c(0, 0.05))+
  scale_x_continuous(breaks = c(1.2,2.2,3.2,4.2,5.2,6.2),
                     labels = levels(all_pcoa_eig$method1))+
  labs(y = '% variance explained') +
  scale_fill_manual(values = col1) +
  theme_bw() +
  # rotate_x_text(angle = 60) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'bottom',
        strip.background = element_rect(color=NA, fill=NA),
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank()) +
  guides(fill=guide_legend(reverse = T, nrow = 1))

p2<- ggplot(data=all_adonis_dis) +
  geom_bar(aes(interaction(method1, method2), y=Fvalue, group = dis1, fill = dis1), 
           stat = "identity",position = "dodge") +
  facet_grid( ~ group1, space = "free_x", scales = "free_x")  +
  guides(x = "axis_nested") +
  scale_y_continuous(expand=c(0, 0.01))+
  labs(y = 'pseudo-F') +
  scale_fill_manual(values = col3) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        strip.background = element_rect(color=NA, fill=NA),
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank()) 

p3<- ggplot(all_mratios_dis, 
            aes(interaction(method1, method2), y=ratio, color=dis1, fill=dis1)) +
  geom_boxplot(alpha = .75) +
  facet_grid( ~ treat1, space = "free_x", scales = "free_x")  +
  geom_hline(yintercept=1, color="black", linetype="dashed") + 
  annotate("text", x=2, y=0.985, label="no diff.") +
  guides(x = "axis_nested") +
  scale_y_continuous(expand=c(0, 0.01))+
  labs(y = 'Ratio of nauplius stage vs. others') +
  scale_fill_manual(values = col3) +
  scale_color_manual(values = col3) +
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.position = 'none',
        strip.background = element_rect(color=NA, fill=NA),
        panel.border = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_blank()) 

p_comp<- ggarrange(p1,p2,p3,ncol = 1, nrow = 3, labels = LETTERS[1:3],
                   align="v", legend = "none",heights=c(4,3,4))

ggsave(filename = "Fig6_beta_comp.pdf", p_comp, width = 10.5, height = 8)


