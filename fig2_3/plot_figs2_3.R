################################################################################ 
# This script has been used to visualize the outputs from dix-seq for Figures 2 and 3 
# in the manuscript "Dix-seq: An integrated pipeline for fast amplicon data analysis".  
# The part of visualization code was directly produced by the Dix-seq 
# Pensgheng Dong @ Zhengzhou, China
# dpsh@henau.edu.cn 
################################################################################

#setting the pathway of the Working Directory
setwd("E:/YOUR_PATHWAY/fig2_3/")

#loading packages and functions
tryCatch({
  library(librarian)
  library(BiocManager)
}, error=function(e) {
  install.packages(c("BiocManager","librarian"))
  library(librarian)
  library(BiocManager)
})

librarian::shelf(reshape2, ggplot2, ggtree, tidyverse, ggstatsplot, palmerpenguins, 
                 RColorBrewer, ggtext, gstat, vegan, ellipse, ggpubr, ggalluvial, 
                 ggrepel, ComplexHeatmap, colorspace, dendextend)

get_pro_ell <- function(ord_in, grp_in, ellipse_pro = 0.97){
  ## from pro$X, pro$Yrot
  require(plyr)
  obs <- data.frame(rbind(pro$X[, 1:2], pro$Yrot[, 1:2]))
  obs$Groups <- grp_in
  names(obs)[1:2] <- c('x', 'y')
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  ell <- ddply(obs, 'Groups', function(x) {
    if(nrow(x) <= 2) {
      return(NULL)
    }
    sigma <- var(cbind(x$x, x$y))
    mu <- c(mean(x$x), mean(x$y))
    ed <- sqrt(qchisq(ellipse_pro, df = 2))
    data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
  })
  names(ell)[2:3] <- c('x', 'y')
  ell <- ddply(ell, .(Groups), function(x) x[chull(x$x, x$y), ])
  ell
}

# reading input data
meta_info<- readRDS('figs2_3_meta_info.rds')
adiv_df<- readRDS('figs2_3_adiv_df.rds')
otu_bc<- readRDS('figs2_3_otu_bc.rds')
otu_wu<- readRDS('figs2_3_otu_wu.rds')
fun_bc<- readRDS('figs2_3_fun_bc.rds')
phylum_top10<- readRDS('figs2_3_phylum_top10.rds')
top20_path<- readRDS('figs2_3_top12_path.rds')
diff_zotu_m2_3<- readRDS('figs2_3_diff_zotu.rds')

# pcoa for zotu and ko 
pcoa_otu <- cmdscale(as.matrix(otu_bc), k=2, eig=T)
pcoa_fun <- cmdscale(as.matrix(fun_bc), k=2, eig=T)

pcoa_otu_df <- data.frame(pc1    = pcoa_otu$points[,1],
                          pc2    = pcoa_otu$points[,2],
                          sample = rownames(pcoa_otu$points),
                          group  = meta_info$group[match(rownames(pcoa_otu$points),meta_info$Run)])
pcoa_otu_df$group1<- factor(pcoa_otu_df$group,level=unique(meta_info$group)[c(2,1,3,4)])
pcoa_otu_eig <- round(pcoa_otu$eig/sum(pcoa_otu$eig)*100, 2)

pcoa_fun_df <- data.frame(pc1    = pcoa_fun$points[,1],
                          pc2    = pcoa_fun$points[,2],
                          sample = rownames(pcoa_fun$points),
                          group  = meta_info$group[match(rownames(pcoa_fun$points),meta_info$Run)])

pcoa_fun_df$group1<- factor(pcoa_fun_df$group,level=unique(meta_info$group)[c(2,1,3,4)])

pcoa_fun_eig <- round(pcoa_fun$eig/sum(pcoa_fun$eig)*100, 2)

ellipse_otu <- NA
segment_otu <- NA
ellipse_fun <- NA
segment_fun <- NA
conf<- 0.95
uniques<- meta_info[,c(1,3)] %>% count(group) %>% filter(n > 1)

centroids_otu<- aggregate(cbind(pc1,pc2)~group1, pcoa_otu_df, mean) %>% filter(group1 %in%  uniques$group)
pcoa_ellipse_otu<- pcoa_otu_df %>% filter(group1 %in%  uniques$group)
ellipse_otu<- do.call(rbind, lapply(unique(pcoa_ellipse_otu$group1), function(t)
  data.frame(group = factor(as.character(t)),
             ellipse(cov(pcoa_ellipse_otu[pcoa_ellipse_otu$group==t,1:2]),
                     centre=as.matrix(centroids_otu[t,2:3]),
                     level=conf), stringsAsFactors=FALSE)))
ellipse_otu$group1 <- factor(ellipse_otu$group, 
                             level=unique(meta_info$group)[c(2,1,3,4)], ordered = TRUE)
segment_otu<- right_join(centroids_otu, pcoa_otu_df, by='group1')
colnames(segment_otu)[2:5] <- c("x","y","xend","yend")
segment_otu$group1 <- factor(segment_otu$group, level=levels(pcoa_otu_df$group1))

centroids_fun<- aggregate(cbind(pc1,pc2)~group1, pcoa_fun_df, mean) %>% filter(group1 %in%  uniques$group)
pcoa_ellipse_fun<- pcoa_fun_df %>% filter(group1 %in%  uniques$group)
ellipse_fun<- do.call(rbind, lapply(unique(pcoa_ellipse_fun$group1), function(t)
  data.frame(group = factor(as.character(t)),
             ellipse(cov(pcoa_ellipse_fun[pcoa_ellipse_fun$group==t,1:2]),
                     centre=as.matrix(centroids_fun[t,2:3]),
                     level=conf), stringsAsFactors=FALSE)))
ellipse_fun$group1 <- factor(ellipse_fun$group, 
                             level=unique(meta_info$group)[c(2,1,3,4)], ordered = TRUE)

#anosim test among different ocean zones
anosim_bc<- anosim(otu_bc, ordered(meta_info$group, level=unique(meta_info$group)))
anosim_wu<- anosim(otu_wu, ordered(meta_info$group, level=unique(meta_info$group)))
anosim_fun<- anosim(fun_bc, ordered(meta_info$group, level=unique(meta_info$group)))

#clustering samples by using weighted unifrac matrixes
wu_dist<- as.dist(otu_wu)
hc1<- hclust(wu_dist,"ave")
clus1<- cutree(hc1, 4)
clus2<- data.frame(label = names(clus1), member = factor(clus1))
design<- merge(clus2, meta_info, by.x = "label", by.y = "Run", all = F)
row.names(design)<- design$label

# annotating taxa information to different zOTUs
order_genus_m2_3<- aggregate(OTU ~ p1, data = diff_zotu_m2_3, FUN = length)
order_genus_m2_3<- order_genus_m2_3[order(order_genus_m2_3$OTU, decreasing = F),]
row.names(order_genus_m2_3)<- order_genus_m2_3$p1
order_genus_m2_3<- order_genus_m2_3[c(order_genus_m2_3$p1[!order_genus_m2_3$p1 %in% c("Unclassied_Bacteria","Unclassied_Archaea","Unclassified", "No_significent")],
                                      c("Unclassied","Unclassied_Archaea","Unclassified", "No_significent")),]

#scale the abundance of top12 KEGG pathways by using Z-Score
path_df<- scale(top20_path[1:12,-1:-4])
row.names(path_df)<- top20_path$level_3[1:12]


# set colors
col_group<- c("#0fb83d", "#c17848", "#0509ef", "#000000") 
names(col_group)<- unique(adiv_df$group)
col_pro<- c(col_group, OTU="#1AA7DA", KO="#CF992B")
col_taxa<-  c('blue', 'orange', 'green', 'yellow', 'red', 'hotpink', 'cyan',"aquamarine1",
              'purple', 'burlywood1', 'skyblue', "darkgreen","magenta","azure", 'gray')
names(col_taxa)<- rev(levels(phylum_top10$Taxonomy))
col_fun<- colorRampPalette(c("navy", "white", "firebrick3"))(100)
col_diff<- c(colorRampPalette(brewer.pal(12,"Paired"))(length(unique(diff_zotu_m2_3$p))),
             "grey")
names(col_diff)<- c(unique(diff_zotu_m2_3$p), "No_significent")

# set groups for comparison test
div_comparisons <- list(c("MG1","MG2"), c("MG2", "MG3"),
                        c("MG3", "MG4"), c("MG2", "MG4"))

#ploting by using ggplot2

div_stat11<- ggbetweenstats(data = subset(adiv_df,div==levels(adiv_df$div)[1]),
                            x = group1, y = index, type = "nonparametric", plot.type = 'violin',
                            ylab = levels(adiv_df$div)[1]) +
  scale_color_manual(values = col_group) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.ticks=element_line(color="black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.title.y = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title=element_blank())

div_stat12<- ggbetweenstats(data = subset(adiv_df,div==levels(adiv_df$div)[2]),
                            x = group1, y = index, type = "nonparametric", plot.type = 'violin',
                            ylab = levels(adiv_df$div)[2]) +
  scale_color_manual(values = col_group) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.ticks=element_line(color="black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.title.y = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title=element_blank())

div_stat13<- ggbetweenstats(data = subset(adiv_df,div==levels(adiv_df$div)[3]),
                            x = group1, y = index, type = "nonparametric", plot.type = 'violin',
                            ylab = levels(adiv_df$div)[3]) +
  scale_color_manual(values = col_group) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.ticks=element_line(color="black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.title.y = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title=element_blank())

div_stat14<- ggbetweenstats(data = subset(adiv_df,div==levels(adiv_df$div)[4]),
                            x = group1, y = index, type = "nonparametric", plot.type = 'violin',
                            ylab = levels(adiv_df$div)[4]) +
  scale_color_manual(values = col_group) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.ticks=element_line(color="black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.title.y = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title=element_blank())

div_stat1<- patchwork::wrap_plots(div_stat11, div_stat12, div_stat13, div_stat14,
                                  nrow = 2)

p_pcoa_otu<- pcoa_otu_df %>% ggplot() +
  geom_polygon(data = ellipse_otu, aes(x=pc1, y=pc2, colour=group1, fill=group1), 
               colour = "black", linetype=3, size =0.2, alpha=0.1, show.legend=F) + 
  geom_point(aes(x=pc1, y=pc2, color=group), size = 3.4, alpha = 0.58)+
  scale_color_manual(values = col_group) +
  scale_fill_manual(values = col_group) +
  xlab(paste("PCoa1 [", pcoa_otu_eig[1]," %]", sep="")) +
  ylab(paste("PCoa2 [", pcoa_otu_eig[2]," %]", sep="")) +
  ggtitle(label= paste0("ANOSIM:\nR2 = ", round(anosim_bc$statistic, 3),"\np-value = ",format.pval(anosim_bc$signif))) +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.ticks=element_line(color="black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title=element_blank())

p_pcoa_fun<- pcoa_fun_df %>% ggplot() +
  geom_polygon(data = ellipse_fun, aes(x=pc1, y=pc2, colour=group1, fill=group1), 
               colour = "black", linetype=3, size =0.2, alpha=0.1, show.legend=F) + 
  geom_point(aes(x=pc1, y=pc2, color=group), size = 3.4, shape = 17, alpha = 0.58)+
  scale_color_manual(values = col_group) +
  scale_fill_manual(values = col_group) +
  xlab(paste("PCoa1 [", pcoa_fun_eig[1]," %]", sep="")) +
  ylab(paste("PCoa2 [", pcoa_fun_eig[2]," %]", sep="")) +
  ggtitle(label= paste0("ANOSIM:\nR2 = ", round(anosim_wu$statistic, 3),"\np-value = ",format.pval(anosim_wu$signif))) +
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.ticks=element_line(color="black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.title=element_blank())

p_otu_fun_pcoa<- ggarrange(p_pcoa_otu, p_pcoa_fun, ncol = 2, legend = "bottom",
                           labels = c("B","C"),widths = c(3.4,3.4), common.legend = T)

p_tree_wu<- ggtree(hc1) %<+% design +
  geom_tippoint(shape=16, aes(colour = factor(group), x=x)) +
  scale_colour_manual(values = col_group) +
  theme_tree2(legend.position='bottom') +
  theme(legend.title=element_blank())

p_phylum_top10<- ggplot(data = phylum_top10,aes(x = samples1, y = 100 * value, alluvium = Taxonomy, stratum = Taxonomy)) +
  geom_alluvium(aes(fill = Taxonomy),alpha = .5, width = 0.6) +
  geom_stratum(aes(fill = Taxonomy, colour = Taxonomy), width = 0.6) +
  scale_y_continuous(expand=c(0, 0.8))+
  labs(y = 'Relative Abundance(%)') +
  scale_fill_manual(values =  col_taxa) +
  scale_colour_manual(values =  rev(col_taxa)) + 
  theme_classic() +
  theme(axis.line.x = element_line(colour = "black"),
        panel.background = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y= element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom") +
  coord_flip() +
  guides(fill=guide_legend(reverse = T),
         color=guide_legend(reverse = T))

p_tree_phylum<- p_phylum_top10 %>% aplot::insert_left(p_tree_wu,width = 0.4) 


p_path_heatmap<- Heatmap(path_df, name = "Zscore abundance",col = col_fun,
                         show_column_names = F,
                         cluster_columns = T,
                         column_title_gp = gpar(col=col_group[c(3,4,1,2)]),
                         row_title_gp = gpar(fontsize = 7.2),
                         row_names_gp = gpar(fontsize = 8),
                         row_title_rot = 0,
                         column_split = meta_info$group,
                         row_split = top20_path$X.level_1,
                         heatmap_legend_param = list(legend_direction = "horizontal",
                                                     title_position = "topcenter"))

pdf("fig3_kegg_top12.pdf", width = 8, height = 5)
draw(p_path_heatmap, heatmap_legend_side = "bottom")
dev.off()
1        

p_diff_m2_3<- ggplot(diff_zotu_m2_3, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color=p2, size=baseMean), alpha = 0.73, shape=19) +
  geom_hline(yintercept=-log10(0.01),
             linetype=4, 
             color = 'black', 
             size = 0.5) +
  geom_vline(xintercept=c(-2,2),
             linetype=4, 
             color = 'black', 
             size = 0.5) +
  geom_text_repel(data = diff_zotu_m2_3,
                  aes(x=log2FoldChange, y=-log10(padj), label=ID), 
                  colour = "black",
                  size=2.5,
                  hjust=0.5, 
                  vjust=0.5) +
  scale_color_manual(name=bquote(paste(FDR <= 0.01 , " and " ,  '|log2FoldChange| >= 2')),
                     values = col_diff[unique(diff_zotu_m2_3$p1)]) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.ticks=element_line(color="black"),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black")) +
  labs(x= bquote(log[2]~FoldChange), 
       y= bquote(-log[10]~padj)) +
  ggtitle(paste("MG2", "MG3", sep = " vs. "))

p_diff_legend<- ggplot(diff_zotu_m2_3, 
                       aes(x = log2FoldChange, 
                           y = -log10(padj), color=p2, size=baseMean)) +
  geom_point(alpha=0.73, shape=19)+
  annotate('rect', xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill="white")+
  scale_color_manual(name=bquote(paste(FDR <= 0.01 , " and " ,  '|log2FoldChange| >= 2')),
                     values = col_diff) +
  theme_classic()  +
  theme(legend.position = c(0.45, 0.45), 
        axis.ticks=element_blank(),
        axis.line = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank()) +
  guides(size=guide_legend(nrow = 1))

p_diff_otu<- ggarrange(p_diff_m2_3, p_diff_legend, ncol = 2,
                       labels = c("A",""))

ggsave(filename = "fig2_div_zotu.pdf", div_stat1, width=12, height=9)
ggsave(filename = "fig2_pcoa_otu_fun.pdf", p_otu_fun_pcoa, width=8, height=4.5)
ggsave(filename = "fig2_phylum_top15.pdf", p_tree_phylum, width=14, height=4.5)
ggsave(filename = "fig3_diff_otu.pdf", p_diff_otu, width=8, height=4.5)


