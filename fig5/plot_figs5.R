################################################################################ 
# This script has been used to visualize the statistical results of Figures5
# in the manuscript "Dix-seq: An integrated pipeline for fast amplicon data analysis".  
# Pensgheng Dong @ Zhengzhou, China
# dpsh@henau.edu.cn 
################################################################################

#setting the pathway of the Working Directory
setwd("E:/YOUR_PATHWAY/fig5/")

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

get_pro_ell <- function(pro, grp_in, ellipse_pro = 0.97){
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
for (i in dir()[grep("rds$",dir())]) {
  t1<- gsub("\\.rds", "", i)
  assign(t1, readRDS(i))
}

# calculating Distance Matrix and making pcoa
for (i in ls()[grep("_bc$",ls())]) {
  tmp<- paste0(gsub("_bc", "", i), "_pcoa")
  assign(tmp, as.data.frame(pcoa(as.dist(get(i)))$vectors))
}

# set colors
col1<- c("#C80000","navy")
names(col1)<- unique(mapping$type)
col2<- c("#f0cba0", "#ff8733","#e38c48","#b35f24", "#582a0f","#C80000",
         "#c5f4fc", "#9dd8fb", "#6fa0f7", "#4653f4", "#351dcf", "navy")
names(col2)<- unique(gsub("[0-9].*$","",mapping$sample))
col3<- c(col1, col2, "grey65", "black")
names(col3)<- c(names(col1), names(col2), c("Tested method", "Wang et al"))

# procrustes test analysis and plot 

for (i in ls()[grep("_pcoa$",ls())][!ls()[grep("_pcoa$",ls())] %in% "q1s97_pcoa"]) {
  t1<- paste0(gsub("_pcoa", "", i), "_pro")
  tx<- paste0(gsub("_pcoa", "", i), "_prox")
  ty<- paste0(gsub("_pcoa", "", i), "_proy")
  tp<- paste0(gsub("_pcoa", "", i), "_plot")
  tl<- paste0(gsub("_pcoa", "", i), "_linetmp")
  td<- paste0(gsub("_pcoa", "", i), "_line")
  t2<- paste0(gsub("_pcoa", "", i), "_test")
  t3<- paste0(gsub("_pcoa", "", i), "_propvar")
  t4<- paste0(gsub("_pcoa", "", i), "_pval")
  t5<- paste0(gsub("_pcoa", "", i), "_m2val")
  t6<- paste0(gsub("_pcoa", "", i), "_ell")
  t7<- paste0(gsub("_pcoa", "", i), "_cor")
  t8<- paste0("p_", gsub("_pcoa", "", i))
  assign(t1, procrustes(get(i), q1s97_pcoa))
  assign(tx, data.frame(get(t1)$X)[,c(1,2)])
  assign(ty, data.frame(get(t1)$Yrot)[,c(1,2)])
  assign(tx, get(tx) %>% mutate(UserName=rownames(get(tx)), method="Tested method"))
  assign(ty, get(ty) %>% mutate(Axis.1=get(ty)[,1],
                                Axis.2=get(ty)[,2],
                                UserName=rownames(get(ty)), 
                                method="Wang et al"))
  assign(tp, bind_rows(get(tx), get(ty)) %>% select(-X1, -X2))
  assign(tp, get(tp) %>% mutate(group=factor(gsub("[0-9].*$","",mapping$sample[match(get(tp)$UserName, mapping$Run)]),
                                             levels = unique(gsub("[0-9].*$","",mapping$sample))),
                                type=mapping$type[match(get(tp)$UserName, mapping$Run)]))
  assign(tl, (get(tx)[rownames(get(tx)),c(1,2)] + get(ty)[rownames(get(ty)),c(1,2)])/2)
  assign(td, rbind(get(tp)[,1:2],get(tl),get(tl)) %>% mutate(method=c(get(tp)$method, rep(c("Tested method","Wang et al"),c(nrow(get(tx)),nrow(get(ty)))))))
  assign(td, get(td) %>% mutate(name=rep(paste(get(tp)$UserName,get(tp)$type,sep = "_"),2)))
  assign(t2, protest(get(i), q1s97_pcoa, perm = 999))
  assign(t3, signif(sqrt(get(t1)$svd$d)/sum(sqrt(get(t1)$svd$d)), 4)*100)
  assign(t4, signif(get(t2)$signif, 1))
  assign(t5, signif(get(t2)$ss, 2))
  assign(t6, get_pro_ell(get(t1), get(tp)$type, 0.95) %>% mutate(title=gsub("_pcoa", "", i)))
  assign(t7, signif(procrustes.syncsa(get(paste0(gsub("_pcoa", "", i), "_bc")), q1s97_bc), 3))
  assign(t8, get(t6) %>% ggplot(aes(x=x, y=y)) +
           geom_line(data = get(td), alpha=0.45, size = 0.58 ,
                     aes(x= Axis.1, y=Axis.2, group=name, color = method)) +
           geom_point(data = get(tp), size = 1.62, alpha=0.75, 
                      aes(x = Axis.1, y = Axis.2, color = group, shape = method)) +
           geom_polygon(aes_string(color = 'Groups', group = 'Groups'), lty=2, fill = NA) + 
           facet_grid( ~ title)+
           scale_color_manual(values = col3) +
           scale_shape_manual(values = c(1,16)) +
           theme_bw() +
           theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 legend.title = element_blank(),
                 panel.border = element_rect(colour = "black"),
                 legend.text = element_text(size=9),
                 legend.position = 'bottom',
                 axis.text = element_text(size=6.7,colour = "black"),
                 axis.title = element_text(size=9,colour = "black"),
                 plot.title = element_text(size = 7),
                 aspect.ratio = 1) +
           guides(color = guide_legend(nrow = 2), 
                  shape = guide_legend(nrow = 2)) +
           theme(strip.background = element_rect(color=NA, fill=NA)) +
           xlab(paste0("Procrustes1 [",get(t3)[1],"%]")) +
           ylab(paste0("Procrustes2 [",get(t3)[2],"%]")) +
           ggtitle(label= paste0("M2 value = ", get(t5),"\nCorrelation value = ",get(t7),"\np-value = ",get(t4)))
         
  )
  
}

# check the m12 sqaured value for different methods
for (i in ls()[grep("m2val", ls())]) { 
  tmp <- get(i); print(tmp)
}

# merge all plots of the different methods into total plot
ls()[grep("p_", ls())]
procruste_plot<- ggarrange(p_dixsq, p_easya, p_dada2, p_deblr, p_q1s99, NA,
                           ncol = 3, nrow = 2, align = "h", labels = LETTERS[1:5],
                           legend = "bottom", common.legend = T)

ggsave(filename = "Fig5_procruste_Bray-Curtis_plot.pdf", procruste_plot, width = 8, height = 8)


