library(tidyverse)
library(gghalves)
library(ggsignif)
library(reticulate)
library(ggsci)

#读入数据
dfs <- py$dfs

#设置颜色
# my_color = 

# 使用for循环遍历命名列表
for(name in names(dfs)) {
    pt = ggplot(dfs[[name]], aes(x = software, y = value, color = software, fill = software)) +
        
        #    scale_y_continuous(limits =c(0, 0.2), breaks = c(0, 0.04, 0.08, 0.12, 0.16, 0.20)) +
        facet_wrap(~alpha_div, scales = "free") +
        theme_bw() +
        theme(aspect.ratio = 1/2.40) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    pt
    
    pp = pt + geom_half_violin(side ='r',
                               position = position_nudge(x = .25, y = 0),
                               cex=0.8) +
        geom_boxplot(fill=NA,
                     #                 color = "gray50",
                     width=0.1, cex=0.5,
                     position = position_nudge(x = .25, y = 0)
        ) +
        geom_jitter(aes(shape = software), size=1.5,
                    position = position_jitter(width = 0.1),
                    show.legend = T) +
        scale_color_npg() +
        scale_fill_npg() +
        scale_y_continuous(expand = expansion(c(0.4, 0.4))) +
        scale_shape_manual(values = c(1, 2, 3, 4, 5, 6))
    
    pp
    
    ggsave(paste0('plots/', 'aplha_div_', name, ".pdf"), width = 12, height = 9)
}

