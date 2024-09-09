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
    pt = ggplot(dfs[[name]], aes(x = depth, y = pct, color = software, fill = software)) +
        
        #    scale_y_continuous(limits =c(0, 0.2), breaks = c(0, 0.04, 0.08, 0.12, 0.16, 0.20)) +
        theme_bw() +
        theme(aspect.ratio = 1/2.40) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    pt
    
    pp = pt + geom_bar(stat = "identity", position = "dodge", width = 0.5) +
        scale_color_npg() +
        scale_fill_npg() +
        scale_y_continuous(expand = expansion(c(0, 0.2)))
    
    pp
    
    ggsave(paste0('plots/', 'uclassified_features_presence_absence_pct_', name, ".pdf"), width = 8, height = 6)
    # ggsave(paste0('plots/', 'uclassified_features_total_reads_pct_', name, ".pdf"), width = 8, height = 6)
}
