library(tidyverse)
library(gghalves)
library(ggsignif)
library(reticulate)
library(ggsci)

#读入数据
# dfs <- py$dfs
df = py$tax_diff_df_sw_rank
# l6最高趋势对的, ez库用错了所以不正常
# df = df %>% filter(rank %in% c("l6"))

#设置颜色
# my_color = 

# 使用for循环遍历命名列表

pt = ggplot(df, aes(x = depth, y = tax_diff_pct, color = rank, fill = rank)) +
    
    #    scale_y_continuous(limits =c(0, 0.2), breaks = c(0, 0.04, 0.08, 0.12, 0.16, 0.20)) +
    facet_wrap(~software, scales = "free") +
    theme_bw() +
    theme(aspect.ratio = 1/1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

pt

pp = pt + geom_step() +
    scale_color_npg() +
    scale_fill_npg() +
    scale_y_continuous(expand = expansion(c(0.05, 0.05)), limits = c(0, 100)) +
    scale_x_reverse(breaks = c(100, 98, 90, 82, 74, 66, 58, 50, 42, 34, 26, 18, 10))

pp

ggsave(paste0('plots/', 'tax_diff_roc_like_plt', ".pdf"), width = 16, height = 12)
