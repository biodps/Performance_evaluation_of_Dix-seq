pacman::p_load(ggplot2, hrbrthemes, ggrepel, ggsci, cowplot, tidyverse, patchwork, reticulate)

dfs <- py$dfs

for(name in names(dfs)) {
    df_up = dfs[[name]] %>% filter(indicator == 'cr_q1_97_otus_count_pct')
    # seqs count considered
    df_down = dfs[[name]] %>% filter(indicator != 'cr_q1_97_otus_count_pct')
    pt_up = ggplot(df_up, aes(x = software, y = value, color = software, fill = software))
    pt_down = ggplot(df_down, aes(x = software, y = value, color = software, fill = software))
    pp_up = pt_up + geom_segment(aes(x = software,
                                     xend = software, 
                                     y = 0, 
                                     yend = value)) +
        geom_point(aes(x = software, y = value, ), size = 4) +
        scale_color_npg() +
        scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25), expand = expansion(c(0, 0.1))) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_text(size = 7),
              aspect.ratio = 1/2) +
        ylab("Pecentage of features mapped to q1_97 OTUs")
    ggsave(paste0('plots/', 'asvs_mapped_to_q1_otus_', name, ".pdf"), pp_up, width = 8, height = 4)
    
    pp_down = pt_down + geom_segment(aes(x = software,
                                         xend = software, 
                                         y = 0, 
                                         yend = -value)) +
        geom_point(aes(x = software, y = -value, ), size = 4) +
        scale_color_npg() +
        #    scale_x_discrete(limits = rev) +
        scale_y_continuous(limits = c(-105, 0), breaks = seq(-100, 0, 25), expand = expansion(c(0.1, 0))) +
        theme_bw() +
        theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.y = element_text(size = 7),
              aspect.ratio = 1/2) +
        ylab("Pecentage of reads mapped to q1_97 OTUs")
    ggsave(paste0('plots/', 'reads_mapped_to_q1_otus_', name, ".pdf"), pp_down, width = 8, height = 4)
}

