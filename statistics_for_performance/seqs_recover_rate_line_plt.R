library(tidyverse)
library(gghalves)
library(ggsignif)
library(reticulate)
library(ggsci)

#读入数据
df <- py$df

#设置颜色
# my_color = 

pt = ggplot(df, aes(x = depth, y = valid_data_recover_rate, color = software, fill = software)) +
    
    scale_y_continuous(limits =c(0,100), breaks = c(0, 20, 40, 60 ,80 ,100)) +
    theme_bw() +
    theme(aspect.ratio = 1/2.40) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

pt

pp = pt + geom_boxplot(fill="white",
                           width=0.5, cex=0.8) +
    scale_color_npg() +
    scale_fill_npg() +
    geom_hline(yintercept = 80, linetype = "dashed", color = "black", size = 0.5)
    # geom_jitter( 
    #             size=1.5,
    #             position = position_jitter(width = 0.1))
    
pp

    # geom_half_violin(aes(fill=factor(Clade,levels = my_sort)),
    #                  side ='r',
    #                  position = position_nudge(x = .25, y = 0),
    #                  cex=0.8)+
    # geom_jitter(aes(color=factor(Clade,levels = my_sort)), 
    #             size=1.5,
    #             position = position_jitter(width = 0.1)) +
    # geom_signif(comparisons = list(c("subtills","cereus"),
    #                                c("cereus","megaterium"),
    #                                c("megaterium","circulans")), 
    #             map_signif_level = function(p) sprintf("p = %.2g", p), 
    #             step_increase = 0.05,
    #             y_position = c(23, 20, 17),
    #             textsize = 5, vjust = -0.2,
    #             tip_length = 0) +
    # scale_fill_manual(values = my_color)+
    # scale_color_manual(values = my_color,guide='none') +
    # theme_test(base_size =15)+
    # labs(x=NULL,y='No.of BGCs / genome')+
    # scale_y_continuous(limits =c(0,25))+
    # theme(legend.title =element_blank(),
    #       axis.title =element_text(size =21, face = "bold"),
    #       axis.text =element_text(color ='black',size = 17,face = "bold"),
    #       axis.text.x =element_text(face ='bold.italic'),
    #       legend.text = element_text(size = 14,face = "bold.italic"))

# ggsave("plots/box-violin-dot.png", width = 8, height = 6, dpi = 600)
ggsave("plots/seqs_recov_rate.pdf", width = 8, height = 6)
