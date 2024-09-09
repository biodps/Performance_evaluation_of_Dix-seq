library("ggplot2")
library("hrbrthemes")
library("ggrepel")
library(reticulate)
library(ggsci)

df = reticulate::py$df
feature_tab_density_plt = df |> ggplot(aes(x=log10(cumulative_abundance), fill = software)) +
    geom_density(color="#e9ecef", alpha=0.7) +
    scale_x_reverse(expand = c(0.2, 0.2)) +
    # geom_vline(xintercept = 1, linetype = 'solid', linewidth = 0.3) +
    # geom_vline(xintercept = 4, linetype = 'dashed', linewidth = 0.3) +
    # geom_vline(xintercept = 8, linetype = 'dotted', linewidth = 0.3) +
    # scale_linetype_manual(name = 'Minimum Abundance Threshould', labels = c('dada2_default_setting', 'unoise3_sensitive_setting', 'unoise3_default_setting')) +
    # scale_fill_npg() +
    xlab("Feature Count") +
    ylab("") + 
    theme_ipsum_tw() +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank() 
    )

feature_tab_density_plt
