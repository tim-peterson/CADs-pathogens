library(ggplot2)
library("scales")


# Figure 4C
ggplot(pathogens_ge0, aes(x=cnt_all, y=cnt_all_0)) + geom_jitter(width = 0.1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2')