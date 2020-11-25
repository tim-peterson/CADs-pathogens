library(ggplot2)
library("scales")


# Figure 4C
ggplot(pathogens_ge0, aes(x=cnt_all, y=cnt_all_0)) + geom_jitter(width = 0.1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2')

# https://stackoverflow.com/questions/51159196/geom-jitter-colour-based-on-values
ggplot(pathogens_ge0, aes(x=cnt_all, y=cnt_all_0)) + geom_jitter(aes(colour = ifelse(cnt_all_0 > 1 | cnt_all_0 < -1, "#3953A4", "#231F20")), width = 0.1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position="none") + scale_x_continuous(trans='log2') + scale_color_identity()