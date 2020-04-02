library(ggplot2)
library("scales")

# needed to invert y-axis while making it log scale
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}


ggplot(pathogens_combined, aes(x=global_rank, y=count)) + geom_jitter(width = 0.1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_x_continuous(trans=reverselog_trans(10))