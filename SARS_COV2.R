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


ggplot(all_drugs_human_DEGs, aes(x=pval_x, y=odds_x, colour = ifelse(is.na(pka_), "#3953A4", "#231F20"))) + geom_jitter(width = 0.1) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + scale_x_continuous(trans=reverselog_trans(10))

df <- all_drugs_human_DEGs

df_nonCADs <- df[is.na(df[, "pka_"]), ]

df_CADs <- df[!is.na(df[, "pka_"]), ]

sars_gt7_gt0 <- read_delim("OneDrive-v3/Data/CADs/v2/PMID_32416070_32747830_human_DESeq2_Enrichr_DSigDB_intersect_pka_gt7_clogP_gt0.csv", 
                           +     "\t", escape_double = FALSE, trim_ws = TRUE)

df_nonCADs_gt70 <- sars_gt7_gt0[is.na(sars_gt7_gt0[, "pka_"]), ]

df_CADs_gt70 <- sars_gt7_gt0[!is.na(sars_gt7_gt0[, "pka_"]), ]


library(reshape2)
df_sars <- melt(data.frame(nonCADs = df_nonCADs['odds_e'], CADs = df_CADs['odds_e']))


write.csv(df_CADs['odds_e'], '/Users/timpeterson/OneDrive-v3/Data/CADs/v2/DrugSigDB_CADs_GE_for_graphing.csv', row.names = FALSE)

write.csv(df_nonCADs['odds_e'], '/Users/timpeterson/OneDrive-v3/Data/CADs/v2/DrugSigDB_nonCADs_GE_for_graphing.csv', row.names = FALSE)


write.csv(df_CADs_gt70['odds_e'], '/Users/timpeterson/OneDrive-v3/Data/CADs/v2/DrugSigDB_CADs_GE_for_graphing_gt_7_0.csv', row.names = FALSE)

write.csv(df_nonCADs_gt70['odds_e'], '/Users/timpeterson/OneDrive-v3/Data/CADs/v2/DrugSigDB_nonCADs_GE_for_graphing_gt_7_0.csv', row.names = FALSE)



ggplot(data = dat, aes (x = variable, y = value, fill = as.factor(iter))) + geom_violin(position = "dodge")