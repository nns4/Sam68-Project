#https://replicationindex.com/2019/03/31/one-tail-or-two-tails-that-is-the-question/
## Calculating P-values
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
## Calculating P-values
options(scipen = 999)

# load the data in S_sam68 variable
S_sam68_RIP <- read.csv("S_sam68_RIP_WC.csv")

# calculating mean of 3 control columns for each gene 
S_sam68_RIP["control_mean"]<-rowMeans(S_sam68_RIP[,4:6])
# calculating sd of 3 control columns for each gene
S_sam68_RIP["control_sd"] <-apply(S_sam68_RIP[,4:6], 1, sd)

# calculating  Right tail p-values pnorm() function 
p_value <- pnorm(S_sam68_RIP$Sam68,mean = S_sam68_RIP$control_mean,sd = S_sam68_RIP$control_sd,lower.tail = FALSE)
S_sam68_RIP["p_value_Sam68_RIP"] <- p_value
S_sam68_RIP["fdr"] <- p.adjust(S_sam68_RIP$p_value)

# subset the gene name
S_sam68_RIP$Gene<-sub("\\.\\d+$", "", S_sam68_RIP$Gene)


# No. of significant gene with p-value less than 0.05
S_sam68_sig <- S_sam68_RIP[S_sam68_RIP$p_value_Sam68_RIP <0.05,]
