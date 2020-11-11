#https://replicationindex.com/2019/03/31/one-tail-or-two-tails-that-is-the-question/
## Calculating P-values
options(scipen = 999)

# load the data in S_sam68 variable
S_sam68 <- read.csv("S_sam68_RIP_WC.csv")
names(S_sam68[,4:6])

# calculating mean of 3 control columns for each gene 
S_sam68["control_mean"]<-rowMeans(S_sam68[,4:6])
# calculating sd of 3 control columns for each gene
S_sam68["control_sd"] <-apply(S_sam68[,4:6], 1, sd)
S_sam68$control_sd<- round(S_sam68$control_sd,3)

# calculating  Right tail p-values pnorm() function 
p_value <- pnorm(S_sam68$Sam68,mean = S_sam68$control_mean,sd = S_sam68$control_sd,lower.tail = FALSE)
S_sam68["p_value"] <- p_value
S_sam68$p_value<- round(S_sam68$p_value,5)

write.csv(S_sam68, file = "Sam68_RIP_results.csv")
