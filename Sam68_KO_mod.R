#"http://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html"
library(limma)
library(edgeR)
count_data <- read.csv("GSE104856_counts.csv", row.names = 1)
head(count_data)
# getting sample name into variable 
snames <- colnames(count_data)
group<-as.factor(substr(snames,start = 12,stop = 13))

dgeFull <- DGEList(count_data, group=group)
dgeFull


pseudoCounts <- log2(dgeFull$counts+1)
head(pseudoCounts)

boxplot(pseudoCounts, col="gray", las=2)
plotMDS(pseudoCounts)

par(mfrow=c(1,2))
## WT1 vs WT2
# A values
avalues <- (pseudoCounts[,1] + pseudoCounts[,2])/2
# M values
mvalues <- (pseudoCounts[,1] - pseudoCounts[,2])
plot(avalues, mvalues, xlab="A", ylab="M", pch=19, main="Control")
abline(h=0, col="red")
## Mt1 vs Mt2
# A values
avalues <- (pseudoCounts[,3] + pseudoCounts[,4])/2
# M values
mvalues <- (pseudoCounts[,3] - pseudoCounts[,4])
plot(avalues, mvalues, xlab="A", ylab="M", pch=19, main="KO")
abline(h=0, col="red")



# Differential expression analysis

dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)
head(dgeFull$counts)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)
dgeFull
dgeTest <- exactTest(dgeFull)
Sam68_KO<-dgeTest$table
Sam68_KO["Gene"] <- rownames(Sam68_KO)
colnames(Sam68_KO) <- c("logFC","logCPM","Pvalue_Sam_KO","Gene")
rownames(Sam68_KO) <- NULL




# Get Top results

resultTOP <- topTags(dgeTest, n=nrow(dgeTest$table))
head(resultTOP)
sum(resultTOP$table$FDR < 0.01)
sigReg <- resultTOP$table[resultTOP$table$FDR<0.01,]
sigReg <- sigReg[order(sigReg$logFC,decreasing = TRUE),]


