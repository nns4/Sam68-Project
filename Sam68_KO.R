#Reference material

#"http://www.nathalievialaneix.eu/doc/html/solution_edgeR-tomato-withcode.html"
library(limma)
library(edgeR)
count_data <- read.csv("GSE104856_counts.csv", row.names = 1)
head(count_data)
# getting sample name into variable 
snames <- colnames(count_data)
group<-as.factor(substr(snames,start = 12,stop = 13))


      
y <- DGEList(count_data, group=group)
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

AveLogCPM <- aveLogCPM(y)
hist(AveLogCPM)
y <- calcNormFactors(y)
y$samples

pch <- c(6,7,9,11)
colors <- rep(c("darkgreen", "red"), 2)

# MDS plot
plotMDS(y, col=colors[group], pch=pch[group])
legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)

plotMD(y, column=1)
abline(h=0, col="red", lty=2, lwd=2)

y <- estimateDisp(y, design)

    
fit <- glmQLFit(y, design)
head(fit$coefficients)


KO.Co <- makeContrasts(KO-Co, levels=design)
res <- glmQLFTest(fit, contrast=KO.Co)

resultTOP <- topTags(res, n=nrow(res$table))
head(resultTOP)
  
is.de <- decideTestsDGE(res)
summary(is.de)

sigDownReg <- resultTOP$table[resultTOP$table$FDR<0.01,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
head(sigDownReg)
sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
head(sigUpReg)


plotMD(res, status=is.de, values=c(1,-1), col=c("red","blue"),
       legend="topright")
