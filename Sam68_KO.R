#Reference material
#"https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html"
#BiocManager::install("limma")
library(limma)
library(edgeR)
count_data <- read.csv("GSE104856_counts.csv", row.names = 1)
head(count_data)

# Creating DGElist object from count matrix
DGE_object<-DGEList(counts = count_data)
DGE_object<-calcNormFactors(DGE_object)
dim(DGE_object)

# Filter low-expressed genes from count matrix 
cutoff <- 1
drop <- which(apply(cpm(DGE_object), 1, max) < cutoff)
DGE_norm <- DGE_object[-drop,] 
dim(DGE_norm) # number of genes left
DGE_norm$counts

# getting sample name into variable 
snames <- colnames(count_data)


group<-as.factor(substr(snames,start = 12,stop = 13))

# Multidimensional scaling (MDS) is a means of visualizing the level of similarity of individual cases of a dataset.
plotMDS(DGE_norm, col = as.numeric(group))


# design a matrix 
#we can think of ~0+x or equally ~x+0 as an equation of the form: y=ax+b. 
#By adding 0 we are forcing b to be zero, that means that we are looking for a line passing the origin (no intercept).
#If we indicated a model like ~x+1 or just ~x, there fitted equation could possibily contain a non-zero term b.
#Equally we may restrict b by a formula ~x-1 or ~-1+x that both mean: no intercept 
#(the same way we exclude a row or column in R by negative index).

model_matrix <- model.matrix(~0 + group)
model <- voom(DGE_norm, model_matrix, plot = T)

#What is voom doing?

#Counts are transformed to log2 counts per million reads (CPM),
#where oper million reads is defined based on the normalization factors we calculated earlier
#A linear model is fitted to the log2 CPM for each gene, and the residuals are calculated
#A smoothed curve is fitted to the sqrt(residual standard deviation) by average expression (see red line in plot above)
#The smoothed curve is used to obtain weights for each gene and sample that are passed into limma along with the log2 CPMs.

# fitting model 
#lmFit fits a linear model using weighted least squares for each gene:
fit <- lmFit(model, model_matrix)
head(coef(fit))

# create groups to comapre for Differential expression
contr <- makeContrasts(KO_Co=groupKO - groupCo, levels = colnames(coef(fit)))
contr

tmp <- contrasts.fit(fit, contr)

tmp <- eBayes(tmp)

# Summarize results
results <- decideTests(tmp)
summary(results)


# sort the gene according to p-values
top.table <- topTable(tmp, sort.by = "P")
head(top.table, 20)

hist(top.table[,"P.Value"])

length(which(top.table$adj.P.Val < 0.05))
write.csv(top.table, file = "Sam68-KO_new.csv")
