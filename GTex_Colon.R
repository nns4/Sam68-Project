
#https://rstudio-pubs-static.s3.amazonaws.com/240657_5157ff98e8204c358b2118fa69162e18.html
#install.packages("WGCNA")
#BiocManager::install("impute")
library(WGCNA)

#BiocManager::install("CePa")
library(CePa)

# read the file in GTex_data it takes time wait for a while
GTex_data <- read.gct("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct")

options(scipen = 999)
GTex_data <- data.frame(GTex_data)

# finding the position of KHDRBS gene
grep("ENSG00000121774",rownames(GTex_data))
# GeneID for KHDRBS gene
GTex_data[grep("ENSG00000121774",rownames(GTex_data)),]


# transpose the data (columns to rows and vice-versa)
GTex_transpose <- transposeBigData(GTex_data)


# add KHDRBS1 gene column to starting of data
GTex_newSub = subset(GTex_transpose, select = -c(`ENSG00000121774.17`) )
ENSG00000121774.17<-GTex_transpose[,"ENSG00000121774.17"]
GTex <- cbind(ENSG00000121774.17,GTex_newSub)
# convert data to numeric type 
GTex <- sapply( GTex, as.numeric )

# Test for Association/Correlation Between Paired Samples (cor.test)
#  create empty dataframe (df)
df <- data.frame(Gene=colnames(GTex), Cor="", P.value="")
estimates = numeric(ncol(GTex))
pvalues = numeric(ncol(GTex))
for (i in 1:ncol(GTex)){
  test <- cor.test(GTex[,1], GTex[,i])
  estimates[i] = test$estimate
  pvalues[i] = test$p.value
}

df$Cor <- estimates
df$P.value <- pvalues
# this will contain correlation and p-value for each gene with KHDRBS1
colnames(df) <- c("Gene","Cor","pvalue_GTex")
df$Gene<-sub("\\.\\d+$", "", df$Gene)

# plot 
hist(df[,"Cor"],col="#7FFFD4",main = "Correlation Plot",xlab = "Cor")
cor.test()

# no. of genes with p-value less than 0.05
dfp <- df[df$pvalue_GTex < 0.05,]
dfp <- na.omit(dfp) ## count for  no. of genes with p-value less than 0.05 

# no. of genes with correaltion greater than 0.5 with sam68
dfcoeg <- dfp[dfp$Cor > 0.5,]
# no. of genes with correaltion less than 0.5 with sam68
dfcoel <- dfp[dfp$Cor < 0.5,]
