library(EnvStats)
library(dplyr)


# find ids from df (from GTex) data frame that match to geneID in Sam68_KO dataframe
x = c()
for (i in Sam68_KO$Gene) {
  idx=which(df$Gene == i)
  if (is.integer(idx)) {
    x =append(x,idx)
  }
  
}

# choose only those geneID which are in Sam68_KO
df <- df[x,]



# find ids from S_sam68_RIP data frame that match to geneID in df (from GTex) dataframe
y = c()
for (i in df$Gene) {
  idy=which(S_sam68_RIP$Gene == i)
  if (is.integer(idy)) {
    y =append(y,idy)
  }
  
}

# choose only those geneID which are in df (from GTex)
S_sam68_RIP <- S_sam68_RIP[y,]


# Due to less rows in df ( GTex data after subset) we have to subset Sam68_KO
z = c()
for (i in df$Gene) {
  idz=which(Sam68_KO$Gene == i)
  if (is.integer(idz)) {
    z =append(z,idz)
  }
  
}

# choose only those geneID which are in df (from GTex)
Sam68_KO <- Sam68_KO[z,]

### Now all 3 dataset have same rows with same geneID

# Generating a data frame with columns GENE, NAME , p_Value_RIP,p_Value_KO,p_Value_GTex
colnames(df)
sam68GTx <- select(df,c("pvalue_GTex"))
colnames(S_sam68_RIP)
sam68RIP <- select(S_sam68_RIP,c("Gene","Name","p_value_Sam68_RIP"))
colnames(Sam68_KO)
sam68KO <- select(Sam68_KO,c("Pvalue_Sam_KO"))

geometric_pvalue <- do.call("cbind", list(sam68RIP,sam68KO,sam68GTx))
write.csv(geometric_pvalue,"geo_pvalue.csv",row.names = FALSE)


# loading the pvalue dataset from 3 subset data of KO,RIP and GTex (run from here )
geo_pvalue <- read.csv("geo_pvalue.csv")

# geometric mean of pvalues from 3 sample

geo_pvalue["geo_mean_pvalue"]<-apply(geo_pvalue[,3:5], 1, geoMean)
geo_pvalue["geo_mean_fdr"] <- p.adjust(geo_pvalue$geo_mean_pvalue)

