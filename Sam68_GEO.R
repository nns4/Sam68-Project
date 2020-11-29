library(magrittr)
library(dplyr)
library(EnvStats)
options(scipen = 999)
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)

# reading in data from deseq2
geo_pvalue <- read.csv("geo_pvalue.csv")


geo_pvalue["geo_mean_pvalue"]<-apply(geo_pvalue[,3:5], 1, geoMean)
geo_pvalue["geo_mean_p.adjust"] <- p.adjust(geo_pvalue$geo_mean_pvalue)

# histogram plot 
hist(geo_pvalue$geo_mean_pvalue,col="red",main = "Hypergeometric Plot",xlab = "Pvalue")

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
allOE_genes <- as.character(geo_pvalue$Gene)

## Extract significant results
sigOE <- filter(geo_pvalue,geo_mean_p.adjust < 0.05)

sigOE_genes <- as.character(sigOE$Gene)

ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
#different plots for enrichment analysis
library(ggupset)
upsetplot(ego,title = "Upset Plot")
barplot(ego, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways")

dotplot(ego,showCategory = 10,title = "Dot Plot")
emapplot(ego,showCategory = 30,title = "Enrichment map")
goplot(ego,title = "Enriched GO induced graph")
ridgeplot(ego)
