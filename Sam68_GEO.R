library(EnvStats)
library(dplyr)


all_RIP.KO_data <- merge(S_sam68_RIP,Sam68_KO,by = "Gene",all.x = TRUE,all.y = TRUE)

all_data <- merge(all_RIP.KO_data,df,by = "Gene",all.x = TRUE,all.y = TRUE)



all_sort_data <- select(all_data,c("Gene","Name","p_value_Sam68_RIP","pvalue_GTex","Pvalue_Sam_KO"))
write.csv(all_sort_data,"geo_pvalue.csv",row.names = FALSE)

geo_pvalue <- read.csv("geo_pvalue.csv")


geo_pvalue["geo_mean_pvalue"]<-apply(geo_pvalue[,3:5], 1, geoMean)
geo_pvalue["geo_mean_fdr"] <- p.adjust(geo_pvalue$geo_mean_pvalue)
write.csv(geo_pvalue, file = "Sam68_GeometricValues.csv")



