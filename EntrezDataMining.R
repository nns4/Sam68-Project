
############################################# Command line Argument Script for Retrive Data from Entrez Databases ##########################################



## Installation of packages required for command line retrieve of data from entrez databases
#install.packages("rentrez")
#install.packages("tibble")
#install.packages("tidyr")
#install.packages('magrittr')



##################################### Run this on terminal ##########################################

########################### Rscript Entrez-retrive-data-Commandline.R dbs query ####################

# 1st argument is the database name e.g. gene (you can use entrez_dbs() function from rentrez package to see all database 
# 2nd argument is the query term protein name or gene name e.g sam68 
# e.g: Rscript Entrez-retrive-data-Commandline.R gds Sam68 can be typed in terminal to retrieve all Sam68 data from GEODatasets.
# Function used has a threshold around 370 so three different statements were used to retrieve data below 370 datasets, above 1000 and bewteen 370 and 1000 data.

######################################################################################################

library(tibble)
library(tidyr)
library(magrittr)
# help functions in rentrez package that will help user learn their way around the NCBI’s databases.
library(rentrez)


# entrez_dbs() function to find the list of available databases
entrez_dbs()


# Provides access to a copy of the command-line arguments given when this script is ran on terminal.  
args <- commandArgs(trailingOnly=TRUE)

# Retrieve information from entrez  databases with query (sam68 or KHDRBS1)
#The entrez database is argument 1 and argument 2 is the gene or protein name.
sam68_search <- entrez_search(db=args[1], term=args[2],retmax= 20,use_history = TRUE) 
# check the max count ids hit from the above line of code
sam68_search$count
length(sam68_search$ids)

sam68_search <- entrez_search(db=args[1], term=args[2],retmax= sam68_search$count,use_history = TRUE,api_key="3f39c1edf9f3a80911635e7336762e802f08")

#retrieve any datasets that is below 370 hits
if ( length(sam68_search$ids) < 370 ) {
  #summary_fourthquartile <- entrez_summary(db="gene",sam68_search$ids[((length(sam68_search$ids)*0.75)+1):(length(sam68_search$ids))])
  summary <- entrez_summary(db=args[1],sam68_search$ids)
  summary
  
  # converting list output(summary_second_half) into data frame
  all_cols <- tibble(all_col = summary)
  data<- all_cols %>% unnest_wider(all_col)
  # some columns are list format in dataframe by converting into character we can then write them into .csv or .tsv format
  data[] <- lapply(data, as.character)
  
  # write dataframe into .csv 
  write.csv(data,paste0(args[1],".csv"),row.names = FALSE)
  
  #this will get any data that is large and above 1000 hits
} else if ( length(sam68_search$ids) > 1000){
  #data is split into small chunks for easy retrieval.
  ids <- sam68_search$ids
  chunk_size <- 300
  ids_chunked <- split(ids, ceiling(seq_along(ids)/chunk_size))
  sam68_search <- entrez_search(db=args[1], term=args[2],retmax=600,use_history = TRUE)
  
  #df_total <- data.frame()
  for (i in 1:length(ids_chunked)) {
    y <- entrez_summary(db="gene",ids_chunked[[i]])
    all_cols <- tibble(all_col = y)
    data<- all_cols %>% unnest_wider(all_col)
    # some columns are in list format in dataframe by converting into character we can then write them into .csv or .tsv format
    data[] <- lapply(data, as.character)
    #df_total <- rbind(df_total,data)
    #print(paste0("nucc",i,".csv"))
    write.csv(data,paste0(args[1],i,".csv"))
  }
  #retrieve datasets that is between 370 and 1000.
} else {
  
  
  # provide the max count hints from entrez_search() function over a query
  sam68_search <- entrez_search(db=args[1], term=args[2],retmax=sam68_search$count,use_history = TRUE)
  
  
  # there is limit for entrez_summary() function so here divide the ids into four quartile for easier processing 
  #summary_firstquartile <- entrez_summary(db="gene",sam68_search$ids[1:(length(sam68_search$ids)*0.25)])
  summary_firstquartile <- entrez_summary(db=args[1],sam68_search$ids[1:(length(sam68_search$ids)*0.25)])
  summary_firstquartile
  
  # converting list output(summary_first_half) into data frame
  all_cols <- tibble(all_col = summary_firstquartile)
  data1<- all_cols %>% unnest_wider(all_col)
  # some columns are list type in dataframe by converting into character we can then write them into .csv or .tsv format
  data1[] <- lapply(data1, as.character)
  
  
  
  #summary_secondquartile <- entrez_summary(db="gene",sam68_search$ids[((length(sam68_search$ids)*0.25)+1):(length(sam68_search$ids)*0.5)])
  summary_secondquartile <- entrez_summary(db=args[1],sam68_search$ids[((length(sam68_search$ids)*0.25)+1):(length(sam68_search$ids)*0.5)])
  summary_secondquartile
  
  # converting list output(summary_second_half) into data frame
  all_cols <- tibble(all_col = summary_secondquartile)
  data2<- all_cols %>% unnest_wider(all_col)
  # some columns are list type in dataframe by converting into character we can then write them into .csv or .tsv format
  data2[] <- lapply(data2, as.character)
  
  #summary_thirdquartile <- entrez_summary(db="gene",sam68_search$ids[((length(sam68_search$ids)*0.5)+1):(length(sam68_search$ids)*0.75)])
  summary_thirdquartile <- entrez_summary(db=args[1],sam68_search$ids[((length(sam68_search$ids)*0.5)+1):(length(sam68_search$ids)*0.75)])
  summary_thirdquartile
  
  # converting list output(summary_second_half) into data frame
  all_cols <- tibble(all_col = summary_thirdquartile)
  data3<- all_cols %>% unnest_wider(all_col)
  # some columns are list type in dataframe by converting into character we can then write them into .csv or .tsv format
  data3[] <- lapply(data3, as.character)
  
  #summary_fourthquartile <- entrez_summary(db="gene",sam68_search$ids[((length(sam68_search$ids)*0.75)+1):(length(sam68_search$ids))])
  summary_fourthquartile <- entrez_summary(db=args[1],sam68_search$ids[((length(sam68_search$ids)*0.75)+1):(length(sam68_search$ids))])
  summary_fourthquartile
  
  # converting list output(summary_second_half) into data frame
  all_cols <- tibble(all_col = summary_fourthquartile)
  data4<- all_cols %>% unnest_wider(all_col)
  # some columns are list type in dataframe by converting into character we can then write them into .csv or .tsv format
  data4[] <- lapply(data4, as.character)
  
  data1 <- data1[,colnames(data2)]
  
  # combine both dataframe
  data <- rbind(data1,data2,data3,data4)
  
  # write dataframe into .csv 
  write.csv(data,paste0(args[1],".csv"),row.names = FALSE)
  
}


############################################ Script below can be used to read the columns from .csv file ####################################################
# After getting the datasets an example.csv file can be loaded into r environment using read.csv() function
#data <-read.csv("example.csv")
# check all column names in data variable
#colnames(data)
# choose columns you want to store in new dataframe
#new_data <- data[c("uid","taxid","project_id","project_acc")]
# Rewrite the new file again to save it.
#write.csv(new_data,"new_data.csv")

####################################################################################################################
