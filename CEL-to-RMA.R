#set working directory
setwd("/home/....")

#load required libraries
library(dplyr)
library(affy)
library(GEOquery)
library(tidyverse)

#set options to increase timeout and download method as wget to prevent download errors
options(timeout = max(300, getOption("timeout"),download.file.method.GEOquery = "wget"))

#Download raw file of a GSE study: eg. GSE201846 [loop can be created to download multiple GSEs]
getGEOSuppFiles("GSE201846")

#untar files
untar("GSE201846/GSE201846_RAW.tar", exdir = 'data/GSE201846/')

#reading in .cel files
raw.data <- ReadAffy(celfile.path = "data/GSE201846/")

#performing RMA normalization
normalized.data <- rma(raw.data)

#get expression estimates
normalized.expr <- as.data.frame(exprs(normalized.data))

#rename column names
names(normalized.expr) <- substring(names(normalized.expr),1,10)

#read matrix files to add clinical features and map probe IDs to gene symbols
gse <- getGEO("GSE201846", GSEMatrix = TRUE)

#fetch feature data to get ID - gene symbol mapping
feature.data <- gse$GSE201846_series_matrix.txt.gz@featureData@data

#subset required columns
feature..data <- feature.data[,c(1,11)]

#merge dataframe to get gene symbols
normalized.expr_final <- normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  inner_join(., feature..data, by = 'ID')

#rearrange dataframe
normalized.expr_final <- normalized.expr_final[,c(1,14,2:13)]

#fetch pheno data
pheno.data<- gse$GSE201846_series_matrix.txt.gz@phenoData@data

#subset required columns
pheno..data <- pheno.data[,c(2,6,8)]

#transpose dataframe
pheno..data<- as.data.frame(t(pheno..data))

#save pheno data
write.table(pheno..data,"data/pheno.txt",sep="",row.names=T, col.names =F)

#save expression data
write.table(normalized.expr_final,"data/norm.txt",sep=" ",row.names=F)

#read pheno and expression data
pheno <- read.table("data/pheno.txt")
norm <- read.table("data/norm.txt", fill=T)

#join both dataframes
GSE201846<- bind_rows(pheno, norm)

#save final merged matrix file
write.table(GSE201846,"data/Final_Matrix_GSE201846.txt",sep=" ",row.names=F, col.names = F)
