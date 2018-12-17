# This script is to produce heatmaps of the original RainDance data

# Reset R

rm(list=ls())

# set working directory

setwd("~/Desktop/DPhil/Experimental/RD Validation/08. RainDance Plots")

#check working directory

getwd()

#Load libraries

library(tidyr)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(pheatmap)

# Read csv file of testes tables and variants. a = variants b = testes table. DIV/0 and N/A are problem-

a = read.csv("fgfr2_s252w.csv", colClasses = c("character", "character", "character", "integer", "double", "double", "double", "double", "integer", "integer"), na.strings = c("#DIV/0!", "#N/A"))

b = read.csv("testestables.csv", colClasses =c("character", "integer", "character", "integer", "character", "character", "integer", "integer"), na.strings = c("#DIV/0!", "#N/A"))

Gene_Variant_Name = "FGFR2 S252W"

# Check the structure of files a and b to make sure they are read okay

str(a)
str(b)

# Join Testes Tables of column and row info with table containing variant info

joined_table = left_join(a,b,by = "Ind.Slide.sample")

my_idv=c(3,4,5,10,14)
my_slice=c("A", "B", "C", "D", "E", "F", "G", "J")


#filtered = joined_table %>% filter(individual==my_idv & Slice=="D")

# Do it in a loop


for(current_individual in my_idv) {
    pdf(paste("Testis", current_individual, Gene_Variant_Name, "plots.pdf", sep="_"))
    for(current_slice in my_slice) {
      filtered = joined_table %>% filter(individual==current_individual & Slice==current_slice)
      if(nrow(filtered) > 0) {
        message('Plotting data for individual: ', current_individual, ' slice: ', current_slice)
        
        filtered %>% acast(Row.x ~ Column.x, value.var = "A") %>%
          pheatmap(color=colorRampPalette((brewer.pal(n=7,name="Reds")))(100),
                   cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
                   main = paste(Gene_Variant_Name, 'Testis', current_individual, 'Slice', current_slice, 'Nucleotide:', "A % "))
        
        filtered %>% acast(Row.x ~ Column.x, value.var = "C") %>%
          pheatmap(color=colorRampPalette((brewer.pal(n=7,name="Reds")))(100),
                   cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
                   main = paste(Gene_Variant_Name, 'Testis', current_individual, 'Slice', current_slice, 'Nucleotide:', "C % "))
        
        filtered %>% acast(Row.x ~ Column.x, value.var = "G") %>%
          pheatmap(color=colorRampPalette((brewer.pal(n=7,name="Reds")))(100),
                   cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
                   main = paste(Gene_Variant_Name, 'Testis', current_individual, 'Slice', current_slice, 'Nucleotide:', "G % "))
        
        filtered %>% acast(Row.x ~ Column.x, value.var = "T") %>%
          pheatmap(color=colorRampPalette((brewer.pal(n=7,name="Reds")))(100),
                   cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
                   main = paste(Gene_Variant_Name, 'Testis', current_individual, 'Slice', current_slice, 'Nucleotide:', "T % "))
        
      } else {
        message('No data for individual: ', current_individual, ' slice: ', current_slice)
      }
    }
  dev.off()
}