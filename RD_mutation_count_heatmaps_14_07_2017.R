# This script is to produce heatmaps of the original RainDance data

# Reset R

rm(list=ls())

# set working directory

setwd("/Users/hannahralph/Desktop/DPhil/Experimental/RD Validation/RainDance_Scripts/mutation_count_heatmap")

#Load libraries

library(tidyr)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(pheatmap)

#

values_to_plot = read.csv("RD_mutation_counts_14-07-2017.csv")

# Read csv file of testes tables b = testes table. DIV/0 and N/A are problem-

b = read.csv("testestables.csv", colClasses =c("character", "integer", "character", "integer", "character", "character", "integer", "integer"), na.strings = c("#DIV/0!", "#N/A"))

# Join Testes Tables of column and row info with table containing variant info

file_name="RainDance_MutationCounts_PerPiece"

my_idv=c(3,4,5,10,14)
my_slice=c("A", "B", "C", "D", "EE", "F", "G", "J")



colours_to_use = colorRampPalette( c('#8E8E8C', 'peachpuff', 'gold2', 'firebrick4'))(5)

scale_to_use = seq(0, 4, length=6) # Needs to be one more than the number of colours I want
# Maximum mutations = 4


legend_breaks_to_use = seq(0, 4, length=5)

legend_labels_to_use = c("0","1", "2", "3", "4")


  joined_table = left_join(values_to_plot,b,by = "Ind.Slide.sample")
  
  for(current_individual in my_idv) {
    fn = paste("output_plots/Testis", current_individual, "plots.pdf", sep="_")
    pdf(fn)
    for(current_slice in my_slice) {
      filtered = joined_table %>% filter(individual==current_individual & Slice==current_slice)
      if(nrow(filtered) > 0) {
        message('Plotting mutation counts for individual: ', current_individual, ' slice: ', current_slice)
        
        # value.var = Column to
            filtered %>% acast(Row ~ Column, value.var = "Count") %>%
              pheatmap(color=colours_to_use,
                       breaks = scale_to_use,
                       cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE,
                       legend_breaks = legend_breaks_to_use, legend_labels = legend_labels_to_use,
                       border_color = FALSE, main = paste('Testis', current_individual, 'Slice', current_slice))
            
           }
        }

  dev.off()
  message(fn)
  }