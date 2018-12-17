# Options
options(stringsAsFactors = FALSE)
options(scipen = 999)
rm(list=ls())

# Variables to update

current_folder = "UMI-test"

# Set working directory

working_directory = "/Users/hannahralph/Desktop/DPhil/Bioinformatics/Other_Scripts/amplicon_pileup/"

setwd(paste(working_directory, current_folder, "output_tables", sep="/"))

#Libraries

library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)

# Input files

all_input_files <- Sys.glob("*formatted-pileup.txt")

for (current_input_file in 1:length(all_input_files)) {
  
  # Get variant table and info
  
  input_file = read.table(all_input_files[current_input_file], header=TRUE)
  current_file_name = gsub("_all-samples_formatted-pileup.txt", "", all_input_files[current_input_file])
  target_position_chr = input_file$chr[1]
  start_position = head(input_file$pos, 1)
  end_position = tail(input_file$pos, 1)
  
  # Melt table, format
  
  melted_table = input_file %>% melt(id.vars=c('Target', 'chr', 'pos', 'ref'))
  
  # Create new variables (sample and COV/A/C/G/T)
  
  melted_table$sample = sapply(strsplit(as.character(melted_table$variable), '_', fixed=TRUE), '[[', 1) 
  melted_table$sample = gsub("\\.", "-", melted_table$sample)
  
  melted_table$variable2 = sapply(strsplit(as.character(melted_table$variable), '_', fixed=TRUE), '[[', 3)
  melted_table$variable2 = gsub("\\.", "", melted_table$variable2)
  melted_table$variable2 = gsub("1", "", melted_table$variable2) # new
  
  melted_table$variable <- NULL
  
  colnames(melted_table) <- c("target", "chr", "pos", "ref_base", "value", "sample", "variable")
  
  melted_table$sample = factor(melted_table$sample)
  
  n_samples = length(levels(melted_table$sample))
  sample_shape_index = factor(rep(1:8, ceiling(n_samples / 8))[1:n_samples])
  melted_table[is.na(melted_table)] <- 0 # Force NA to 0 so that average can be calculated
  average_covs = melted_table %>% filter(variable == 'Coverage') %>% group_by(sample) %>% summarise(average=mean(value))
  
  # Graphs 
  
  graph_title = paste(current_folder, ":", current_file_name)
  
  file_name = paste("PileupPlot_", current_folder, "_", current_file_name, "_", Sys.Date(), ".pdf", sep="")
  
  pdf(paste(working_directory, current_folder, "/output_graphs/",file_name, sep=""), width=12)
  
  # Nucleotide Plot
  plot(ggplot(melted_table %>% filter(variable != 'Coverage' & variable != ref_base), 
              aes(pos, value, colour=variable, shape=ref_base)) 
       + xlab(paste("Position in Amplicon (", target_position_chr, ":", start_position, "-", end_position,")", sep="")) # X Label. Sep "" to remove spaces
       + ylab("Percentage non-reference (%)") # Y Label
       + ylim(0,0.6)
       + labs(title=graph_title) # Title
       + theme(legend.position = "top")
       + scale_colour_discrete(name="Non-Reference Nucleotide")
       + scale_shape_discrete(name="Reference Nucleotide")
       + geom_point())
  # Nucleotide Log Plot
  suppressWarnings(plot(ggplot(melted_table %>% filter(variable != 'Coverage' & variable != ref_base), 
                               aes(pos, value, colour=variable, shape=ref_base)) 
                        + xlab(paste("Position in Amplicon (chr", target_position_chr, ":", start_position, "-", end_position,")", sep="")) # X Label. Sep "" to remove spaces
                        + ylab("Percentage non-reference (%)") # Y Label
                        + labs(title=paste(graph_title, "(Log)")) # Title
                        + theme(legend.position = "top")
                        + scale_colour_discrete(name="Non-Reference Nucleotide")
                        + scale_shape_discrete(name="Reference Nucleotide")
                        + geom_point()
                        + scale_y_log10()))
  # Nucleotide Plot by Sample
  plot(ggplot(melted_table %>% filter(variable != 'Coverage' & variable != ref_base), 
              aes(pos, value, shape=sample, colour=sample)) 
       + xlab(paste("Position in Amplicon (chr", target_position_chr, ":", start_position, "-", end_position,")", sep="")) # X Label. Sep "" to remove spaces
       + ylab("Percentage non-reference (%)") # Y Label
       + labs(title=paste(graph_title, "(By Sample)")) # Title
       + labs(ref_base = "Reference Nucleotide")
       + geom_point()
       + scale_shape_manual(values=sample_shape_index[order(levels(melted_table$sample))])
       + theme(legend.position = "top"))
  # Coverage Plot
  plot(ggplot(average_covs, aes(sample, average)) 
       + xlab("Sample") # X Label
       + ylab("Average Coverage Across Amplicon") # Y Label
       + labs(title=paste(graph_title, "Coverage"))
       + theme(axis.text.x=element_text(angle=90))
       + geom_bar(stat='identity'))
  
  #Close file
  dev.off()
  
  message(paste("PDF created for ", current_file_name))
  
}
