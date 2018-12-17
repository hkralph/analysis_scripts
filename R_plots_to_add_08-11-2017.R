# Make RD plots to add to Nils protein diagrams

# Global
options(stringsAsFactors = FALSE)
rm(list=ls())

# library

library(ggplot2)
library(dplyr)
library(ggrepel)

# Working Directory. 


# Deva

#working_directory = "/Users/hannahralph/Desktop/DPhil/Experimental/RD Validation/RainDance_Scripts/PAE_hotspot_graphs/transcript_plots"

working_directory = "/t1-data/user/hralph/plots"

setwd(working_directory)

# Nils code. 

if (working_directory == "/t1-data/user/hralph/plots") {
  source('/t1-data/user/koelling/_share/data.R')
  source('/t1-data/user/koelling/_share/uniprot_plot.R')
}

# Tables

germline_mutations <- read.csv("selfish_mutation_codons_HR.csv")

RD_af_positive_subset <- read.csv("RD_AF_positive_subset_paper.csv")

RD_median_coverage_positive_subset <- read.csv("RD_median_coverage_positive_subset.csv")

# Get gene list

gene_list = unique(RD_af_positive_subset$Gene)
number_of_genes = length(gene_list)

message(paste("There are", number_of_genes, "genes to plot"))

# Make graphs

for (current_gene in 1:number_of_genes) {
  
  current_gene_name = gene_list[current_gene]
  # Testing
  #current_gene_name = "PTPN11"

    RD_cov_data = filter(RD_median_coverage_positive_subset, Gene == current_gene_name)
    RD_af_data = filter(RD_af_positive_subset, Gene == current_gene_name)
    
    current_transcript = unique(RD_cov_data$Transcript)
    
    message(paste("Processing", current_gene_name, "transcript", current_transcript))
    
    # Lists to hold graphs
    
    RD_variant_plots=list()
    amplicon_median_coverage_plots = list()
    germline_mutation_plots=list()
    
    # RD Variant Plot
    
    # Get protein annotations
    
    unique_protein_table = data.frame("Protein"=RD_af_data$Protein, "Codon"=RD_af_data$Codon)
    unique_protein_table = distinct(unique_protein_table)
    max_rate_vector <- c()
    
    for (current_unique_protein in 1:nrow(unique_protein_table)) {
      current_protein = unique_protein_table$Protein[current_unique_protein]
      # value_to_return = with(LookUpTable, ColumnToReturn[as.character(Columnto Search) == SearchKey])
      max_rate = with(RD_af_data, RD_af_data$AF[as.character(RD_af_data$Protein) == current_protein])
      max_rate = max(max_rate)
      max_rate_vector <- c(max_rate_vector, max_rate)
    }
    
    unique_protein_table$Max = max_rate_vector
    
    # Variants above 1% present in FGFR2, KRAS, NF1, PTPN11
    
    if (current_gene_name == "FGFR2" | current_gene_name == "KRAS" | current_gene_name == "NF1" | current_gene_name == "PTPN11") {
      y_limit = 3
    } else if (current_gene_name != "FGFR2") {
      y_limit = 1
    }
    
    RD_variant_plots[['Variants']] =
      
      ggplot() +
      geom_point(data=RD_af_data, 
                 aes(x=RD_af_data$Codon, 
                     y=RD_af_data$AF, 
                     colour=factor(RD_af_data$Indiv)),shape=4) +
      labs(x="Codon", # Or amplicon
           y="Observed (%)",
           title=current_gene_name) +
      ylim(0,y_limit) +
      guides(colour=guide_legend(title="Individual")) +
      geom_text_repel(data=unique_protein_table, 
                      aes(x=unique_protein_table$Codon,
                          y=unique_protein_table$Max), 
                      label=unique_protein_table$Protein,
                      size=3, nudge_y = 0.05, nudge_x=0.05)
    
    # Amplicon Coverage Plot
    
    amplicon_median_coverage_plots[['Coverage']] =
      
      ggplot() + 
      geom_segment(RD_cov_data, 
                   mapping=aes(x=StartCodon,
                               xend=EndCodon,
                               y=CovMedian,
                               yend=CovMedian)) +
     # scale_color_identity(name="Amp", guide="legend", breaks=c("black"), labels= c("Amplicon")) +
      scale_y_continuous(limits=c(0, 35000)) +
      labs(x="Codon",y="Median",title=current_gene_name) +
      geom_hline(aes(lty="amplicon", yintercept = 10000), colour="#990000", show.legend = TRUE) +
      scale_linetype_manual(name="legend", values=c("amplicon"="dashed"))
    
    # geom_errorbar(RD_cov_data, 
    #           mapping=aes(x=c(MiddleCodon), 
    #                       ymin=c(CovMin),
    #                       ymax=c(CovMax))) 
    
    # Get germline track
    
    if (current_gene_name == "FGFR2" | current_gene_name == "PTPN11") {
      
      current_germline_data = filter(germline_mutations, germline_mutations$Impacted.Gene == current_gene_name)
      current_label_data = filter(current_germline_data, current_germline_data$Reports >= 5)
      
      current_label_data$Reports <- NULL
      current_label_data$Reference <- NULL
      current_label_data$Syndrome <- NULL
      
      current_label_data = distinct(current_label_data)
      reports_df = data.frame()
      max_reports_df = data.frame()
      
      for (current_row in 1:nrow(current_label_data)) {
        current_protein_codon = current_label_data$ProteinCodon[current_row] 
        current_label_height = with(current_germline_data, Reports[as.character(ProteinCodon) == current_protein_codon])
        current_label_height = current_label_height[1]
        reports_df = rbind(reports_df, current_label_height)
      }
      
      current_label_data = bind_cols(current_label_data, reports_df)
      colnames(current_label_data) <- c("Impacted.Gene", "ProteinCodon", "CodonNumber", "Reports")
      
      # Only RD overlap
      
      RD_variants = data.frame(RD_af_data$Codon)
      RD_variants = distinct(RD_variants)
      colnames(RD_variants) = c("CodonNumber")
      
      joined_var = inner_join(RD_variants, current_germline_data, by="CodonNumber")
      
      germline_label_data = data.frame("Label"=joined_var$CodonLabel,
                                       "Codon"=joined_var$CodonNumber)
      germline_label_data = distinct(germline_label_data)
      
      for (current_row in 1:nrow(germline_label_data)) {
        current_protein_codon = germline_label_data$Codon[current_row] 
        max_reports = with(joined_var, Reports[as.character(CodonNumber) == current_protein_codon])
        max_reports = max(current_label_height)
        max_reports_df = rbind(max_reports_df, max_reports)
      }
      
      germline_labels  = bind_cols(germline_label_data,max_reports_df)
      colnames(germline_labels) <- c("Label", "Codon", "Max")
      
      summarised_germline_data = current_germline_data %>%
        group_by(CodonNumber, Syndrome) %>%
        summarise(syndromeReports = sum(Reports)) %>%
        group_by(CodonNumber) %>%
        summarise(sumReports = sum(syndromeReports),
                  commonSyndrome = Syndrome[which.max(syndromeReports)])
      
      germline_mutation_plots[['Germline']] =
        
        ggplot() +
        
        geom_bar(data=summarised_germline_data,
                 aes(x=CodonNumber,
                     y=sumReports,
                     fill=commonSyndrome), stat="identity", width=3) +
        scale_y_log10() +
        scale_fill_brewer(palette='Set2', type='qual')
      
      
    }
    
    
    # Add plots to Nils
    
    # add_line_for_cosmic
    # highlight_positions
    # title
    
    pdf_file_name = paste(working_directory, "/", current_gene_name, "_", Sys.Date(), "_protein_diagram.pdf", sep="")
    
    graph_title = paste(current_gene_name, ": ", sep="")
    
    if (working_directory == "/t1-data/user/hralph/plots" & current_gene_name == "FGFR2" | current_gene_name == "PTPN11") {
      
      pdf(pdf_file_name, width = 16, height = 8) #output filename
      print(plot_uniprot(current_transcript, show_tracks = c('protein','cosmic_all'), 
                         extra_tracks = c(RD_variant_plots[1], amplicon_median_coverage_plots[1], germline_mutation_plots[1]), 
                         extra_heights=c(3,1,2), 
                         use_new_cosmic = TRUE,
                         add_line_for_cosmic=FALSE,
                         title=graph_title
      ))
      dev.off()
    }
    
    if (working_directory == "/t1-data/user/hralph/plots" & current_gene_name != "FGFR2") {
      
      pdf(pdf_file_name, width = 16, height = 8) #output filename
      print(plot_uniprot(current_transcript, show_tracks = c('protein','cosmic_all', 'clinvar'), 
                         extra_tracks = c(RD_variant_plots[1], amplicon_median_coverage_plots[1]), 
                         extra_heights=c(3,1), 
                         cosmic_version = 82,
                         use_new_cosmic = TRUE,
                         add_line_for_cosmic=FALSE,
                         title=graph_title
      ))
      dev.off()
    }
    
    
    # if (working_directory == "/t1-data/user/hralph/plots") {
    # 
    #   pdf(pdf_file_name, width = 16, height = 8) #output filename
    #   ##print(plot_uniprot(current_transcript, show_tracks = c('exac_common', 'exac_rare', 'exac_doubletons', 'cosmic_all', 'protein'), extra_tracks = coverage_plots[current_gene], extra_heights=1, highlight_positions=252))
    #   #print(plot_uniprot(current_transcript, show_tracks = c('protein','cosmic_all', 'clinvar'), 
    #   print(plot_uniprot(current_transcript, show_tracks = c('protein','cosmic_all'), 
    #                      extra_tracks = c(RD_variant_plots[1], amplicon_median_coverage_plots[1], germline_mutation_plots[1]), 
    #                      #extra_tracks = c(germline_mutation_plots[1], amplicon_coverage_plots[1]),
    #                      extra_heights=c(3,1,2), 
    #                      #extra_heights=c(3,0.1), 
    #                      use_new_cosmic = TRUE,
    #                      add_line_for_cosmic=FALSE,
    #                      title=graph_title
    #                      ))
    #   dev.off()
    
    
    
  } # end of for current gene
  
  
    
    
    
  
  