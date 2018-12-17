# Z score Ranking

# Global
options(stringsAsFactors = FALSE)
rm(list=ls())

# Libraries

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

# Working Directory. 

setwd("/Users/hannahralph/Desktop/DPhil/Bioinformatics/Other_Scripts/smMIP_analysis/clonalHaem-ZscoreRanking")

# Load input files

long_pileups_subset = fread("full_table_joined_noSNPs.csv")

long_pileups_subset[is.na(long_pileups_subset)] <- 0

# Calculate z-scores

message(paste("Calculating z-scores"))

# Add z scores

long_pileups_subset = long_pileups_subset %>% 
  group_by(individual,NucleotideChange) %>%
  mutate(z = (AF - mean(AF)) / sd(AF)) %>%
  ungroup


# Calculate P-values

message(paste("Calculating p-values"))

results = long_pileups_subset %>%
  group_by(position, NucleotideChange) %>%
  mutate(
    mean_z_for_pos_change = mean(z),
    sd_z_for_pos_change = sd(z),
    p = pnorm(z, mean = mean_z_for_pos_change, sd = sd_z_for_pos_change, lower.tail = FALSE)
  ) %>%
  ungroup

# NAdjust P values

results$bonferroni = p.adjust(results$p, method= "bonferroni")

results$fdr = p.adjust(results$p, method= "fdr")

# Remove 0 AF

results_non0AF = filter(results, AF != 0)

# Annovar Output

write.csv(results, "results_padj_.csv", row.names = FALSE)

# Plot all Pvalues - histograms

hist(results$p)
hist(results_non0AF$p)

hist(results$bonferroni)
hist(results_non0AF$bonferroni)

hist(results$fdr)
hist(results_non0AF$fdr)

hist(results$AF)
hist(results_non0AF$AF)

# Plot all Pvalues - QQ_plot

p_values = results$p
p_values_sorted = sort(p_values)
p_values_sorted_expected = (1:(length(p_values_sorted))/length(p_values_sorted))
qqplot(p_values_sorted,p_values_sorted_expected)

p_values = results$bonferroni
p_values_sorted = sort(p_values)
p_values_sorted_expected = (1:(length(p_values_sorted))/length(p_values_sorted))
qqplot(p_values_sorted,p_values_sorted_expected)

p_values = results_non0AF$p
p_values_sorted = sort(p_values)
p_values_sorted_expected = (1:(length(p_values_sorted))/length(p_values_sorted))
qqplot(p_values_sorted_expected,p_values_sorted)

p_values = results$fdr
p_values_sorted = sort(p_values)
p_values_sorted_expected = (1:(length(p_values_sorted))/length(p_values_sorted))
qqplot(p_values_sorted,p_values_sorted_expected)

# Plot by position

results_to_plot = filter(results, position == "chr10:123276893" & NucleotideChange == "A>T")

mean_z_at_pos = results_to_plot$mean_z_for_pos_change[1]
SD_z_at_pos = results_to_plot$sd_z_for_pos_change[1]

hist(results_to_plot$z)
abline(v=mean_z_at_pos,col="blue")
abline(v=(mean_z_at_pos+SD_z_at_pos),col="blue")
abline(v=(mean_z_at_pos-SD_z_at_pos),col="blue")



