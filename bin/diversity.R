#!/usr/bin/env Rscript
library(tidyverse)
library(vegan)
library(ggplot2)
library(optparse)
library(cowplot)

boxplots <- function(df_div, level,dict) {
  groups <- sapply(colnames(df_div), function(x) dict[[x]])
  df_t <- as.data.frame(t(as.matrix(df_div)))
  colnames(df_t) <- df_div[, 1]
  shannon_diversity <- round(diversity(df_t, index = "shannon"),2)
  simpson_diversity <- round(diversity(df_t, index = "simpson"),2)
  chao1_diversity <- round((estimateR(df_t))["S.chao1",],2)

  df_index <- (data.frame(colnames(df_div), groups, shannon_diversity,simpson_diversity,chao1_diversity )) %>% 
    rename(
      sample = colnames.df_div.
      )

  counts <- table(df_index$group)
  valid_groups <- names(counts[counts >= 3])
  df_to_plot <- df_index[df_index$group %in% valid_groups, ]
  output_file <- paste(level, "_diversity_index.csv",sep = "")

  write.csv(df_to_plot, output_file,row.names=FALSE) 

  ### boxplot
  shannon_plot <- ggplot(df_to_plot, aes(x=groups, y=shannon_diversity, color=groups)) + 
      geom_boxplot(alpha=0.2) +
      geom_jitter(color="#999999", size=0.4, alpha=0.8) +
      theme(legend.position="none",,panel.background = element_blank(),
      axis.text.x = element_text(angle = 65, hjust = 1)) +
      scale_fill_brewer(palette="BuPu") +
      xlab("group") + ylab("Shannon Index")  #+
      #geom_signif(test="wilcox.test", comparisons = list(c("HC", "INR")), map_signif_level=TRUE) 

  simpson_plot <- ggplot(df_to_plot, aes(x=groups, y=simpson_diversity, color=groups),col = "white") + 
      geom_boxplot(alpha=0.2) +
      geom_jitter(color="#999999", size=0.4, alpha=0.8) +
      theme(legend.position="none",panel.background = element_blank(),
      axis.text.x = element_text(angle = 65, hjust = 1)) +
      xlab("group") + ylab("Simpson Index") 
      
  chao1_plot <- ggplot(df_to_plot, aes(x=groups, y=chao1_diversity, color=groups),col = "white") + 
      geom_boxplot(alpha=0.2) +
      geom_jitter(color="#999999", size=0.4, alpha=0.8) +
      theme(legend.position="none",panel.background = element_blank(),
      axis.text.x = element_text(angle = 65, hjust = 1)) +
      xlab("group") + ylab("Chao1 Index") 

  output_file <- paste(level, "_diversity_boxplot", ".pdf", sep = "")
  pdf(width=10, height=4,file=output_file)
  g = plot_grid(nrow=1,shannon_plot, simpson_plot, chao1_plot, labels = "AUTO")
  print(g)
  dev.off()
}

args <- commandArgs(trailingOnly=TRUE)
groups <- args[2]
groups <- groups[1]
groups<- sub("\\[", "", groups)
groups <-sub("\\]", "", groups)
groups <-sub(":::", "", groups)
groups <-sub("::", "", groups)
groups <-sub(":$", "", groups)

pairs <- strsplit(groups, ", ")[[1]]
key_value_list <- strsplit(pairs, ":")
group_dict <- list()
subgroup_dict <- list()
subsubgroup_dict <- list()

for (pair in key_value_list) {
  value <- pair[1]
  key <- pair[2]
  group_dict[[value]] <- key
  if(length(pair) == 4){
  subgroup_dict[[pair[1]]] <- pair[3]
  subsubgroup_dict[[pair[1]]] <- pair[4]
  }else if(length(pair) == 3){
  subgroup_dict[[pair[1]]] <- pair[3]
  }
}

df <- read.csv(args[1],check.names=FALSE)
samples <- names(df)
df$species <- NULL

boxplots(df,"groups",group_dict)
if(length(subgroup_dict)>0){
  boxplots(df,"subgroups",subgroup_dict)
}
if(length(subsubgroup_dict)>0){
  boxplots(df,"subsubgroups",subsubgroup_dict)

}