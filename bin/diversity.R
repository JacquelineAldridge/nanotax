#!/usr/bin/env Rscript
library(tidyverse)
library(vegan)
library(ggplot2)
library(optparse)
library(cowplot)

args <- commandArgs(trailingOnly=TRUE)
groups <- args[2]
groups <- groups[1]
groups<- sub("\\[", "", groups)
groups <-sub("\\]", "", groups)

pairs <- strsplit(groups, ", ")[[1]]

key_value_list <- strsplit(pairs, ":")

group_dict <- list()

for (pair in key_value_list) {
  value <- pair[1]
  key <- pair[2]
  group_dict[[value]] <- key
}


df <- read.csv(args[1],check.names=FALSE)
# names(df) <- gsub("_T1$", "", names(df))
samples <- names(df)
df$species <- NULL

groups <- sapply(colnames(df), function(x) group_dict[[x]])
df_t <- as.data.frame(t(as.matrix(df)))
colnames(df_t) <- df[, 1]
shannon_diversity <- round(diversity(df_t, index = "shannon"),2)
simpson_diversity <- round(diversity(df_t, index = "simpson"),2)
chao1_diversity <- round((estimateR(df_t))["S.chao1",],2)

df_index <- (data.frame(colnames(df), groups, shannon_diversity,simpson_diversity,chao1_diversity )) %>%  #groups
  rename(
    sample = colnames.df.
    )

print(df_index)

counts <- table(df_index$group)

valid_groups <- names(counts[counts >= 3])

filtered_df_index <- df_index[df_index$group %in% valid_groups, ]

write.csv(filtered_df_index, "diversity_index.csv",row.names=FALSE) 





### boxplot
shannon_plot <- ggplot(filtered_df_index, aes(x=groups, y=shannon_diversity, color=groups)) + 
    geom_boxplot(alpha=0.2) +
    geom_jitter(color="#999999", size=0.4, alpha=0.8) +
    theme(legend.position="none",,panel.background = element_blank()) +
    scale_fill_brewer(palette="BuPu") +
    xlab("group") + ylab("Shannon Index")  #+
    #geom_signif(test="wilcox.test", comparisons = list(c("HC", "INR")), map_signif_level=TRUE) 

simpson_plot <- ggplot(filtered_df_index, aes(x=groups, y=simpson_diversity, color=groups),col = "white") + 
    geom_boxplot(alpha=0.2) +
    geom_jitter(color="#999999", size=0.4, alpha=0.8) +
    theme(legend.position="none",panel.background = element_blank()) +
    xlab("group") + ylab("Simpson Index") 
    
chao1_plot <- ggplot(filtered_df_index, aes(x=groups, y=chao1_diversity, color=groups),col = "white") + 
    geom_boxplot(alpha=0.2) +
    geom_jitter(color="#999999", size=0.4, alpha=0.8) +
    theme(legend.position="none",panel.background = element_blank()) +
    xlab("group") + ylab("Chao1 Index") 
pdf(width=10, height=4,file=paste("diversity_boxplot",".pdf", sep=""))
plot_grid(nrow=1,shannon_plot, simpson_plot, chao1_plot, labels = "AUTO")
dev.off()
