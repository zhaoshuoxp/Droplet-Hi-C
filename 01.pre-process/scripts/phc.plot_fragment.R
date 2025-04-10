#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
library(scales)
library(grid)

args <- commandArgs(trailingOnly = TRUE)
pairs <- args[1] # .pairs or .pairs to plot
pairs_name <- basename(pairs)
if (substr(pairs, nchar(pairs)-1, nchar(pairs)) == "gz") {
	# skip <- system(paste0("zcat ", pairs, " | grep '^#' - | wc -l"), intern = T)
	system(paste0("zcat ", pairs, " | awk '{if(substr($1,1,1)==\"#\"){next}{print}}' - | awk 'BEGIN {srand()} !/^$/ {if (rand() <= .001) print $0}' > ", pairs, ".tmp"))
} else {
	# skip <- system(paste0("grep '^#' ", pairs, " | wc -l"), intern = T)
	system(paste0("cat ", pairs, " | awk '{if(substr($1,1,1)==\"#\"){next}{print}}' - | awk 'BEGIN {srand()} !/^$/ {if (rand() <= .001) print $0}' > ", pairs, ".tmp"))
}
opairs <- fread(paste0(pairs, ".tmp"), sep="\t", header = F) #, skip = as.numeric(skip))
opairs <- opairs[opairs$V2 != "*" & opairs$V4 != "*", ]
total <- nrow(opairs)

if (dim(opairs)[2] < 10){
	warning("Barcodeds columns not detected in Hi-C pairs. Summarize fragments regardless of barcodes identity.\n", immediate. = T, call. = F)
} else {
	inter_cells <- length(which(opairs$V9 != opairs$V10))
	opairs <- opairs[which(opairs$V9 == opairs$V10),]
}
trans <- length(which(opairs$V2 != opairs$V4))
opairs <- opairs[which(opairs$V2 == opairs$V4),]
flen.median <- median(opairs$V5 - opairs$V3)
flen.mean <- mean(opairs$V5 - opairs$V3)
flen.25 <- quantile((opairs$V5 - opairs$V3), probs = 0.25)
flen.75 <- quantile((opairs$V5 - opairs$V3), probs = 0.75)
flen.max <- max(opairs$V5 - opairs$V3)
# my_anno = grid.text(paste0("sample:", pairs, "\n", "total pairs: ", total, "\n", "trans pair: ", trans, "\n", "meidna fragment size: ", flen.median, "\n", "mean fragment size: ", flen.mean), 
# 	x=Inf, y = Inf, vjust=1, hjust=1,
# 	gp=gpar(col="#c3533d", fontsize=10, fontface="bold"))

flen.hist <- as.data.frame(opairs$V5 - opairs$V3) %>% rename("flen" = 1) %>%
ggplot(aes(x=flen+1)) +
geom_histogram(bins=50, fill="#477a96", color="white") + 
xlab("Fragment length (bp)") + 
ylab("Fragment counts") + 
scale_x_log10(limits = c(1,10^8), labels = trans_format("log10", math_format(10^.x))) +
theme_bw() + 
annotate("text", label = paste0("Sample: ", pairs_name, "\n", "Median fragment size: ", as.integer(flen.median), "\n", "Mean fragment size: ", as.integer(flen.mean), "\n", "25% Quantile: ", flen.25, "\n", "75% Quantile: ", flen.75), 
	x = max(flen.max)/10, y = Inf, vjust = 1.5, hjust = 0)

ggsave(flen.hist, filename = paste0(pairs, ".fragmentLen.png"), dpi = 300, width = 16, height = 8)
# system(paste0("rm ", pairs, ".tmp"))
