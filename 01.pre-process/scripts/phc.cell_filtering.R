#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
source("/nfs/baldar/quanyiz/app/Droplet-Hi-C/02.analysis/scripts/phc_help.R")
source("/nfs/baldar/quanyiz/app/Droplet-Hi-C/02.analysis/scripts/DPT_help.R")
source("/nfs/baldar/quanyiz/app/Droplet-Hi-C/02.analysis/scripts/basics.R")

options<-commandArgs(trailingOnly = T)

stat_file<-paste0('mapping',options[1],'_mm10.PairCount.stat.csv')
    ### read count summary per cells, calculated during pre-processing
stat <- read.csv(stat_file, sep = "\t", row.names = 1) %>% 
    setNames(c("total","mapped","unmapped","duplicate","cis","cis_1kb-","cis_1kb+","cis_10kb+","trans"))

valid <- PHCrankPair(obj = stat, prefix = options[1])
write.table(valid, file = paste0(options[1],"_cutoff.cells.txt"), row.names = F, col.names = F, sep = "\t", quote = F)

sum_file<-stat_file<-paste0('mapping',options[1],'_mm10.sc.pairdedup.summary.tx')
data <- read.table(sum_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(data) <- c("Category", "Count", "Fraction")

data$Group <- case_when(
  data$Category %in% c("Mapped_Read_Pairs", "Unmapped") ~ "Sequencing",
  data$Category %in% c("Intra_chromosomal", "Inter_chromosomal") ~ "Chromosomal",
  data$Category %in% c("short_range_1k", "long_range_1k") ~ "Range",
  TRUE ~ NA_character_ 
)

data_filtered <- data %>% filter(!is.na(Group))
data_filtered$Group <- factor(data_filtered$Group, levels = c("Sequencing", "Chromosomal", "Range"))
options(repr.plot.width=8, repr.plot.height=4)
ggplot(data_filtered, aes(x = Group, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "",
       x = "",
       y = "Read Count") +
  theme_minimal()+
	theme(
		axis.text.y = element_text(size=14,colour = 'black'),
		axis.title.y = element_text(size=14,colour = 'black'),
		axis.text.x = element_text(size=14,colour = 'black'),
        legend.text = element_text(size=14), 
		legend.title = element_blank(), 
	)->p

png(file=paste0(options[1],'_range_qc.png',),height = 4, width = 4, res=600, units = "in", family="Arial")
p
dev.off()