#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
    stop("Unable to determine the script directory.")
}
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."))
source(file.path(repo_dir, "02.analysis", "scripts", "phc_help.R"))
source(file.path(repo_dir, "02.analysis", "scripts", "DPT_help.R"))
source(file.path(repo_dir, "02.analysis", "scripts", "basics.R"))

options <- commandArgs(trailingOnly = TRUE)
if (length(options) < 1 || !nzchar(options[1])) {
    stop("Sample name is required.")
}
sample <- options[1]
genome <- if (length(options) >= 2 && nzchar(options[2])) options[2] else "mm10"
if (!genome %in% c("mm10", "hg38")) {
    stop(paste0("Unsupported genome: ", genome))
}

stat_file <- file.path("mapping", paste0(sample, "_", genome, ".PairCount.stat.csv"))
    ### read count summary per cells, calculated during pre-processing
stat <- read.csv(stat_file, sep = "\t", row.names = 1) %>%
    setNames(c("total","mapped","unmapped","duplicate","cis","cis_1kb-","cis_1kb+","cis_10kb+","trans"))

valid <- PHCrankPair(obj = stat, prefix = sample)
write.table(valid, file = paste0(sample,"_cutoff.cells.txt"), row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

sum_file <- file.path("mapping", paste0(sample, "_", genome, ".sc.pairdedup.summary.txt"))
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
ggplot(data_filtered, aes(x = Group, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "",
       x = "",
       y = "Read Count") +
  theme_minimal()+
	theme(
		axis.text.y = element_text(size=14,colour = 'black'),
		axis.title.y = element_text(size=14,colour = 'black'),
		axis.text.x = element_text(angle=45, hjust=1,size=14,colour = 'black'),
        legend.text = element_text(size=14),
		legend.title = element_blank(),
	)->p

png(file=paste0(sample,'_range_qc.png'),height = 6, width = 6, res=600, units = "in", family="Arial")
p
dev.off()
