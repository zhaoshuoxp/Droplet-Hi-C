#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(scales))

args <- commandArgs(trailingOnly = TRUE)
pairs <- args[1]
if (is.na(pairs) || !file.exists(pairs)) {
    stop("Error: Please provide a valid .pairs or .pairs.gz file.")
}
pairs_name <- basename(pairs)
input_cmd <- if (endsWith(pairs, ".gz")) {
    paste("gzip -dc --", shQuote(pairs))
} else {
    paste("cat --", shQuote(pairs))
}

# Peek at the first data line to check the number of columns
peek_cmd <- paste(input_cmd, "| awk '/^[^#]/ {print; exit}'")
first_line <- fread(cmd = peek_cmd, header = FALSE, sep = "\t")
if (nrow(first_line) == 0) {
    warning("No pair records found; fragment plot was not generated.", call. = FALSE)
    quit(save = "no", status = 0)
}

sel_cols <- c(2, 3, 4, 5)
if (ncol(first_line) >= 10) {
    sel_cols <- c(2, 3, 4, 5, 9, 10)
}

awk_cmd <- "awk 'BEGIN{srand()} /^#/{next} {seen++; if(seen == 1 || rand() <= 0.001) print}'"
read_cmd <- paste(input_cmd, "|", awk_cmd)

cat("Reading 0.1% sample of valid pairs directly into memory...\n")
opairs <- fread(cmd = read_cmd, header = FALSE, sep = "\t", select = sel_cols)

if (ncol(opairs) == 6) {
    setnames(opairs, c("chr1", "pos1", "chr2", "pos2", "cell1", "cell2"))
} else {
    setnames(opairs, c("chr1", "pos1", "chr2", "pos2"))
}

opairs <- opairs[chr1 != "*" & chr2 != "*" & chr1 != "!" & chr2 != "!"]

if ("cell1" %in% names(opairs)) {
    opairs <- opairs[cell1 == cell2]
} else {
    warning("Barcode columns not detected in Hi-C pairs. Summarize fragments regardless of barcode identity.\n", immediate. = TRUE, call. = FALSE)
}

opairs <- opairs[chr1 == chr2]

if (nrow(opairs) == 0) {
    warning("No same-cell cis pairs were sampled; fragment plot was not generated.", call. = FALSE)
    quit(save = "no", status = 0)
}

opairs[, flen := abs(pos2 - pos1)]

flen.median <- median(opairs$flen, na.rm = TRUE)
flen.mean   <- mean(opairs$flen, na.rm = TRUE)
flen.25     <- quantile(opairs$flen, probs = 0.25, na.rm = TRUE)
flen.75     <- quantile(opairs$flen, probs = 0.75, na.rm = TRUE)
flen.max    <- max(opairs$flen, na.rm = TRUE)
annotation.x <- max(flen.max / 10, 1)

cat("Generating plot...\n")
flen.hist <- ggplot(opairs, aes(x = flen + 1)) +
    geom_histogram(bins = 50, fill = "#477a96", color = "white") +
    xlab("Fragment length (bp)") +
    ylab("Fragment counts") +
    scale_x_log10(limits = c(1, 10^8), labels = trans_format("log10", math_format(10^.x))) +
    theme_bw() +
    annotate("text",
        label = paste0("Sample: ", pairs_name, "\n",
                       "Median fragment size: ", as.integer(flen.median), "\n",
                       "Mean fragment size: ", as.integer(flen.mean), "\n",
                       "25% Quantile: ", as.integer(flen.25), "\n",
                       "75% Quantile: ", as.integer(flen.75)),
        x = annotation.x, y = Inf, vjust = 1.5, hjust = 0)

out_png <- paste0(pairs, ".fragmentLen.png")
ggsave(flen.hist, filename = out_png, dpi = 300, width = 16, height = 8)
cat("Plot saved to", out_png, "\n")
