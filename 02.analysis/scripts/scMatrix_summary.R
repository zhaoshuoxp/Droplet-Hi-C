#!/projects/ps-renlab/y2xie/anaconda3/envs/seurat/bin Rscript
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

OP2 <- function(x) {### sparse matrix handling
    nms <- colnames(x)
    uniquenms <- unique(nms)
    sparseMatrix(i = x@i + 1, j = match(nms, uniquenms)[x@j + 1], x = x@x, 
        dimnames = list(rownames(x), uniquenms), repr = "C")
}

args <- commandArgs(trailingOnly = TRUE)
file <- args[1] # plain text matrix folder
cnv <- args[2] # GAM residuals
outpath <- args[3] # output directory

### Generate sparse matrix
### read colnames / rownames
name <- read.table(paste0(file, "/barcode.tsv"), header = T)
name <- name %>% 
mutate(bin = paste0(chrom, ":", start, "-", end)) %>% 
tibble::rownames_to_column("index")

### read count
dataset <- read.table(paste0(file, "/matrix.mtx"), header = T)
dataset <- dataset %>% 
mutate(bin1_id = paste0(chrom1, ":", start1, "-", end1)) %>% 
mutate(bin2_id = paste0(chrom2, ":", start2, "-", end2))
dataset$bin1_index <- name[match(dataset$bin1_id, name$bin), "index"]
dataset$bin2_index <- name[match(dataset$bin2_id, name$bin), "index"]

### coarse to sparse
mat <- sparseMatrix(i = as.integer(dataset[,"bin1_index"]), 
                    j = as.integer(dataset[,"bin2_index"]), 
                    x = dataset[,"count"])

rownames(mat) <- name$bin[1:nrow(mat)]
colnames(mat) <- name$bin[1:ncol(mat)]

tmp <- mat
colnames(tmp) <- gsub("^([^:]+):([^:]+)$", "\\1", colnames(tmp))
tmp <- as(tmp, "dgTMatrix")
obj_mtx_collapse <- OP2(tmp)

### summarize contact matrix
mat_contacts <- colSums(mat) %>% as.data.frame %>% setNames("contacts")
mat_contacts$contacts_bin <- colSums(mat != 0)

### intra
intra <- list()
mat_idx <- match(gsub("^([^:]+):([^:]+)$", "\\1", rownames(obj_mtx_collapse)), colnames(obj_mtx_collapse))

for (i in 1:length(mat_idx)){
    intra[[i]] <- obj_mtx_collapse[i, mat_idx[i]]
}

mat_contacts$intra <- unlist(intra)

### inter
mat_contacts$inter <- rowSums(obj_mtx_collapse) - mat_contacts$intra

#### get bin (binarize first)
tmp <- mat
tmp@x[tmp@x > 0] <- 1
colnames(tmp) <- gsub("^([^:]+):([^:]+)$", "\\1", colnames(tmp))
tmp <- as(tmp, "dgTMatrix")
obj_mtx_collapse <- OP2(tmp)
intra <- list()
mat_idx <- match(gsub("^([^:]+):([^:]+)$", "\\1", rownames(obj_mtx_collapse)), colnames(obj_mtx_collapse))

for (i in 1:length(mat_idx)){
    intra[[i]] <- obj_mtx_collapse[i, mat_idx[i]]
}

mat_contacts$intra_bin <- unlist(intra)

mat_contacts$inter_bin <- rowSums(obj_mtx_collapse) - mat_contacts$intra_bin

### add GAM residual
cov <- read.table(paste0(cnv), header = F) %>% setNames(c("chr", "start", "end", "residuals"))
cov <- cov %>% 
mutate(bin_id = paste0(chr, ":", start, "-", end)) %>%
tibble::column_to_rownames("bin_id")

mat_contacts <- merge(mat_contacts, cov[, "residuals", drop = F], by = 0) %>% tibble::column_to_rownames("Row.names")

### export
write.table(mat_contacts, paste0(outpath, "/contacts_summary.txt"), row.names = T, col.names = T, sep = "\t", quote = F)
