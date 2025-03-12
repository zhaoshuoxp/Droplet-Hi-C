
## 2022-12-19
## Yang Xie (y2xie@ucsd.edu)
## Functions for Droplet Paired-Tag analysis in R

### set color ###
library(Matrix)
library(dichromat)
library(viridis)
source("/projects/ps-renlab2/y2xie/scripts/basics.R")

tori <- c("#477a96", "#c3533d", "#85212b", "#bcc9d1", "#6a9d84")
colfunc2 <- colorRampPalette(c(tori))

### mouse to human gene conversion or vice versa
convert_mouse_to_human <- function(gene_list, direction = "m2h"){
    mouse_human_genes = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt", sep = "\t")
    output = c()
    if(direction == "m2h"){
        qry <- "mouse, laboratory"
        ref <- "human"
    }else if(direction == "h2m"){
        ref <- "mouse, laboratory"
        qry <- "human"
    }else{
        stop("direction can only be m2h or h2m.")
    }
    for(gene in gene_list){
        class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name==qry))[['DB.Class.Key']]
        if(!identical(class_key, integer(0)) ){
            human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name==ref))[,"Symbol"]
            for(human_gene in human_genes){
                output = append(output,human_gene)
            }
        }
    }
    return(output)
}

### read in fragment counts to calculate frip ###
ImportArcFRiP <- function(raw_count, frip_count){
    tmp1 <- read.table(paste0(raw_count), header = T, row.names = 1)
    tmp2 <- read.table(paste0(frip_count), header = T, row.names = 1)
    frip <- merge(x = tmp1, y = tmp2, by = 0) %>% setNames(c("bc", "counts_total", "counts_unique", "counts_in_peaks", "unique_in_peaks"))
    frip$frip <- frip$unique_in_peaks / frip$counts_unique
    return(frip)
}

### Find FRiP cutoff ###
CutoffArcFRiP <- function(frip, xcut_low = 100, ycut_low = 0.05, ycut_high = 0.8, prefix){

}

### Plot calculated frip ###
PlotArcFRiP <- function(frip, xcut_low = 100, ycut_low = 0.05, ycut_high = 0.8, prefix){
    valid <- frip[frip$frip > ycut_low & frip$frip < ycut_high & frip$counts_unique > xcut_low, ]
    label <- data.frame(anno = paste0("Unique reads cutoff: ", xcut_low, "\n",
                                      "FRiP cutoff: ", ycut_low, "-", ycut_high, "\n", 
                                      "PF cells: ", nrow(valid), "\n")) 
                                      # "Median total reads: ", as.integer(median(valid$counts_total))))
    t1 <- frip %>%
    ggplot(aes(x = counts_unique, y = 100*frip)) + 
    geom_point(size = 0.5, alpha = 0.1, color = "grey") + 
    geom_vline(xintercept = xcut_low, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = 100*ycut_low, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = 100*ycut_high, color = colfunc2(1), linetype="dashed") +
    geom_point(data = valid, aes(x = counts_unique, y = 100*frip),
               size = 0.1, alpha = 0.5, color = colfunc2(3)[[2]]) + 
    geom_text(data = label, aes(x = max(frip$counts_unique), y = max(100*as.numeric(frip$frip)), 
                                hjust = 1, vjust = 1, label = anno), size = 2.5) +  
    theme_Publication() + xlab("log10(Fragments)") + ylab("% FRiP") + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +   
    scale_x_log10()
    
    write.table(valid, file = paste0(prefix, "_PF_cells.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
    ggsave(t1, filename = paste0(prefix, "_FRiP_macs2.png"), height = 6, width = 6, dpi = 300)
    return(valid)
}

### Plot calculated frip, color by local jitter density for better visualization ###
PlotArcFRiP_d <- function(frip, xcut_low = 100, ycut_low = 0.05, ycut_high = 0.8, prefix){
    filt_data <- frip
    ### assign color density 
    x <- densCols(log10(filt_data$counts_unique), 100*filt_data$frip, colramp=colorRampPalette(viridis(12)))
    frip$dens <- col2rgb(x)[1,] + 1L

    ### calculate density again to exclude low counts population
    filt_data <- filt_data %>% filter(counts_unique > 100)
    x <- densCols(log10(filt_data$counts_unique), 100*filt_data$frip, colramp=colorRampPalette(viridis(12)))
    frip[frip$counts_unique > 100, ]$dens <- (col2rgb(x)[1,] + 1L)

    valid <- frip[frip$frip > ycut_low & frip$frip < ycut_high & frip$counts_unique > xcut_low, ]
    label <- data.frame(anno = paste0("Unique reads cutoff: ", xcut_low, "\n",
                                      "FRiP cutoff: ", ycut_low, "-", ycut_high, "\n", 
                                      "PF cells: ", nrow(valid), "\n")) 
                                      # "Median total reads: ", as.integer(median(valid$counts_total))))
    t1 <- frip %>%
    ggplot(aes(x = counts_unique, y = 100*frip)) + 
    geom_point(size = 0.5, alpha = 0.1, aes(color = dens)) + 
    geom_vline(xintercept = xcut_low, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = 100*ycut_low, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = 100*ycut_high, color = colfunc2(1), linetype="dashed") +
    geom_text(data = label, aes(x = max(frip$counts_unique), y = max(100*as.numeric(frip$frip)), 
                                hjust = 1, vjust = 1, label = anno), size = 2.5) +  
    theme_Publication() + xlab("log10(Fragments)") + ylab("% FRiP") + 
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10)) +   
    scale_x_log10() + 
    theme(legend.position = "") + 
    scale_color_viridis(trans = "log")
    
    write.table(valid, file = paste0(prefix, "_PF_cells.txt"), row.names = F, col.names = F, sep = "\t", quote = F)
    ggsave(t1, filename = paste0(prefix, "_FRiP_macs2.png"), height = 6, width = 6, dpi = 300)
    return(valid)
}

### Plot CUT&TAG library fragment length pattern ###
PlotArcFragment <- function(frag_path, prefix){
    suppressMessages(library(data.table))
    frag <- fread(frag_path)
    frag$len <- frag$V3 - frag$V2
    frags <- frag[sample(nrow(frag), 10000, replace = F), ]
    frags <- frags[, c("len", "V5")]
    frags_pv <- frags[rep(1:nrow(frags), frags[["V5"]]), ]
    t6 <- frags_pv %>% dplyr::filter(len < 1000) %>%
    ggplot(aes(x = len)) + 
    geom_histogram(color = "white", bins = 100) +
    theme_Publication() + 
    scale_x_continuous(breaks = seq(0,1000,100)) + 
    xlab("Fragment length")
    ggsave(t6, filename = paste0(prefix, "_fragments.pdf"), height = 6, width = 6, dpi = 300)
    return(t6)
}

### After in-house processing of 10X Multiome, merge both modalities for plotting ###
PairArc <- function(dmat, rmat, names = "atac"){
    translate <- read.table("/projects/ps-renlab/y2xie/projects/genome_ref/arc_bc-translation.txt", header = T, row.names = 1, sep = "\t")
    translate$atac <- paste0(translate$atac, "-1")
    translate$rna <- paste0(translate$rna, "-1")
    # rownames(translate) <- paste0(translate$rna, "-1")
    if (names == "atac"){ ### cells name in dna
        idx <- match(colnames(rmat), translate$rna)
        colnames(rmat) <- translate[idx, names]
        colnames(dmat) <- colnames(dmat)
    }else if (names == "rna"){ ### cells name in rna
        idx <- match(colnames(dmat), translate$atac)
        colnames(dmat) <- translate[idx, names]
        colnames(rmat) <- colnames(rmat)
    }
    dna <- colSums(dmat) %>% as.data.frame() %>% setNames("dna") %>% tibble::rownames_to_column("bc")
    rna <- colSums(rmat) %>% as.data.frame() %>% setNames("rna") %>% tibble::rownames_to_column("bc")
    merge1 <- merge(x = dna, y = rna, by = "bc")
    return(merge1)
}

### Plot reads count for each pair of barcode. Save plot and segmented cells with name prefix ###
PlotArcPair <- function(pair, dcutoff = 100, rcutoff = 100, prefix, names = "atac"){
    # if(!identical(colnames(pair), c("bc", "dna_count", "rna_count"))){
    if((!all(c("bc", "dna_count", "rna_count") %in% colnames(pair)))){
        stop('colnames of input dataframe should be c("bc", "dna_count", "rna_count")')
    }
    translate <- read.table("/projects/ps-renlab/y2xie/projects/genome_ref/arc_bc-translation.txt", header = T, row.names = 1, sep = "\t")
    translate$atac <- paste0(translate$atac, "-1")
    translate$rna <- paste0(translate$rna, "-1")
    
    pair_pf <- pair[pair$dna_count > dcutoff & pair$rna_count > rcutoff, ]
    label <- data.frame(anno = paste0("DNA cutoff: ", min(pair_pf$dna_count), "\n", "DNA PF cells: ", nrow(pair[pair$dna_count > dcutoff, ]), "\n", 
                                      "RNA cutoff: ", min(pair_pf$rna_count), "\n",  "RNA PF cells: ", nrow(pair[pair$rna_count > rcutoff, ]), "\n", 
                                      "Both PF cells: ", nrow(pair_pf)))
    t1 <- pair %>% dplyr::filter(dna_count > 5 & rna_count > 5) %>%
    ggplot(aes(x = dna_count, y = rna_count)) + 
    geom_point(size = 0.1, color = "grey") + 
    geom_vline(xintercept = dcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = rcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_point(data = pair_pf, aes(x = dna_count, y = rna_count), size = 0.1, color = colfunc2(3)[[2]]) + 
    theme_Publication() + 
    geom_text(data = label, aes(x = max(pair$dna_count), y = min(pair$rna_count), hjust = 1, vjust = -0.5, label = anno), size = 2.5) + 
    ylab("Unique UMI (RNA)") + xlab("Unique fragemnts (DNA)") + 
    scale_x_log10() + 
    scale_y_log10()
    ggsave(t1, filename = paste0(prefix, "_valid_cells.png"), dpi = 300, height = 6, width = 6)
    
    ### export names
    idx <- match(pair_pf$bc, translate[, names])
    valid_bc <- translate[idx, c("atac", "rna")]
    valid_bc <- merge(x = valid_bc, y = pair_pf, by.x = names, by.y = "bc") ### atac, rna, bc, dna_count, rna_count
    write.table(valid_bc, file = paste0(prefix, "_valid_cells.xls"), row.names = F, col.names = T, sep = "\t", quote = F)
    return(valid_bc)
}

### Plot reads count for each pair of barcode, color by scatter density. Save plot and segmented cells with name prefix ###
PlotArcPair_d <- function(pair, dcutoff = 100, rcutoff = 100, prefix, names = "atac"){
    if((!all(c("bc", "dna_count", "rna_count") %in% colnames(pair)))){
        stop('colnames of input dataframe should be c("bc", "dna_count", "rna_count")')
    }
    translate <- read.table("/projects/ps-renlab/y2xie/projects/genome_ref/arc_bc-translation.txt", header = T, row.names = 1, sep = "\t")
    translate$atac <- paste0(translate$atac, "-1")
    translate$rna <- paste0(translate$rna, "-1")
    
    pair_pf <- pair[pair$dna_count > dcutoff & pair$rna_count > rcutoff, ]
    label <- data.frame(anno = paste0("DNA cutoff: ", min(pair_pf$dna_count), "\n", "DNA PF cells: ", nrow(pair[pair$dna_count > dcutoff, ]), "\n", 
                                      "RNA cutoff: ", min(pair_pf$rna_count), "\n",  "RNA PF cells: ", nrow(pair[pair$rna_count > rcutoff, ]), "\n", 
                                      "Both PF cells: ", nrow(pair_pf)))

    filt_data <- pair %>% dplyr::filter(dna_count > 30 & rna_count > 30)
    ### assign color density 
    x <- densCols(log10(filt_data$dna_count), log10(filt_data$rna_count), colramp=colorRampPalette(viridis(12)))
    pair$dens <- min(col2rgb(x)[1,] + 1L)
    pair[pair$dna_count > 30 & pair$rna_count > 30, ]$dens <- col2rgb(x)[1,] + 1L

    t1 <- pair %>% dplyr::filter(dna_count > 5 & rna_count > 5) %>%
    ggplot(aes(x = dna_count, y = rna_count)) + 
    geom_point(size = 0.1, aes(color = dens)) + 
    geom_vline(xintercept = dcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = rcutoff, color = colfunc2(1), linetype="dashed") + 
    theme_Publication() + 
    geom_text(data = label, aes(x = max(pair$dna_count), y = min(pair$rna_count), hjust = 1, vjust = -0.5, label = anno), size = 2.5) + 
    ylab("Unique UMI (RNA)") + xlab("Unique fragemnts (DNA)") + 
    scale_x_log10() + 
    scale_y_log10() + 
    theme(legend.position = "") + 
    scale_color_viridis(trans = "log")

    ggsave(t1, filename = paste0(prefix, "_valid_cells.png"), dpi = 300, height = 6, width = 6)
    
    ### export names
    idx <- match(pair_pf$bc, translate[, names])
    valid_bc <- translate[idx, c("atac", "rna")]
    valid_bc <- merge(x = valid_bc, y = pair_pf, by.x = names, by.y = "bc") ### atac, rna, bc, dna_count, rna_count
    write.table(valid_bc, file = paste0(prefix, "_valid_cells.xls"), row.names = F, col.names = T, sep = "\t", quote = F)
    return(valid_bc)
}

### Plot reads count for each pair of barcode by providing valid barcode ###
PlotArcPair2 <- function(pair, dvalid, rvalid, prefix, names = "atac"){
    if(!identical(colnames(pair), c("bc", "dna_count", "rna_count"))){
        stop('colnames of input dataframe should be c("bc", "dna_count", "rna_count")')
    }
    translate <- read.table("/projects/ps-renlab/y2xie/projects/genome_ref/arc_bc-translation.txt", header = T, row.names = 1, sep = "\t")
    translate$atac <- paste0(translate$atac, "-1")
    translate$rna <- paste0(translate$rna, "-1")
    pair$rna_bc <- translate[match(pair$bc, translate$atac), "rna"]
    validd <- which(pair$bc %in% dvalid)
    validr <- which(pair$rna_bc %in% rvalid)
    pair_pf <- pair[intersect(validd, validr), ]
    dcutoff <- min(pair_pf$dna_count)
    rcutoff <- min(pair_pf$rna_count)
    label <- data.frame(anno = paste0("DNA cutoff: ", dcutoff, "\n", "DNA PF cells: ", length(dvalid), "\n", 
                                      "RNA cutoff: ", rcutoff, "\n",  "RNA PF cells: ", length(rvalid), "\n", 
                                      "Both PF cells: ", nrow(pair_pf)))
    t1 <- pair %>% dplyr::filter(dna_count > 5 & rna_count > 5) %>%
    ggplot(aes(x = dna_count, y = rna_count)) + 
    geom_point(size = 0.1, color = "grey") + 
    geom_vline(xintercept = dcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_hline(yintercept = rcutoff, color = colfunc2(1), linetype="dashed") + 
    geom_point(data = pair_pf, aes(x = dna_count, y = rna_count), size = 0.1, color = colfunc2(3)[[2]]) + 
    theme_Publication() + 
    geom_text(data = label, aes(x = max(pair$dna_count), y = min(pair$rna_count), hjust = 1, vjust = -0.5, label = anno), size = 2.5) + 
    ylab("Unique UMI (RNA)") + xlab("Unique fragemnts (DNA)") + 
    scale_x_log10() + 
    scale_y_log10()
    ggsave(t1, filename = paste0(prefix, "_valid_cells.png"), dpi = 300, height = 6, width = 6)
    
    valid_bc <- pair_pf[, c("bc", "rna_bc", "dna_count", "rna_count")] %>% setNames(c("atac", "rna", "dna_count", "rna_count"))
    # write.table(valid_bc, file = paste0(prefix, "_valid_cells.xls"), row.names = F, col.names = T, sep = "\t", quote = F)
    return(valid_bc)
}



### Fast run Seurat ###
RunRNA <- function(obj, reduction = "pca", var = "none", batch.label = "none"){
    obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2500)
    if (var != "none") {
        obj <- ScaleData(obj, vars.to.regress = var)
    } else {
        obj <- ScaleData(obj)
    }
    obj <- RunPCA(obj, features = VariableFeatures(object = obj), verbose = F)
    if (batch.label != "none") {
        obj <- RunHarmony(obj, group.by.vars = batch.label)
        reduction <- "harmony"
    }
    obj <- RunUMAP(obj, dims = 1:30, min.dist = 0.1, seed.use=131, reduction = reduction,
                   n.components = 2L, umap.method = "uwot", n.neighbors = 25, uwot.sgd=TRUE, verbose=F)
    return(obj)
}

### Fast run Signac ###
RunDNA <- function(obj){
    obj <- RunTFIDF(obj, method = 3)
    obj <- FindTopFeatures(obj, min.cutoff = 'q15')
    obj <- RunSVD(obj)
    obj <- RunUMAP(object = obj, seed.use = 2022, reduction = 'lsi', dims = 2:30, 
                   min.dist = 0.01, n.neighbors = 15, verbose = F)
    # obj <- FindNeighbors(object = obj, reduction = 'lsi', dims = 2:30, verbose=F)
    # obj <- FindClusters(object = obj, resolution = 0.3)
    return(obj)
}

### sparse matrix: collapse columns with same colname, return sparse still ###
OP2 <- function (x) {
    nms <- colnames(x)
    uniquenms <- unique(nms)
    sparseMatrix(i = x@i + 1, j = match(nms, uniquenms)[x@j + 
        1], x = x@x, dimnames = list(rownames(x), uniquenms), 
	dims = c(nrow(x), length(uniquenms)),
        repr = "C")
}

### calculate cpm/rpkm by cluster (aggreagated cells for calculation) ###
    ### for histone, nCount_histone is used. for RNA, nCount_RNA is used. ###
ArcXPM <- function (obj_mtx, meta, method = "CPM", group.by, gname = "gene_id", gene_length = "/projects/ps-renlab/y2xie/projects/genome_ref/mm10-2020-A_build/mm10_gene_allowlist.bed") {
    colnames(obj_mtx) <- meta[colnames(obj_mtx), group.by]
    # obj_mtx <- as(obj_mtx, "dgTMatrix")
    obj_mtx <- as(obj_mtx, "TsparseMatrix")
    obj_mtx_collapse <- OP2(obj_mtx)
    spars <- length(obj_mtx_collapse@x)/obj_mtx_collapse@Dim[1]/obj_mtx_collapse@Dim[2]
    cat(paste0("sparsity: ", spars, "\n"))
    if (spars > 0.2) {
        obj_mtx_collapse <- as(obj_mtx_collapse, "matrix")
        cat("coarse dgTMatrix into Matrix.\n")
    }
    if (method == "CPM") {
        readSums <- aggregate(meta$nCount_histone, list(meta[, group.by]), sum)
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ", length(names(readSums)), "\n")
        cat("check obj_mtx_collapse: ", length(colnames(obj_mtx_collapse)), 
            "\n")
        readSums <- readSums[order(match(names(readSums), colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^6/readSums)
    }
    else if (method == "RPKM") {
        ### mm10: /projects/ps-renlab/y2xie/projects/genome_ref/mm10-2020-A_build/mm10_gene_allowlist.bed
        ### gname needs to be either ensembl or gene_id
        length <- read.table(file = gene_length, header = F)
        length <- length %>% setNames(c("chr", "start", "end", "strand", "ensembl", "gene_id"))
        length$length <- length$end - length$start
        len_mtx <- as.data.frame(length[match(rownames(obj_mtx_collapse), length[,gname]), ])
        readSums <- aggregate(meta$nCount_RNA, list(meta[, group.by]), sum)
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ", length(names(readSums)), "\n")
        cat("check obj_mtx_collapse: ", length(colnames(obj_mtx_collapse)), 
            "\n")
        readSums <- readSums[order(match(names(readSums), colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^9/readSums)/len_mtx$length
    }
    return(obj_collapse_XPM)
}

ArcXPM2 <- function (obj_mtx, meta, method = "CPM", group.by, gname = "gene_id", gene_length = "/projects/ps-renlab/y2xie/projects/genome_ref/mm10-2020-A_build/mm10_gene_allowlist.bed") {
    colnames(obj_mtx) <- meta[colnames(obj_mtx), group.by]
    obj_mtx <- as(obj_mtx, "dgTMatrix")
    uniqnms <- unique(colnames(obj_mtx))
    obj_mtx_collapse <- list()
    for (col_name in uniqnms) {
      cta <- which(colnames(obj_mtx) == col_name)
      obj_mtx_collapse[[col_name]] <- rowSums(obj_mtx[, cta])
    }
    obj_mtx_collapse <- do.call(cbind, obj_mtx_collapse)
    spars <- 1 - sum(obj_mtx_collapse != 0) / prod(dim(obj_mtx_collapse))
    cat(paste0("sparsity: ", spars, "\n"))
    if (method == "CPM") {
        readSums <- aggregate(meta$nCount_histone, list(meta[, group.by]), sum)
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ", length(names(readSums)), "\n")
        cat("check obj_mtx_collapse: ", length(colnames(obj_mtx_collapse)), 
            "\n")
        readSums <- readSums[order(match(names(readSums), colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^6/readSums)
    }
    else if (method == "RPKM") {
        ### mm10: /projects/ps-renlab/y2xie/projects/genome_ref/mm10-2020-A_build/mm10_gene_allowlist.bed
        ### gname needs to be either ensembl or gene_id
        length <- read.table(file = gene_length, header = F)
        length <- length %>% setNames(c("chr", "start", "end", "strand", "ensembl", "gene_id"))
        length$length <- length$end - length$start
        len_mtx <- as.data.frame(length[match(rownames(obj_mtx_collapse), length[,gname]), ])
        readSums <- aggregate(meta$nCount_RNA, list(meta[, group.by]), sum)
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ", length(names(readSums)), "\n")
        cat("check obj_mtx_collapse: ", length(colnames(obj_mtx_collapse)), 
            "\n")
        readSums <- readSums[order(match(names(readSums), colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^9/readSums)/len_mtx$length
    }
    return(obj_collapse_XPM)
}

### calculate big cor: propagate::bigcor ### 
bigcor <- function(x, y, chunks = 2000, method = "pearson"){
    cor <- list()
    nchunks <- ceiling(seq_along(1:ncol(x))/chunks)
    range_chunks <- split(seq_len(ncol(x)), nchunks)
    for(chunk in names(range_chunks)){
        x_qry <- x[,range_chunks[[chunk]]]
        cor[[chunk]] <- cor(x_qry, y, method = method)
    }
    cor <- do.call(rbind, cor)
    return(cor)
}

### calculate cluster overlap score after integration; from Yang Li ### 
cal_ovlpScore <- function(t1, t2){ 
  t1.table <- table(t1)
  t2.table <- table(t2)
  t1.pct <- apply(t1.table, 2, function(x){x/sum(x)})
  t2.pct <- apply(t2.table, 2, function(x){x/sum(x)})
  t1.labels <- colnames(t1.pct)
  t2.labels <- colnames(t2.pct)
  ovlpScore.df <- data.frame(anno1=as.character(), anno2=as.character(), ovlpScore=as.numeric())
  for(t1.label in t1.labels){
    for(t2.label in t2.labels){
      t1.pct.df <- data.frame(t1.pct[,t1.label])
      colnames(t1.pct.df) <- "t1"
      t1.pct.df$ident <- rownames(t1.pct.df)
      t2.pct.df <- data.frame(t2.pct[,t2.label])
      colnames(t2.pct.df) <- "t2"
      t2.pct.df$ident <- rownames(t2.pct.df)
      comp.df <- dplyr::full_join(t1.pct.df, t2.pct.df, by="ident")
      comp.df[is.na(comp.df)] <- 0
      comp.df$ident <- NULL
      comp.df <- t(comp.df)
      ovlpScore <- sum(apply(comp.df, 2, min))
      out <- data.frame(anno1=t1.label, anno2=t2.label, ovlpScore=ovlpScore)
      ovlpScore.df <- rbind(ovlpScore.df, out)
    }
  }
  return(ovlpScore.df)
}


### calculate number of PC used for clustering; this is an implementation of ALLCools significant_pc_test function
significant_pc_test <- function(pc, p_cutoff = 0.1, min_pc = 15, downsample = 50000){
    if(nrow(pc) > downsample){
        pc <- pc[sample(nrow(pc), size = downsample), ]
        print(paste0("Downsample PC matrix to ", downsample, " cells to calculate significant PC components"))
    }
    for (i in 1:(ncol(pc)-1)){
        cur_pc = pc[, i]
        next_pc = pc[, i + 1]
        test_result <- ks.test(cur_pc, next_pc)
        if (test_result$p.value > p_cutoff){
            break
        }
    }
    n_components <- min(i + 1, ncol(pc))
    min_pc <- min(15, ncol(pc))
    if (n_components < min_pc){
        print(paste0('only ', n_components, ' passed the P cutoff. In order to proceed following analysis, will use first ', min_pc, ' PCs'))
        n_components <- min_pc
    } 
    return(n_components)
}

### calculate enrichment of gene with aggregated cell by group matrix
### can be used to calculate QC score with specific set of gene
ssgsea <- function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
    row_names = rownames(X)
    num_genes = nrow(X)
    gene_sets = lapply(gene_sets, function(genes) {which(row_names %in% genes)})

    # Ranks for genes
    R = matrixStats::colRanks(X, preserveShape = T, ties.method = 'average')

    # Calculate enrichment score (es) for each sample (column)
    es = apply(R, 2, function(R_col) {
        gene_ranks = order(R_col, decreasing = TRUE)

        # Calc es for each gene set
        es_sample = sapply(gene_sets, function(gene_set_idx) {
            # pos: match (within the gene set)
            # neg: non-match (outside the gene set)
            indicator_pos = gene_ranks %in% gene_set_idx
            indicator_neg = !indicator_pos

            rank_alpha  = (R_col[gene_ranks] * indicator_pos) ^ alpha

            step_cdf_pos = cumsum(rank_alpha)    / sum(rank_alpha)
            step_cdf_neg = cumsum(indicator_neg) / sum(indicator_neg)

            step_cdf_diff = step_cdf_pos - step_cdf_neg

            # Normalize by gene number
            if (scale) step_cdf_diff = step_cdf_diff / num_genes

            # Use ssGSEA or not
            if (single) {
                sum(step_cdf_diff)
            } else {
                step_cdf_diff[which.max(abs(step_cdf_diff))]
            }
        })
        unlist(es_sample)
    })

    if (length(gene_sets) == 1) es = matrix(es, nrow = 1)

    # Normalize by absolute diff between max and min
    if (norm) es = es / diff(range(es))

    # Prepare output
    rownames(es) = names(gene_sets)
    colnames(es) = colnames(X)
    return(es)
}

### implement of onestep_clust() from hicat
seurat_onestep_clust <- function(mmg, reduction = "pca", var = "none", batch.label = "none", npc = "none", res = "none"){
    cat("Perform log CPM normalization...\n")
    mmg <- NormalizeData(mmg, normalization.method = "LogNormalize", scale.factor = 1e6)
    ### to accomodate seurat5...
    mmg[["RNA"]] <- as(object = mmg[["RNA"]], Class = "Assay")
    norm.dat <- mmg[["RNA"]]@data
    
    cat("Perform variable features selection using Brennecke's method...\n")
    ### pre-filtering on gene used
    filt_genes <- row.names(norm.dat)[which(Matrix::rowSums(norm.dat > 1) >= 10)]
    norm.count <- norm.dat
    norm.count@x = 2^(norm.count@x) - 1
    vg <- find_vg(norm.count[filt_genes,], rescaled = F, return_type = "data", verbose = F)
    sg <- vg[which(vg$loess.padj < 0.5 | vg$dispersion > 3), "gene"]
    cat(paste0("number of variable genes calculated by hicat method: ", length(sg), "\n"))
    
    if (length(sg) > 3000){
        sg <- head(sg[order(vg[sg, "loess.padj"], -vg[sg, "z"])], 3000)
    }

    Misc(object = mmg, slot = "hicat_vg") <- vg
    VariableFeatures(mmg) <- sg

    cat("Perform data scaling and PCA...\n")
    if (var != "none") {
        mmg <- ScaleData(mmg, vars.to.regress = var)
    }else {
        mmg <- ScaleData(mmg)
    }
    mmg <- RunPCA(mmg, features = VariableFeatures(object = mmg), npcs = 50, verbose = F)
    
    ### select top pc
    if (npc == "none") {
        npc <- significant_pc_test(mmg@reductions$pca@cell.embeddings, min_pc = 15)
    }
    cat(paste0("Number of PC selected for downstream analysis: ", npc, "\n")) 
    
    if (batch.label != "none") {
        library(harmony)
        mmg <- RunHarmony(mmg, group.by.vars = batch.label, dims.use = 1:npc)
        reduction <- "harmony"
    }
    
    if(ncol(mmg) >= 500000){
    k <- 50
    }else if(ncol(mmg) >= 200000 & ncol(mmg) < 500000){
        k <- 25
    }else{
        k <- 15
    }
    cat(paste0("Number of neighbors selected for downstream analysis: ", k, "\n"))
    
    cat("Perform UMAP embedding...\n")
    mmg <- RunUMAP(mmg, dims = 1:npc, seed.use = 131, reduction = reduction, 
                   n.neighbors = k, min.dist = 0.1, n.components = 2L, 
                   umap.method = "uwot", uwot.sgd = TRUE, verbose = F)

    if (res == "none") {
        res <- 0.3
    }
    cat(paste0("resolution used for leiden clustering: ", res, "\n")) 

    mmg <- FindNeighbors(object = mmg, k.param = k, reduction = reduction, dims = 1:npc, 
                     nn.method = "annoy", annoy.metric = "cosine", verbose = F)
    
    ### in hicat the resolution for one step clustering is set to be 0.01. 
    ### Disable this setting since it does not give any information.
    ### update 240206: Extremely slow when using leidenalg. Use leidenbase instead for faster clustering.
    suppressWarnings(mmg <- FindClusters(object = mmg, algorithm = 4, method = "igraph", weights = T, resolution = res, verbose = F))

    ### update 240208: do not have an efficient way of dealing with singleton. From now on export RNA_snn and use leidenalg for clustering.

#     cat("Convering SNN graph to igraph...\n")
#     graph <- mmg@graphs$RNA_snn
#     adj.matrix <- Matrix::sparseMatrix(i = graph@i+1, p = graph@p, x = graph@x, dims = graph@Dim, dimnames = graph@Dimnames)
#     sgraph <- igraph::graph_from_adjacency_matrix(adjmatrix = adj.matrix, weighted = TRUE, diag = FALSE)
#     sgraph <- igraph::simplify(sgraph, remove.multiple = T)
    
#     cat("Performing leiden clustering with leidenbase...\n")
#     lres <- leidenbase::leiden_find_partition(sgraph, resolution_parameter = res, num_iter = 10, partition_type = "RBConfigurationVertexPartition", seed = 921)
#     cat(paste0("Modulaity: ", lres$modularity, "\n"))
#     mmg$leiden <- lres$membership
    return(mmg)
}

### Preprocess with onestep_clust()
ProcessRNA_hicat <- function(rna_path, valid){
    ### rna_path: path to 10X generated RNA raw counts
    ### valid: valid cell table, contains at least column recording RNA barcode information, with name "rna"
    mmg <- Read10X(rna_path)
    if(!("rna" %in% colnames(valid))){
        stop("rna barcode information is not found in valid barcodes file. Make sure column with name 'rna' exists.")
    }
    if ("atac_bc" %in% colnames(valid)){
        cat("The provided barcode file seems to be a multiplexed reaction. Determine multiplet based on the observed collision rate...\n")
        nondup <- valid$rna[!duplicated(valid$rna) & !duplicated(valid$rna, fromLast = TRUE)] ### RNA barcodes appear once
        dup <- valid$rna[duplicated(valid$rna) | duplicated(valid$rna, fromLast = TRUE)] %>% unique ### RNA barcodes appear multiple times
        cat("Number of RNA barcodes matched with multiple ATAC barcodes: ", length(unique(dup)), "\n")
        observed_collide <- 100 * (length(unique(valid$rna)) - length(nondup))/length(unique(valid$rna))
        cat(paste0("observed collision rate:", round(observed_collide, 2), "%\n"))
        # estimated_collide <- 100 * (length(unique(valid$rna)) * 0.08)/10000
        # collide <- estimated_collide - observed_collide
        preindex <- substr(valid$atac_bc, 1, (nchar(valid$atac_bc[1]) - 16))
        collide_ratio <- sum(table(preindex)^2)/(sum(table(preindex))^2 - sum(table(preindex)^2))
        collide <- observed_collide * collide_ratio
        if(is.na(collide)){
            collide <- 0
        }
        part <- max(0.2, collide/100) ### update 240709: change to default
        cat(paste0("Estimated doublets rate after removing explict multiplet: ", round(collide, 2), "%\n"))
    }else if (length(unique(valid$rna)) != nrow(valid)){
        stop("Duplicated RNA barcodes detected, but the valid barcode file is not for multiplexed sample.")
    }else{
        nondup <- unique(valid$rna)
        collide <- 100 * (length(nondup) * 0.08)/10000
        part <- max(0.20, collide/100)
    }
    mmg <- mmg[, nondup]
    mmg <- CreateSeuratObject(mmg)
    
    cat(paste0("cells number after paired segment: ", length(nondup), "\n"))
    cat(paste0("start preprocessing...\n"))
    mmg$percent.mt <- PercentageFeatureSet(mmg, pattern = "^MT-")
    mmg$percent.ribo <- PercentageFeatureSet(mmg, pattern = "(^RPL|^RPS|^MRP)")
    gene_lcf <- median(log10(mmg$nFeature_RNA)) - 3*mad(log10(mmg$nFeature_RNA))
    gene_hcf <- median(log10(mmg$nFeature_RNA)) + 3*mad(log10(mmg$nFeature_RNA))
    count_hcf <- median(log10(mmg$nCount_RNA)) + 3*mad(log10(mmg$nCount_RNA))
    
    cat(paste0(#"gene lowest cutoff: ", max(300, round(10^gene_lcf)), 
        "\ngene highest cutoff: ", round(10^gene_hcf), 
        "\nreads highest cutoff: ", min(50000, round(10^count_hcf))), "\n")

    ### update 240709: to consider very low count population observed in MERFISH
    mmg <- subset(mmg, subset = nCount_RNA < min(50000, 10^count_hcf) & nFeature_RNA < 10^gene_hcf & percent.mt < 10) # & nFeature_RNA > max(100, 10^gene_lcf)
    cat(paste0("cells number after quality segment: ", ncol(mmg), "\n"))
    cat(paste0("start processing RNA data...\n"))
    mmg <- seurat_onestep_clust(mmg)
    
    cat(paste0("start identifying doublets...\n"))
    cat(paste0("estimated doublets rate: ", collide, "%\n"))

    ### update 240709: perform 50 iteration to identify constant doublets
    library(foreach)
    library(doParallel)
    cl <- makeCluster(12)
    registerDoParallel(cl)
    dscores <- foreach(i = 1:50, .combine = 'c') %dopar% {
        library(Seurat)
        source("/projects/ps-renlab2/y2xie/scripts/hicat_doubletfinder_mo.R")
        source("/projects/ps-renlab2/y2xie/scripts/hicat_vg.R")
        source("/projects/ps-renlab2/y2xie/scripts/DPT_help.R")
        score <- doubletfinder(mmg, VariableFeatures(mmg), proportion.artificial = part, plot = F)
        return(score)
    }
    stopCluster(cl)
    
    even_indices <- seq(from = 1, to = length(dscores), by = 2)
    dconstant <- Reduce(intersect, lapply(dscores[even_indices], function(x) names(x)[which(x > 0.2)]))
    avscore <- do.call(cbind, lapply(dscores[even_indices], function(x) as.data.frame(x))) %>%
    setNames(c(1:10)) %>% 
    rowMeans() %>% as.data.frame %>% 
    setNames("avg_dscore") %>% tibble::rownames_to_column("bc") %>%
    mutate(group = ifelse(bc %in% dconstant, "Constant", "Varied")) %>%
    mutate(group = factor(group, levels = c("Varied", "Constant")))

    # doublet <- doubletfinder(mmg, VariableFeatures(mmg), proportion.artificial = part, plot = T)
    # cat(paste0("calculated doublet rate: ", 100*(table(doublet[[1]] > 0.25)/length(doublet[[1]]))[2], "%\n"))

    mmg$doublet_score <- avscore$avg_dscore #doublet[[1]]
    mmg$doublet_ident <- ifelse(colnames(mmg) %in% dconstant, "Doublet", "Singlet") #ifelse(doublet[[1]] > 0.25, "Doublet", "Singlet")
    return(mmg)
}

### Find DE genes for clusters, adopted from Zemke 2023 Nature
### need to be tested
run.edgeR <- function(counts, sample, group, prefix, pairwise = FALSE){
    library("edgeR")

    ### counts: a gene by group matrix. Can obtained through AggregateExpression() 
    ### sample: a dataframe recording sample and group information. 
    ### group: name of group information in sample (group %in% colnames(sample))
    ### prefix: saving file as...
    cat("now creating DGE class for edgeR...\n")
    group.list <- factor(sample[, group])
    dge = DGEList(counts = counts, group = group.list, samples = sample)
    keep = filterByExpr(dge)
    dge = dge[keep, ,keep.lib.sizes=FALSE]
    dge = calcNormFactors(dge, 'TMM')
    design = model.matrix(~0 + group, data = dge$samples)
    colnames(design) <- make.names(gsub(group, "", colnames(design)))
    dge = estimateDisp(dge, trend.method = 'locfit', design = design)
    fit = glmFit(dge, design = design)
    
    cat("now performing glmLRT for individual clusters...\n")
    coef <- unique(group.list)
    for (f in coef){
        ### Perform an ANOVA like test for all groups
        cat(paste0("cluster: ", f, "\n"))
        one_v_rest = c(rep(-1/(length(coef)-1), times = length(coef)))
        one_v_rest[which(coef == f)] = 1  
        vrest = glmLRT(fit, contrast = one_v_rest)
        write.table(topTags(vrest, n = "Inf")$table, file = paste0(prefix, "_edgeR_", f, "_vs_all_groups.xls"),
                    sep = "\t", quote = F)
        cat(paste0("cluster: ", f, " finish.\n"))
        if(pairwise){
            cat(paste0("Now calculate pairwise test for cluster: ", f, "\n"))
            for (s in coef){
            ### Perform pairwise test for each group
                if (s == f){
                        next
                }
                one_v_one = c(rep(0, times = length(coef)))
                one_v_one[which(coef == f)] = 1  
                one_v_one[which(coef == s)] = -1 
                vtest = glmLRT(fit, contrast = one_v_one)
                write.table(topTags(vtest, n = "Inf")$table, file = paste0(prefix, "_edgeR_", f, "_vs_", s, ".xls"),
                            sep = "\t", quote = F)
            }
        }
    }
}

### Calculate tissue / cell type specificity of marker genes
### adopted from hicat (https://rdrr.io/github/AllenInstitute/scrattch.hicat/src/R/util.R#sym-calc_tau)
calc_tau <- function(m) {  
  row_maxes <- apply(m, 1, max)
  
  m <- m / row_maxes
  tau <- Matrix::rowSums(1 - m) / (ncol(m) - 1)
  tau[is.na(tau)] <- 0
  return(tau)
}
