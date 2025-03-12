library(Matrix)
suppressPackageStartupMessages(library(dplyr))
library(ggplot2)
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggbreak))
library(reshape2)
# source("/projects/ps-renlab/y2xie/scripts/basics.R")

get.elbow.points.indices <- function(x, y, threshold = 0.0005) {
    d1 <- diff(y) / diff(x) # first derivative
    d2 <- diff(d1) / diff(x[-1]) # second derivative
    indices <- which(abs(d2) > threshold)  
    return(indices)
}
  
get.cells.10X <- function(y, n = 4000) { # n: expected number of cells
    d1 <- quantile(y[1:n], probs = 0.99)
    threshold <- 0.1*d1
    indices <- which.min(y > threshold)
    return(indices)
}

PHCspeciesPlot <- function(obj, valid_cells = "none", prefix){
    ### provided reads mapped to two species (mm10 & hg38) for selecting valid cells ###
    ### /home/y2xie/scratch/33.Paired_HiC_NovaSeq_221206/scripts/phc.count_human_and_mouse_reads.py ###
    if(valid_cells[1] == "none"){
        valid_cells <- rownames(obj)
    }
    obj <- obj %>% mutate(mm10_ratio = mm10 / (mm10 + GRCh38)) %>% mutate(GRCh38_ratio = GRCh38 / (mm10 + GRCh38))
    pobj <- obj[valid_cells,]
    ss_cells <- rownames(pobj[pobj$mm10_ratio > 0.8 | pobj$GRCh38_ratio > 0.8,])
    label <- data.frame(anno = paste0("PF cells: ", length(valid_cells), "\n", 
                                      "Cells assigned to one species: ", length(ss_cells), "\n",
                                      "% Cells assigned to one species: ", 100*round(length(ss_cells)/length(valid_cells), 4)))
    fig <- pobj %>%
    ggplot(aes(x = mm10, y = GRCh38)) + 
    geom_point(size = 0.5, color = "grey") + 
    geom_point(data = pobj[pobj$mm10_ratio > 0.8, ], aes(x = mm10, y = GRCh38), size = 0.5, color = "#477A96") + 
    geom_point(data = pobj[pobj$GRCh38_ratio > 0.8, ], aes(x = mm10, y = GRCh38), size = 0.5, color = "#AE4237") +
    theme_classic() + 
    ylab("Reads mapped to GRCh38") + 
    xlab("Reads mapped to mm10") + 
    xlim(0, min(max(pobj$mm10), max(pobj$GRCh38))) + 
    ylim(0, min(max(pobj$mm10), max(pobj$GRCh38))) +
    geom_text(data = label, size = 3, aes(x = min(max(pobj$mm10), max(pobj$GRCh38)), y = min(max(pobj$mm10), max(pobj$GRCh38)), 
                  hjust = 1, vjust = 1, label = anno))
    # ggsave(fig, filename = paste0(prefix, "_species.banyard.pdf"), 
    #        height = 5, width = 5, dpi = 300)
    ggsave(fig, filename = paste0(prefix, "_species.banyard.png"), 
           height = 5, width = 5, dpi = 300)
    write.table(as.data.frame(ss_cells), file = paste0(prefix, "_species_clean.cells.txt"),
                row.names = F, col.names = F, sep = "\t", quote = F)
    return(fig)
}

PHCrankPair <- function(obj, count_cutoff = "none", valid_cells = "none", prefix, method = "elbow"){
    ### provided pair count results for selecting valid cells ###
    ### /projects/ps-renlab/y2xie/scripts/Paired-HiC/phc.count_pairs_sc.py
    ### PairCount.stat.csv: c("total", "mapped", "unmapped", "duplicate", "cis", "cis_1kb-", "cis_1kb+", "cis_10kb+", "trans")
    obj_rank <- obj[, "total", drop = F] %>% as.data.frame() %>% setNames("count") %>% 
    arrange(desc(count)) %>% 
    tibble::rownames_to_column("bc") %>% 
    tibble::rowid_to_column("rank")
    if(count_cutoff != "none" & valid_cells != "none"){
        stop("only specify valid cells group or cutoff, but not both.")
    }
    if(count_cutoff == "none" & valid_cells == "none"){ ### automatically find cutoff!
        if(method == "elbow"){
            obj_rank_subset <- obj_rank[50:50000,]
            elbow <- get.elbow.points.indices(obj_rank_subset$rank, obj_rank_subset$count, threshold = 0.0005)
            cutoff <- min(round(median(elbow)), round(mean(elbow)))
            count_cutoff <- obj_rank_subset[cutoff, "count"]
        }else if(method == "10X"){
                obj_rank_subset <- obj_rank[50:50000,]
                cutoff <- get.cells.10X(obj_rank_subset$count, n = 7000)
                count_cutoff <- obj_rank_subset[cutoff, "count"]
            }       
    }
    valid <- obj_rank[obj_rank$count > count_cutoff, "bc"]
    if(valid_cells != "none"){ ### disable finding cutoff
        valid <- valid_cells
    }
    
    label <- data.frame(anno = paste0("Cutoff:", count_cutoff, "\n", 
                                      "PF cells: ", nrow(obj_rank[obj_rank$bc %in% valid,]), "\n", 
                                      "Median pairs: ", as.integer(median(obj_rank[obj_rank$bc %in% valid, "count"]))))
    fig <- obj_rank %>% ggplot(aes(x = rank, y = count)) +
    geom_point(size = 0.5, color = "grey") + 
    geom_point(data = obj_rank[obj_rank$bc %in% valid,], aes(x = rank, y = count), 
               color = colfunc2(3)[[2]], size = 0.5) + 
    theme_bw() + 
    geom_text(data = label, size = 3, aes(x = min(obj_rank$rank), y = min(obj_rank$count), 
                  hjust = 0, vjust = -0.2, label = anno)) + 
    scale_x_log10() +
    scale_y_log10()
    # ggsave(fig, filename = paste0(prefix, "_cells.select.pdf"), height = 4, width = 4, dpi = 300)
    ggsave(fig, filename = paste0(prefix, "_cells.select.png"), height = 4, width = 4, dpi = 300)
    return(valid)
}

PHCplotPair <- function(obj, count_min = 0, count_max = 5000000, prefix){
    ### provided pair count results for plot ###
    ### PairCount.stat.csv: c("total", "mapped", "unmapped", "duplicate", "cis", "cis_1kb-", "cis_1kb+", "cis_10kb+", "trans")
    t1 <- obj %>% filter(total > count_min & total < count_max) %>% 
    mutate("cis1k-%" = 100* `cis_1kb-` / total, 
           "trans%" = 100* trans / total, 
           "cis1k+%" = 100* `cis_1kb+` / total,
           "cis10k+%" = 100* `cis_10kb+` / total)
    
    p1 <- t1 %>% dplyr::select(c("cis1k-%", "trans%", "cis1k+%", "cis10k+%")) %>% 
    as.matrix() %>% 
    reshape2::melt() %>%
    ggplot(aes(x=Var2, y=value)) + 
    geom_violin(aes(fill = Var2, color = Var2), alpha = 0.5) + 
    geom_boxplot(width=0.1, color="black", alpha=0.8, coef = 2, outlier.shape = NA) +
    scale_color_manual(values = colfunc(6)) + 
    scale_fill_manual(values = colfunc(6)) + 
    xlab("")
    
    p2 <- t1 %>% dplyr::select(c("cis_1kb-", "trans", "cis_1kb+", "cis_10kb+")) %>% 
    as.matrix() %>% 
    reshape2::melt() %>%
    ggplot(aes(x=Var2, y=log10(value))) + 
    geom_violin(aes(fill = Var2, color = Var2), alpha = 0.5) + 
    geom_boxplot(width=0.1, color="black", alpha=0.8, coef = 2, outlier.shape = NA) +
    scale_color_manual(values = colfunc(6)) + 
    scale_fill_manual(values = colfunc(6)) + 
    xlab("")
    
    fig <- p1 + p2 + plot_layout(ncol = 2)
    ggsave(fig, filename = paste0(prefix, "_pairs.summary.pdf"), height = 6, width = 12, dpi = 300)
    # ggsave(fig, filename = paste0(prefix, "_pairs.summary.png"), height = 6, width = 12, dpi = 300)
    return(fig)
}

### seruat-like findmarker for single-cell compartment score
diffComp <- function(data, meta, name_col = "barcode", idents_col = "cluster", 
                     ident.1 = NULL, ident.2 = NULL, p_adj = 0.05){
    ### two columns meta: cells, cluster
    new_name <- meta[match(colnames(data), meta[, name_col]), idents_col]
    if(!is.null(new_name)){
        colnames(data) <- new_name
        if(!is.null(ident.2)){
            data <- data[, colnames(data) %in% c(ident.1, ident.2)]
        }
    }else{
        print("name column are different from data colnames.")
        break
    }
    # rtk <- (apply(data, MARGIN = 1, function(x) max(x) - min(x)) > 0.5) ### filter variable bin!
    rtk <- apply(data, MARGIN = 1, function(x) (sd(x) / mean(x)) * 100)
    cutoff <- median(log10(rtk), na.rm = TRUE)
    data <- data[!is.na(rtk) & (rtk > 10^cutoff), ]
    if(is.null(ident.2)){
        diffcomp <- list()
        for (f in unique(colnames(data))){ ### this requires colnames to be corced
            diffb <- list()
            print(paste0("calculated differential compartment for cluster ", f))
            qry <- data[, colnames(data) %in% f, drop = F] 
            ref <- data[, colnames(data) != f, drop = F]
            for (r in 1:nrow(data)){
                tqry <- unlist(qry[r,])
                tref <- unlist(ref[r,])
                trank <- wilcox.test(tqry, tref)
                tdiff_median <- median(tqry) - median(tref)
                tdiff_mean <- mean(tqry) - mean(tref)
                diffb[[r]] <- data.frame(bin = rownames(data)[r],
                                         p_value = trank$p.value,
                                         # p_value_adj = nrow(data)*(trank$p.value),
                                         median_diff = tdiff_median,
                                         mean_diff = tdiff_mean,
                                         cluster = f)
            }
            diffb <- do.call(rbind, diffb) 
            diffb$p_value_adj <- p.adjust(diffb$p_value, method = "BH")
            diffb <- diffb %>% 
            dplyr::filter(p_value < 0.05)
            diffcomp[[f]] <- diffb
        }
        diffcomp <- do.call(rbind, diffcomp)
    }
    else{
        print(paste0("calculated differential compartment for cluster ", ident.1, " vs cluster ", ident.2))
        diffb <- list()
        qry <- data[, colnames(data) %in% c(ident.1), drop = F] 
        ref <- data[, colnames(data) %in% c(ident.2), drop = F]
        for (r in 1:nrow(data)){
            tqry <- unlist(qry[r,])
            tref <- unlist(ref[r,])
            trank <- wilcox.test(tqry, tref)
            tdiff_median <- median(tqry) - median(tref)
            tdiff_mean <- mean(tqry) - mean(tref)
            diffb[[r]] <- data.frame(bin = rownames(data)[r],
                                     p_value = trank$p.value,
                                     # p_value_adj = nrow(data)*(trank$p.value),
                                     median_diff = tdiff_median,
                                     mean_diff = tdiff_mean)
        }
        diffb <- do.call(rbind, diffb) 
        diffb$p_value_adj <- p.adjust(diffb$p_value, method = "BH")
        diffcomp <- diffb # %>% dplyr::filter(p_value < 0.05)
    }
    return(diffcomp)
}


### find significant trans interaction give the query cooler dump matrix
### 'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2'
inter_p <- function(tmtx, i){
    M <- nrow(tmtx) ### 1/M is the probs of trans interaction assuming visibility is the same
    n = sum(tmtx$count) ### n is the total number of observed reads between two chrom
    k = as.integer(tmtx[i, "count"]) ### k is the observed number of contacts between two bin
    p_val <- dbinom(x = k:n, size = n, prob = 1/M) %>% sum
    return(p_val)
}

inter_p_cnv <- function(tmtx, i){
    ### perform cnv normalization with cnv table
    M <- sum(tmtx$residuals) ### 1/M is the probs of trans interaction assuming visibility depends on copy number (log2 by default)
    n = sum(tmtx$count) ### n is the total number of observed reads between two chrom
    k = as.integer(tmtx[i, "count"]) ### k is the observed number of contacts between two bin
    probs = tmtx$residuals[i]/M
    p_val <- dbinom(x = k:n, size = n, prob = probs) %>% sum
    return(p_val)
}

FindTrans <- function(tmtx, p_cutoff = 0.01, chrom_norm = TRUE, 
                          chromsize = "/projects/ps-renlab2/y2xie/projects/genome_ref/hg38.main.chrom.sizes") {
    ### input: cooler dump contacts info, with chrom1-start1 being the query bin. 
    ### will calculate and test for all chromosomes
    chromsize = read.table(chromsize)
    qbin <- paste0(unique(tmtx$chrom1), ":", unique(tmtx$start1))
    print(paste0("find significant trans interaction for bin ", qbin))
    print(paste0("p value cutoff after correction: ", p_cutoff))
    tmtx <- tmtx %>% mutate(bin = paste0(chrom2, ":", start2, "-", end2))

    tmtx$p_value <- NA
    tmtx$p_adj <- NA
    tmtx$sig_trans <- "False"

    print("start calculation...")
    for (f in unique(tmtx$chrom2)){
        ttmtx <- tmtx[tmtx$chrom2 == f, ]
        # ccnv <- cnv[cnv$chrom == f, ]
        p_list <- list()
        for (i in 1:nrow(ttmtx)){
            p_list[[i]] <- inter_p(ttmtx, i)
        }
        tmtx[tmtx$chrom2 == f, ]$p_value <- unlist(p_list)
        if(chrom_norm == TRUE){
            len1 <- chromsize[chromsize$V1 == unique(tmtx$chrom1), "V2"] %>% as.numeric
            len2 <- chromsize[chromsize$V1 == f, "V2"] %>% as.numeric
            max1 <- chromsize[chromsize$V1 == "chr1", "V2"] %>% as.numeric
            max2 <- chromsize[chromsize$V1 == "chr2", "V2"] %>% as.numeric
            ### normalize for chrom size
            tmtx[tmtx$chrom2 == f, ]$p_adj <- p.adjust(unlist(p_list), method = "BH")*len1*len2/(max1*max2)         
        }else{
            tmtx[tmtx$chrom2 == f, ]$p_adj <- p.adjust(unlist(p_list), method = "BH")
        }        
        if(length(which(tmtx[tmtx$chrom2 == f, ]$p_adj < p_cutoff)) > 0){
            tmtx[tmtx$chrom2 == f & tmtx$p_adj < p_cutoff, ]$sig_trans <- "True"
        }
        print(paste0("chromosome: ", f, ", significant trans interactions with ", qbin, ": ", nrow(tmtx[tmtx$chrom2 == f & tmtx$p_adj < p_cutoff, ])))
    }
    print("finish calculation...")
    print(paste0("Total significant trans interaction with ", qbin, " :", nrow(tmtx[tmtx$sig_trans == "True",])))
    return(tmtx)
}

FindTransNorm <- function(tmtx, p_cutoff = 0.01, chrom_norm = TRUE, 
                          chromsize = "/projects/ps-renlab2/y2xie/projects/genome_ref/hg38.main.chrom.sizes") {
    chromsize = read.table(chromsize)
    ### input: cooler dump contacts info, with chrom1-start1 being the query bin. 
    ### will calculate and test for all chromosomes
    qbin <- paste0(unique(tmtx$chrom1), ":", unique(tmtx$start1))
    print(paste0("find significant trans interaction for bin ", qbin, " with cnv normalization"))
    print(paste0("p value cutoff after correction: ", p_cutoff))
    if (!("residuals" %in% colnames(tmtx))) {
        stop(paste("Warning: the residuals column is not found in the input."))
    }

    tmtx <- tmtx %>% mutate(bin = paste0(chrom2, ":", start2, "-", end2))

    tmtx$p_value <- NA
    tmtx$p_adj <- NA
    tmtx$sig_trans <- "False"

    print("start calculation...")
    for (f in unique(tmtx$chrom2)){
        ttmtx <- tmtx[tmtx$chrom2 == f, ]
        # ccnv <- cnv[cnv$chrom == f, ]
        p_list <- list()
        for (i in 1:nrow(ttmtx)){
            p_list[[i]] <- inter_p_cnv(ttmtx, i)
        }
        tmtx[tmtx$chrom2 == f, ]$p_value <- unlist(p_list)
        if(chrom_norm == TRUE){
            len1 <- chromsize[chromsize$V1 == unique(tmtx$chrom1), "V2"] %>% as.numeric
            len2 <- chromsize[chromsize$V1 == f, "V2"] %>% as.numeric
            max1 <- chromsize[chromsize$V1 == "chr1", "V2"] %>% as.numeric
            max2 <- chromsize[chromsize$V1 == "chr2", "V2"] %>% as.numeric
            ### normalize for chrom size
            tmtx[tmtx$chrom2 == f, ]$p_adj <- p.adjust(unlist(p_list), method = "BH")*len1*len2/(max1*max2)         
        }else{
            tmtx[tmtx$chrom2 == f, ]$p_adj <- p.adjust(unlist(p_list), method = "BH")
        }        
        if(length(which(tmtx[tmtx$chrom2 == f, ]$p_adj < p_cutoff)) > 0){
            tmtx[tmtx$chrom2 == f & tmtx$p_adj < p_cutoff, ]$sig_trans <- "True"
        }
        print(paste0("chromosome: ", f, ", significant trans interactions with ", qbin, ": ", nrow(tmtx[tmtx$chrom2 == f & tmtx$p_adj < p_cutoff, ])))
    }
    print("finish calculation...")
    print(paste0("Total significant trans interaction with ", qbin, " :", nrow(tmtx[tmtx$sig_trans == "True",])))
    return(tmtx)
}

ssgsea = function(X, gene_sets, alpha = 0.25, scale = T, norm = F, single = T) {
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
