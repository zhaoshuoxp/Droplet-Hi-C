
## set colors
library(dichromat)
# library(wesanderson)
library(scales)

#   
kanagawa <- c("#c1a06e", "#d9d1ba", "#82b0b2", "#346b84", "#1a3657")
tori <- c("#477a96", "#c3533d", "#85212b", "#bcc9d1", "#6a9d84")
koshu <- c("#224D4D", "#A8B38D", "#FFD37C", "#EFD6A7", "#609092")
zero2 <- c("#B14846", "#98565E", "#76545D", "#FDC9C9", "#EDECEC")
ichigo <- c("#253245", "#37517C", "#7FA1BD", "#95C2DF", "#D7D4E4")
ukiyo <- c("#6f4066", "#985e8b", "#fefbe4", "#1A3657", "#346B84", "#BCC9D1", "#8c9c9e", "#224D4D", "#515547", "#6A9D84", "#A8B38D", "#FFD37C", "#e7ba69", "#c56742", "#cbb698", "#af7d7e", "#B14846", "#701b28")

colfunc <- colorRampPalette(c(ukiyo))
colfunc2 <- colorRampPalette(c(tori))
rred <- colorRampPalette(c(zero2))
bblue <- colorRampPalette(c(ichigo))

### formatting theme
### https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=14, base_family="Helvetica") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.margin = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)

}

### Paired-Tag DNA matrix processing
RunDNAMask <- function(x, maskMT=TRUE, maskBL=TRUE, maskSex=FALSE){
  if(maskMT){
    cat("there are",nrow(x),"bins before mt bins masking\n", sep = " ")
    mt_bins <- grep(rownames(x), pattern = "^chrM"); x <- x[-mt_bins,]
    cat("there are",nrow(x),"bins after mt bins masking\n", sep = " ")
  } else{cat("no mt masking is done.\n")}
  if(maskBL){
    cat("read in encode blacklist bed file...\n")
    black_list = read.table("/projects/ps-renlab/y2xie/projects/genome_ref/mm10.blacklist.bed")
    black_list.gr = GenomicRanges::GRanges(black_list[,1], IRanges::IRanges(black_list[,2], black_list[,3]))
    idy = S4Vectors::queryHits(GenomicRanges::findOverlaps(GenomicRanges::GRanges(rownames(x)), black_list.gr))
    cat("there are",nrow(x),"bins before blacklist masking\n", sep = " ")
    if(length(idy) > 0){x <- x[-idy,]}
    cat("there are",nrow(x),"bins after blacklist masking\n", sep = " ")
  } else{cat("no encode blacklist masking is done.\n")}
  if(maskSex){
    cat("there are",nrow(x),"bins before sex chr bins masking\n", sep = " ")
    sex_bins <- grep(rownames(x), pattern = paste("^chrX","^chrY",sep="|")); x <- x[-sex_bins,]
    cat("there are",nrow(x),"bins after sex chr bins masking\n", sep = " ")
  } else{cat("no sex chromosome masking is done.\n")}
  return(x)
}

### simple Seurat
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

### Paired-Tag manual filtering 
read<-function(dna, rna, genome){
  dna.m<-readMM(paste(path,dna,"_",genome,"_sorted_rmdup_mtx2/matrix.mtx", sep=""))
  dna.bc<-read.csv(paste(path,dna,"_",genome,"_sorted_rmdup_mtx2/barcodes.tsv", sep=""), sep="\t", head=F)
  sdna<-colSums(dna.m)
  names(sdna)<-dna.bc[,1]
  rna.m<-readMM(paste(path,rna,"_",genome,"_sorted_rmdup_mtx2/matrix.mtx", sep=""))
  rna.bc<-read.csv(paste(path,rna,"_",genome,"_sorted_rmdup_mtx2/barcodes.tsv", sep=""), sep="\t", head=F)
  srna<-colSums(rna.m)
  names(srna)<-rna.bc[,1]
  m<-merge(sdna,srna,by=0)
  
  gc()
  return(m)
}

###########################################################
#### barcodes selection mathods
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

get.dist.intersect <- function(x){ # find distribution intersection of cells and empty barcodes
    if(min(x) < 1){ # no way
      cutoff = 1
    }else{cutoff = min(x)}
    mixx <- normalmixEM(x[x > cutoff])
    d1 <- mixx$lambda[1] * dnorm(x, mean = mixx$mu[1], sd = mixx$sigma[1])
    d2 <- mixx$lambda[2] * dnorm(x, mean = mixx$mu[2], sd = mixx$sigma[2])
    range <- which(x > min(mixx$mu) & x < max(mixx$mu))
    indices <- which.min(abs(d1[range] - d2[range]))
    return(indices)
}
  ###########################################################
  
plot_mtx<-function(m,dna.cutoff=100,rna.cutoff=100,dna.u=dna.cutoff*20,rna.u=rna.cutoff*20){
  plot(m[,2:3], pch=19, cex=0.25, col="grey", log="xy", xlab="# of DNA reads", ylab="# of RNA reads", main=paste(dna,rna,sep="_"))
  m.f<-m[m[,2]>dna.cutoff & m[,3]>rna.cutoff & m[,2]<dna.u & m[,3]<rna.u,]
  ncell<-dim(m.f)[1]
  points(m.f[,2:3], pch=19, cex=0.25, col="firebrick")
  med.dna<-median(m.f[,2])
  med.rna<-median(m.f[,3])
  lines(c(dna.cutoff,dna.cutoff),c(1,1000000), lwd=2, lty=2, col="grey")
  lines(c(1,1000000),c(rna.cutoff,rna.cutoff), lwd=2, lty=2, col="grey")
  lines(c(dna.u,dna.u),c(1,1000000), lwd=2, lty=2, col="grey")
  lines(c(1,1000000),c(rna.u,rna.u), lwd=2, lty=2, col="grey")
  legend("bottomright", legend=c(paste("Median DNA reads: ", med.dna, sep=""), paste("Median RNA reads: ", med.rna, sep=""), paste("# of cells: ", ncell, sep="")),bty="n")
  return(m.f)
}

write_PF<-function(n, prefix="NA"){
  id<-n[,1]
  cid<-paste(prefix, id, sep=":")
  out<-cbind(id, cid, n[,2:3])
  colnames(out)<-c("Cell_ID_raw", "Cell_ID_cov", "DNA_reads", "RNA_reads")
  write.table(out, sep="\t", col.names=T, row.names = F, quote=F, file=paste(R_path, dna,"_",rna,"_PF_cells.xls", sep=""))
}

get.elbow.points.indices <- function(x, y, threshold) { # x is rank / indices, y is reads count 
    d1 <- diff(y) / diff(x) # first derivative
    d2 <- diff(d1) / diff(x[-1]) # second derivative
    indices <- which(abs(d2) > threshold)  
    return(indices)
  }

### to collapase sparse matrix based on group labels
library(Matrix)
OP2 <- function(x) {
    nms <- colnames(x)
    uniquenms <- unique(nms)
    # build the sparseMatrix again: x's with same index values are automatically
    # added together, keeping in mind that indexes stored from 0 but built from 1
    sparseMatrix(i = x@i + 1, 
                 j = match(nms, uniquenms)[x@j + 1],
                 x = x@x,
                 dimnames = list(rownames(x), uniquenms),
                 repr = "C")
}

RunXPM <- function(obj, meta, method = "CPM", group.by = "Stage", group.by2 = NULL) 
{
    cat("Check whether nCount_histone and nCount_RNA exists in metadata before running! \n")
    obj_mtx <- GetAssayData(obj, slot = "counts")
    mark <- meta[colnames(obj_mtx)[1],]$Target
    if (is.null(group.by2)){
       colnames(obj_mtx) <- obj@meta.data[, group.by] 
    } else {
       colnames(obj_mtx) <- paste0(obj@meta.data[,group.by], "_", obj@meta.data[,group.by2]) 
    }
    if (method == "CPM") {meta <- meta[meta$Target == mark, , drop=F]}
    obj_mtx <- as(obj_mtx, "dgTMatrix")
    obj_mtx_collapse <- OP2(obj_mtx)
    spars <- length(obj_mtx_collapse@x)/obj_mtx_collapse@Dim[1]/obj_mtx_collapse@Dim[2]
    cat(paste0("sparsity: ", spars, "\n"))
    if (spars > 0.2) {
        obj_mtx_collapse <- as(obj_mtx_collapse, "matrix")
        cat("coarse dgTMatrix into Matrix.\n")
    }

    if (method == "CPM") { # oh please I need total reads
        if (is.null(group.by2)){
            readSums <- aggregate(meta$nCount_histone, list(meta[, group.by]), sum)
        } else{
            readSums1 <- aggregate(meta$nCount_histone, list(meta[, group.by], meta[, group.by2]), sum)
            readSums <- data.frame(tmp = paste0(readSums1$Group.1, "_", readSums1$Group.2),
                      sums = readSums1$x) 
        }
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ",length(names(readSums)),"\n")
        cat("check obj_mtx_collapse: ",length(colnames(obj_mtx_collapse)),"\n")
        readSums <- readSums[order(match(names(readSums),colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^6/readSums)
        
    }
    else if (method == "RPKM") {
        length <- read.table(file = "/projects/ps-renlab/y2xie/projects/genome_ref/Paired-Tag/mm10/mm10.PairedTag.txt", 
            header = F)
        len_mtx <- as.data.frame(length[match(rownames(obj_mtx_collapse), 
            length$V5), ])
        len_mtx$length <- len_mtx$V3 - len_mtx$V2
        if (is.null(group.by2)){
            readSums <- aggregate(meta$nCount_RNA, list(meta[, group.by]), sum)
        } else{
            readSums1 <- aggregate(meta$nCount_RNA, list(meta[, group.by], meta[, group.by2]), sum)
            readSums <- data.frame(tmp = paste0(readSums1$Group.1, "_", readSums1$Group.2),
                      sums = readSums1$x) 
        }
        colnames(readSums) <- c("tmp", "sums")
        readSums <- setNames(as.numeric(readSums$sums), readSums$tmp)
        cat("check readSums: ",length(names(readSums)),"\n")
        cat("check obj_mtx_collapse: ",length(colnames(obj_mtx_collapse)),"\n")        
        readSums <- readSums[order(match(names(readSums),colnames(obj_mtx_collapse)))]
        obj_collapse_XPM <- t(t(obj_mtx_collapse) * 10^9/readSums)/len_mtx$length
    }
    return(obj_collapse_XPM)
}

### Taken from SnapATAC, clustering using jacaard index
# calculate Jaccard index and normalization (for coverage bias), borrowed from SnapATAC. Also embedded in Paired_Map
# calculate Jaccard index without downsampling
calJaccard<-function(X_i, X_j){
  ### input m = features
  ### input n = cells
  ### transpose here
  X_i<-t(X_i)
  X_j<-t(X_j)
  A = Matrix::tcrossprod(X_i, X_j)
  bi = Matrix::rowSums(X_i)
  bj = Matrix::rowSums(X_j)
  jmat = as.matrix(A / (replicate(ncol(A), bi) + t(replicate(nrow(A), bj)) - A))
  ### output m = features
  ### output n = cells
  ### transpose here
  jmat<-t(jmat)
  return(jmat)
}

# calculate expected JC matrix based on read depth: the expected JC index
# between two random cells i and j can be calculated as pi*pj/(pi+pj-pi*pj)
.normOVE <- function(p1, p2){
  pp = tcrossprod(p1, p2);
  ss = matrix(rep(p1,each=length(p2)), ncol=length(p2), byrow=TRUE) +  matrix(rep(p2, each=length(p1)), ncol=length(p2), byrow=FALSE)
  ee = pp/(ss - pp)
  ee<-t(ee)
  return(ee)
}

trainRegressModel <- function(rmat, p1, p2){ ### adopted from SnapATAC
  # remove the diag elements in the jmat
  idx.pairwise = which(rmat == 1, arr.ind=TRUE);
  # calculate expected JC matrix based on read depth
  emat = .normOVE(p1, p2);
  # estimate the global scaling factor
  scale.factor = mean(rmat / emat);
  # fill the missing value for the diagnoal elements
  rmat[idx.pairwise] = scale.factor * emat[idx.pairwise];
  data = data.frame(x=c(emat), y=c(rmat));
  # regression model to predict real JC matrix based on expected matrix. Real data indicates that the relationship is non-linear
  model <- lm(y ~ x + I(x^2), data);
  return(model);
}

normObservedJmat <- function(rmat, model, p1, p2, method){ ### adopted from SnapATAC

  idx.pairwise = which(rmat == 1, arr.ind=TRUE);
  emat = .normOVE(p1, p2);
  scale.factor = mean(rmat / emat);
  rmat[idx.pairwise] = scale.factor * emat[idx.pairwise];
  preds = predict(model, data.frame(x=c(emat)), se.fit = TRUE)
  if(method == "zscore"){
    norm = (c(rmat) - preds$fit) / (preds$se.fit);
  }else if(method == "residual"){
    norm = c(rmat) -  preds$fit;
  }
  nmat = matrix(norm, nrow(emat), ncol(emat));
  return(nmat);
}

### read transpoased matrix into R
read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  ## empty strings are converted to NA
  out = read.csv(text=x, sep=sep, header=header, na.strings = "", ...)
  return(out)
  
}
