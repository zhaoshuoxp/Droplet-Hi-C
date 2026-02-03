# Droplet Hi-C
Droplet Hi-C is a cheap, convenient and scalable method for chromatin architecture profiling in single cells, based on the widely available 10X Chromium Single Cell ATAC platform.

🍹 This repository contains scripts and notebook to reproduce the results for our [manuscript](https://www.nature.com/articles/s41587-024-02447-1): **Droplet Hi-C enables scalable, single-cell profiling of chromatin architecture in heterogeneous tissues**



## Modified Libraries & Environment

To support modern Python environments (Python 3.12+) and ensure stability with the latest scientific stack, we have included **custom-modified versions** of core libraries. These modifications address critical compatibility issues (e.g., NumPy 2.0 ABI changes, SciPy sparse matrix API updates).

### 1. schicluster (v0.1 modified)

We provide a patched version of `schicluster that is detached from git-scm versioning and explicitly adapted for:

- **NumPy >= 2.0.1**
- **SciPy >= 1.15.3**
- **pandas >= 2.3.3**
- **cooler >= 0.10.4**

**Installation:** Please install this modified version from the local source (do not install from PyPI):

https://github.com/zhaoshuoxp/scHiCluster/tree/master

### 2. hictools (Modified)

A patched version of `hictools` is provided to ensure compatibility .

**Installation:**

Bash

```
cd hictools_v2
g++ -O3 -o hictools hictools.cpp cxstring.cpp
```

> **Note:** We recommend installing these packages in **Editable Mode (`-e`)** inside a fresh conda environment to prevent conflicts with system-level packages.



## Abstract

Comprehensive analysis of chromatin architecture is crucial for understanding the gene regulatory programs during development and in disease pathogenesis, yet current methods often inadequately address the unique challenges presented by analysis of heterogeneous tissue samples. Here, we introduce Droplet Hi-C, which employs a commercial microfluidic device for high-throughput, single-cell chromatin conformation profiling in droplets. Using Droplet Hi-C we mapped the chromatin architecture at single-cell resolution from the mouse cortex and analyzed gene regulatory programs in major cortical cell types. Additionally, we used this technique to detect copy number variation (CNV), structural variations (SVs) and the extrachromosomal DNA (ecDNA) in cancer cells, revealing clonal dynamics and other oncogenic events during treatment. We further refined this technique to allow for joint profiling of chromatin architecture and transcriptome in single cells, facilitating a more comprehensive exploration of the links between chromatin architecture and gene expression in both normal tissues and tumors. Thus, Droplet Hi-C not only addresses critical gaps in chromatin analysis of heterogeneous tissues but also emerges as a versatile tool enhancing our understanding of gene regulation in health and disease.

![DHC_abstract](./images/abstract.png)


## [01.pre-process](https://github.com/Xieeeee/Droplet-Hi-C/tree/main/01.pre-process)
This directory contains scripts and analysis notebook for pre-processing (including mapping, contacts extraction, cells filtering) Droplet Hi-C data. 

## [02.analysis](https://github.com/Xieeeee/Droplet-Hi-C/tree/main/02.analysis)
This directory contains scripts and analysis notebook for re-producing figures in our manuscript. 

## 03.protocol
Protocol can be found at [protocol.io](https://www.protocols.io/view/droplet-hi-c-for-fast-and-scalable-profiling-of-ch-dpxe5pje).
