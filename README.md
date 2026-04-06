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

> **Note:** We recommend installing these packages in **Editable Mode (`-e`)** inside a fresh conda environment to prevent conflicts with system-level packages.



### 2. hictools (C++)

A patched version of `hictools` is provided to ensure compatibility with NovaSeq X and improved performance .

- **Zero-Copy Parsing (C++17):** Replaced all `std::string` allocations (`substr`, `split`, concatenation) with  `std::string_view`.

- **Pipeline to `kseq.h`:** Bypassed the single-threaded `zlib` by piping `pigz` output directly into Heng Li's  `kseq.h` C library.

- **OpenMP Parallelism:**  a batched-processing model using `#pragma omp parallel for`. Now reads chunks of 100,000 SAM lines and distributes workload uniformly across all available logical CPU cores.

- **Double Buffering:** Implemented large, custom memory buffers (up to 4MB) to aggregate output data.

- **Smart Quality of Life:** Added automatic fallback detection for `.fq.gz` vs. `.fastq.gz` extensions and exposed dynamic thread count parameters to the CLI.

**To compile:**

```
cd hictools_v2
g++ -O3 -std=c++17 -fopenmp hictools.cpp -o hictools
```


### 3. hictools (Rust)

The Rust version provides a modernized, memory-safe alternative that delivers near C++ performance without the complexity of manual memory management or external C libraries.

**To compile Rust version hictools**:

1. install Rust:
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```
2. compile:
```
cd hictools_rs
cargo build --release
```
3. now hictools is in `target/release/hictools`

   


## Abstract

Comprehensive analysis of chromatin architecture is crucial for understanding the gene regulatory programs during development and in disease pathogenesis, yet current methods often inadequately address the unique challenges presented by analysis of heterogeneous tissue samples. Here, we introduce Droplet Hi-C, which employs a commercial microfluidic device for high-throughput, single-cell chromatin conformation profiling in droplets. Using Droplet Hi-C we mapped the chromatin architecture at single-cell resolution from the mouse cortex and analyzed gene regulatory programs in major cortical cell types. Additionally, we used this technique to detect copy number variation (CNV), structural variations (SVs) and the extrachromosomal DNA (ecDNA) in cancer cells, revealing clonal dynamics and other oncogenic events during treatment. We further refined this technique to allow for joint profiling of chromatin architecture and transcriptome in single cells, facilitating a more comprehensive exploration of the links between chromatin architecture and gene expression in both normal tissues and tumors. Thus, Droplet Hi-C not only addresses critical gaps in chromatin analysis of heterogeneous tissues but also emerges as a versatile tool enhancing our understanding of gene regulation in health and disease.

![DHC_abstract](./images/abstract.png)


## [01.pre-process](https://github.com/Xieeeee/Droplet-Hi-C/tree/main/01.pre-process)
This directory contains scripts and analysis notebook for pre-processing (including mapping, contacts extraction, cells filtering) Droplet Hi-C data. 

## [02.analysis](https://github.com/Xieeeee/Droplet-Hi-C/tree/main/02.analysis)
This directory contains scripts and analysis notebook for re-producing figures in our manuscript. 

## 03.protocol
Protocol can be found at [protocol.io](https://www.protocols.io/view/droplet-hi-c-for-fast-and-scalable-profiling-of-ch-dpxe5pje).
