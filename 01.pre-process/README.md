## Preprocessing workflow
We used [scHiCluster](https://github.com/zhoujt1994/scHiCluster) to impute single cell contacts matrices, and performed domain calling as well as loop calling. In details, 
1. Fastq are processed, mapped and convert to pairs file using [04.proc_paired_hic_v2.1.sh](https://github.com/Xieeeee/Droplet-Hi-C/tree/main/01.pre-process/scripts/04.proc_paired_hic_v2.1.sh).
2. High quality cells are selected using the notebook [01.cell_filtering.ipynb](https://github.com/Xieeeee/Droplet-Hi-C/tree/main/01.pre-process/01.cell_filtering.ipynb).
3. Contacts information from individual cells are extracted using scripts [phc.sc_splitPairs_v2.py](https://github.com/Xieeeee/Droplet-Hi-C/tree/main/01.pre-process/scripts/phc.sc_splitPairs_v2.py).
4. Imputation is carried out for individual cells using [schic_impute_v2.sh](https://github.com/Xieeeee/Droplet-Hi-C/tree/main/01.pre-process/scripts/schic_impute_v2.sh).
5. Single cell GAD score is calculated using [schic_gene_score_v2.sh](https://github.com/Xieeeee/Droplet-Hi-C/tree/main/01.pre-process/scripts/schic_gene_score_v2.sh). This allows us to perform embedding and cell types annotation.
6. Single cell insulation score and domain calling are done using [schic_domain_v2.sh](https://github.com/Xieeeee/Droplet-Hi-C/tree/main/01.pre-process/scripts/schic_domain_v2.sh).
7. Loops calling for each cell cluster is done using [schic_loops_v3.sh](https://github.com/Xieeeee/Droplet-Hi-C/tree/main/01.pre-process/scripts/schic_loops_v3.sh).

For Paired Hi-C barcodes mapping between Hi-C and RNA modalities, matching relationship can be found in supp.
