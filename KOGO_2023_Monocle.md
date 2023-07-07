KOGO_Monocle
================
2023-06-26

### **Library**

We obtained SingleCellExperiment(SCE) objects without estimated droplets
from previous section. In this section, we will remove low quality cells
with large mitochondrial proportions, low RNA contents or low detected
gene numbers.

``` r
library(Seurat)
library(magrittr)
library(monocle)
library(SingleCellExperiment)
library(scran)
library(dplyr)
```

### **Data loading and pre-processing**

Here, we describe a brief trajectory analysis of T cell subset using
monocle. This dataset has various celltypes including T cell.

``` r
load(paste0(rdatadir, 'seurat.RData'))

p <- DimPlot(seurat, group.by = 'celltype',label=T)
ggsave(filename = paste0(rdatadir, 'dimplot.png'), plot = p)
```

Since we will draw trajectory graph of T cell, extract the T cell
population from whole dataset and re-normalize.

``` r
subset <- subset(seurat, cells = colnames(seurat)[seurat$celltype == 'T cell'])
SCEset <- as.SingleCellExperiment(subset, assay = 'RNA')

clusters <- quickCluster(SCEset)
SCEset <- computeSumFactors(SCEset, clusters = clusters)
SCEset <- logNormCounts(SCEset, pseudo_count = 1)
```

Highly variable genes are also changed specifically in T cell population
and it will be used for further analysis.

``` r
dec <- modelGeneVar(SCEset)
hvg.t <- getTopHVGs(dec, fdr.threshold = 0.05)

T_cell <- as.Seurat(SCEset)
```

### **Transfer into cds**

Monocle package uses differently structured object named cell_data_set
(cds), so normalized expressions, metadata for cells, and metadata for
genes shoud be recombined for creating cds object.

``` r
cds <- as.CellDataSet(T_cell)
```

### **Monocle process**

Next, we’ll perform some normalisation and variance estimation steps,
which will be used in the differential expression analyses later on.

``` r
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
```

Typically, to order cells by progress, we want to reduce the number of
genes analyzed. So we can select for a subset of genes that we believe
are important in setting said ordering, such as overdispersed genes. So
we can use the top 15 most overdispersed genes to order our cells.

``` r
cds <- reduceDimension(cds, num_dim = 15, reduction_method = 'tSNE')
cds <- monocle::clusterCells(cds)
```

    ## Distance cutoff calculated to 3.222216

``` r
p <- monocle::plot_cell_clusters(cds, color = "Cluster")
ggsave(filename = paste0(rdatadir, 'monocle_cluster.png'), plot = p)
```

### **Constructing Single Cell Trajectories**

In Monocle, a single cell trajectory is the inferred developmental
timeline of single cells.

The trajectory analysis consists of three stages: 1. Choose genes that
define progress 2. Reduce the dimensionality of the data 3. Order cells
in pseudotime

``` r
clustering_DEG_genes <- differentialGeneTest(cds, fullModelFormulaStr = '~Cluster')
ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)
cds <- reduceDimension(cds, num_dim = 15, method = 'DDRTree')

cds <- orderCells(cds)
p <- plot_cell_trajectory(cds, color_by = "Cluster")
ggsave(filename = paste0(rdatadir, 'monocle_Cluster.png'), plot = p)
p <- plot_cell_trajectory(cds, color_by = "Pseudotime")
ggsave(filename = paste0(rdatadir, 'monocle_Pseudotime.png'), plot = p)
```

### **Finding Genes that Change as a Function of Pseudotime**

Once we have a trajectory, we can use differentialGeneTest() to find
genes that have an expression pattern that varies according to
pseudotime.

``` r
pseudotime_de <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)")

pseudotime_de %>% arrange(qval) %>% head()
```

    ##          status           family          pval          qval gene_short_name
    ## IGKC         OK negbinomial.size  0.000000e+00  0.000000e+00            IGKC
    ## IGHG3        OK negbinomial.size 4.650128e-266 4.678958e-262           IGHG3
    ## STMN1        OK negbinomial.size 1.907198e-218 1.279348e-214           STMN1
    ## HMGB2        OK negbinomial.size 5.077172e-163 2.554325e-159           HMGB2
    ## MTRNR2L1     OK negbinomial.size 1.970460e-152 7.930707e-149        MTRNR2L1
    ## IGHG4        OK negbinomial.size 1.624966e-146 5.450136e-143           IGHG4
    ##          use_for_ordering
    ## IGKC                 TRUE
    ## IGHG3                TRUE
    ## STMN1                TRUE
    ## HMGB2                TRUE
    ## MTRNR2L1             TRUE
    ## IGHG4                TRUE

``` r
pseudotime_de %>% arrange(qval) %>% head(n = 100) %>% select(gene_short_name) -> pseudotime_gene
pseudotime_gene <- pseudotime_gene$gene_short_name

pseudotime_gene
```

    ##   [1] "IGKC"      "IGHG3"     "STMN1"     "HMGB2"     "MTRNR2L1"  "IGHG4"    
    ##   [7] "HBB"       "RPL41"     "MKI67"     "HIST1H4C"  "RPL36A"    "RPS29"    
    ##  [13] "RPL34"     "HLA-DRB1"  "RPL39"     "CXCL13"    "RPS14"     "GZMB"     
    ##  [19] "RPS12"     "RPS27"     "CENPF"     "HMGB1"     "RPL30"     "RPS28"    
    ##  [25] "KIAA0101"  "MT2A"      "RPL21"     "RPLP2"     "CCL3"      "RPL31"    
    ##  [31] "TUBB"      "RPS25"     "HLA-DRA"   "FABP5"     "IGLC3"     "NCL"      
    ##  [37] "RPL32"     "RPS27A"    "GAPDH"     "IGLC2"     "TYMS"      "EEF1A1"   
    ##  [43] "RPL36"     "NUSAP1"    "RPS10"     "RPL9"      "UBE2C"     "NKG7"     
    ##  [49] "IGHGP"     "RPL37"     "RPL10"     "RPS13"     "TOP2A"     "PTPRC"    
    ##  [55] "HMGN2"     "RPL23A"    "RPL37A"    "TPT1"      "H2AFZ"     "TUBA1B"   
    ##  [61] "RNF19A"    "RPS3A"     "RPL26"     "CD74"      "IGHG1"     "HBA1"     
    ##  [67] "RPS18"     "RPS26"     "RPL38"     "RPL35A"    "RPS21"     "ASPM"     
    ##  [73] "KLRB1"     "RPL27A"    "CREM"      "RPS15A"    "RPS23"     "RRM2"     
    ##  [79] "MT1E"      "RPL13"     "FOS"       "CDK1"      "RPL11"     "SMC4"     
    ##  [85] "NEAT1"     "CA1"       "RPS20"     "FOSB"      "RPL18A"    "HSPA6"    
    ##  [91] "RPS19"     "RPL27"     "PTTG1"     "RPS17"     "HBA2"      "RPL35"    
    ##  [97] "HLA-DRB5"  "HSPD1"     "MTRNR2L12" "RPS6"

We can see that the expression patterns of the genes change
significantly with pseudotime.

``` r
p <- plot_genes_in_pseudotime(cds[c('MKI67', 'TOP2A', 'STMN1'),])
ggsave(filename = paste0(rdatadir, 'monocle_proliferating.png'), plot = p)

p <- plot_genes_in_pseudotime(cds[c('CCR7', 'IL7R'),])
ggsave(filename = paste0(rdatadir, 'monocle_naive.png'), plot = p)

p <- plot_genes_in_pseudotime(cds[c('PDCD1', 'CTLA4', 'LAG3', 'TIGIT'),])
ggsave(filename = paste0(rdatadir, 'monocle_exhausted.png'), plot = p)
```

### **Reference**

Butler, A., Hoffman, P., Smibert, P., Papalexi, E. & Satija, R. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat. Biotechnol. 36, 411–420 (2018).

Cao, J. et al. The single-cell transcriptional landscape of mammalian organogenesis. Nature 566, 496–502 (2019).

Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018). ‘Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors.’ Nat. Biotechnol., 36(5), 421-427. doi: 10.1038/nbt.4091

L. Ma, M.O. Hernandez, Y. Zhao, M. Mehta, B. Tran, M. Kelly, Z. Rae, J.M. Hernandez, J.L. Davis, S.P. Martin, D.E. Kleiner, S.M. Hewitt, K. Ylaya, B.J. Wood, T.F. Greten, X.W. Wang. Tumor cell biodiversity drives microenvironmental reprogramming in liver cancer. Canc. Cell, 36 (4): 418-430 (2019) 

Lun, A. T., McCarthy, D. J. & Marioni, J. C. A step-by-step workflow for low-level analysis of single-cell RNA-seq data with Bioconductor. F1000Res 5, 2122 (2016). 

McCarthy, D. J., Campbell, K. R., Lun, A. T. & Wills, Q. F. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. Bioinformatics 33, 1179–1186 (2017)

Cole Trapnell and Davide Cacchiarelli et al (2014): The dynamics and regulators of cell fate decisions are revealed by pseudo-temporal ordering of single cells. Nature Biotechnology

Xiaojie Qiu, Andrew Hill, Cole Trapnell et al (2017):Single-cell mRNA quantification and differential analysis with Census. Nature Methods 

Xiaojie Qiu, Cole Trapnell et al (2017): Reverse graph embedding resolves complex single-cell developmental trajectories. BioRxiv
