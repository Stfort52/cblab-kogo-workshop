QC1-DropletUtils
================
2023-07-07

### **Library**

Here, we describe a short but to-the-point quality control process using
DropletUtils to remove unusable or low quality cells.

``` r
library(Seurat)
library(scater)
library(DropletUtils)
library(scran)
library(ggplot2)
```

### **Set directory path**

Set the path to load or save data.

``` r
rdatadir = '/BiO/home/data/QC/'
ranaldir = '/BiO/home/edu03/QC_1/'  # /home/username/QC_1/
```

make analysis directory if it does not exist
```r
if (!dir.exists(ranaldir)) { dir.create(ranaldir, recursive = TRUE)}
```

### **Load Data**

The fastq file obtained after scRNA-seq. can be mapped using CellRanger
to obtain a count matrix. Then you can load the data in R.

``` r
rawsce_1 <- read10xCounts(paste0(rdatadir, "20094_0001_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_2 <- read10xCounts(paste0(rdatadir, "20094_0002_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_3 <- read10xCounts(paste0(rdatadir, "20094_0003_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_4 <- read10xCounts(paste0(rdatadir, "20094_0004_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_5 <- read10xCounts(paste0(rdatadir, "20094_0005_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_6 <- read10xCounts(paste0(rdatadir, "20094_0006_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_7 <- read10xCounts(paste0(rdatadir, "20094_0007_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_8 <- read10xCounts(paste0(rdatadir, "20094_0008_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_9 <- read10xCounts(paste0(rdatadir, "20094_0009_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
rawsce_12 <- read10xCounts(paste0(rdatadir, "20094_0012_A_B/raw_feature_bc_matrix"), type = "sparse", compressed = TRUE)
```

### **Check Data**

``` r
rawsce_1
```

    ## class: SingleCellExperiment 
    ## dim: 36604 6794880 
    ## metadata(1): Samples
    ## assays(1): counts
    ## rownames(36604): ENSG00000243485 ENSG00000237613 ... Htag2 Htag3
    ## rowData names(3): ID Symbol Type
    ## colnames: NULL
    ## colData names(2): Sample Barcode
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):

### **Run DropletUtils**

To remove empty droplets, use DropletUtils for each data.

``` r
rawsce_list <- c('rawsce_1', 'rawsce_2', 'rawsce_3', 'rawsce_4', 'rawsce_5',
                 'rawsce_6', 'rawsce_7', 'rawsce_8', 'rawsce_9', 'rawsce_12')

for (i in 1:length(rawsce_list)) {
  rawsce <- get(rawsce_list[i])
  
  br.out <- barcodeRanks(counts(rawsce))
  
  p <- ggplot(data.frame(br.out), aes(x = rank, y = total)) + 
    geom_point() + 
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    ggtitle(rawsce_list[i]) +
    theme_classic()
  
  o <- order(br.out$rank)
  
  set.seed(2023)
  e.out <- emptyDrops(counts(rawsce))  ## Cells that have UMI counts lower than 100 are empty cells.
  table(Sig=e.out$FDR <= 0.05, Limited=e.out$Limited)
  is.cell <- e.out$FDR <= 0.05
  
  print(sum(is.cell, na.rm=TRUE))
  
  p <- p + geom_hline(yintercept=min(br.out$fitted[o], na.rm=TRUE), color="red", linetype="dashed")
  p <- p + geom_hline(yintercept=metadata(br.out)$knee, color="dodgerblue", linetype="dashed")
  p <- p + geom_hline(yintercept=min(metadata(br.out)$inflection), color="forestgreen", linetype="dashed")

  ggsave(filename = paste0(ranaldir, 'DropletUtils_', rawsce_list[i], '.png'), plot = p)
  
  colnames(rawsce) = colData(rawsce)$Barcode
  rawsce <- rawsce[,which(e.out$FDR <= 0.05)]
  
  assign(paste0('DropletUtils_', rawsce_list[i]), rawsce)
}
```

    ## [1] 12252
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 7280
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 6502
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 4487
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 4095
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 2659
    ## 
    ##   FALSE    TRUE 
    ## 6794877       3

    ## [1] 2001
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 5123
    ## 
    ##   FALSE    TRUE 
    ## 6794879       1

    ## [1] 2838
    ## 
    ##   FALSE 
    ## 6794880

    ## [1] 7178
    ## 
    ##   FALSE 
    ## 6794880

### **Check DropletUtils data**

``` r
DropletUtils_rawsce_1
```

    ## class: SingleCellExperiment 
    ## dim: 36604 12252 
    ## metadata(1): Samples
    ## assays(1): counts
    ## rownames(36604): ENSG00000243485 ENSG00000237613 ... Htag2 Htag3
    ## rowData names(3): ID Symbol Type
    ## colnames(12252): AAACCCAAGGCTGAAC-1 AAACCCAAGGGTGAAA-1 ...
    ##   TTTGTTGTCTTCGACC-1 TTTGTTGTCTTGGTCC-1
    ## colData names(2): Sample Barcode
    ## reducedDimNames(0):
    ## mainExpName: NULL
    ## altExpNames(0):

### **Reference**

Lun, A. T., Riesenfeld, S., Andrews, T., Gomes, T., & Marioni, J. C.
(2019). EmptyDrops: distinguishing cells from empty droplets in
droplet-based single-cell RNA sequencing data. Genome biology, 20(1),
1-9.

Pekayvaz, K., Leunig, A., Kaiser, R., Joppich, M., Brambs, S., Janjic,
A., … & Nicolai, L. (2022). Protective immune trajectories in early
viral containment of non-pneumonic SARS-CoV-2 infection. Nature
communications, 13(1), 1-21.
