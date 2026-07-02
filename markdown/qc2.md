# QC2: Remove low quality cells

2026-07-21


```R
options(repr.plot.format = "png", jupyter.plot_mimetypes = c("image/png", "image/jpeg"))
```

## Remove low-quality cells

Removing low-quality cells is crucial for downstream analysis steps.
Each dataset has its own QC thresholds, and these thresholds should be determined based on reasonable evidence.
If you encounter any low-quality cells during downstream analysis, re-run the QC with adjusted thresholds!

We obtained SingleCellExperiment (SCE) objects without estimated droplets from the previous section.
In this section, we will remove low-quality cells that have large mitochondrial proportions, low RNA content, or low numbers of detected genes.


```R
library(DropletUtils)
library(dplyr)
library(scater)

set.seed(42) # for reproducibility

data_path <- file.path("data/", "qc2/")
save_path <- file.path("qc2")

if (!dir.exists(save_path)) {
  dir.create(save_path)
}

```

### Load droplet-filtered SCE objects


```R
sample_numbers <- c(1:9, 12)
sce_list <- list()

for (i in sample_numbers) {
    file_name <- paste0(data_path, "DropletUtils_filtered_sce_", i, ".RData")
    tmpenv <- new.env()
    load(file_name, envir = tmpenv)
    sce_list[[paste0("sce_", i)]] <- tmpenv[[paste0("DropletUtils_rawsce_", i)]]
}

```

### Process a single sample as an example

A SCE object contains:
- a gene-by-cell counts as a sparse matrix
- gene data (gene annotation, etc)
- cell data (sample information, experimental condition information, etc). 

Gene and cell informations are stored in `rowData(SCE)` and `colData(SCE)`, respectively. Here, Ensembl gene ids are translated into gene symbol for easier further analysis.

Since there are 10 samples, we will process an example sample step-by-step and make a function `preprocess()` to apply the same process to other samples.

The example SCE object is shown below.


```R
sce <- sce_list$sce_5 # 5 for example
sce
```


    class: SingleCellExperiment 
    dim: 36604 4096 
    metadata(1): Samples
    assays(1): counts
    rownames(36604): ENSG00000243485 ENSG00000237613 ... Htag2 Htag3
    rowData names(3): ID Symbol Type
    colnames(4096): AAACCCACATCCTAAG-1 AAACCCATCAAGCTGT-1 ...
      TTTGTTGGTGGGATTG-1 TTTGTTGTCTTCGACC-1
    colData names(2): Sample Barcode
    reducedDimNames(0):
    mainExpName: NULL
    altExpNames(0):


To remove low quality cells, we use the `addPerCellQC()` function from the `scater` R package. This function calculates several quality control (QC) metrics for each cell, including:
- number of unique molecular identifiers (UMIs) per cell
- number of genes detected per cell
- percentage of UMIs assigned to mitochondrial (MT) genes


```R
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

mtgenes <- rowData(sce)[grep("^MT-", rowData(sce)$Symbol),]$Symbol
is.mito <- rownames(sce) %in% mtgenes
print(table(is.mito))
```

    is.mito
    FALSE  TRUE 
    36591    13 



```R
sce <- addPerCellQC(
  sce,
  subsets = list(MT=mtgenes),
  percent_top = c(50, 100, 200, 500), 
  detection_limit = 0
)

sce$log10_sum <- log10(sce$sum + 1)
sce$log10_detected <- log10(sce$detected + 1)

sce <- sce[,sce$sum!=0]
```

We define low-quality cells with:
-  &lt;500 UMIs
-  &gt;15% MT gene percents
-  &lt;100 detected genes

Criteria is visualized as histogram below.


```R
umi <- 500
mtpct <- 15
detect <- 100

ggplot(data.frame(colData(sce)), aes(x = log10_sum)) + 
  geom_histogram(bins = 150) +
  theme_bw() +
  geom_vline(xintercept = log10(umi), color="red", linetype="dashed")
  
ggplot(data.frame(colData(sce)), aes(x = subsets_MT_percent)) + 
  geom_histogram(bins = 150) +
  theme_bw() +
  geom_vline(xintercept = mtpct, color="red", linetype="dashed")

ggplot(data.frame(colData(sce)), aes(x = detected)) + 
  geom_histogram(bins = 150) +
  theme_bw() +
  geom_vline(xintercept = detect, color="red", linetype="dashed") + 
  ylim(0, 500)

```


    
![png](qc2_files/qc2_14_0.png)
    


    Warning message:
    “[1m[22mRemoved 1 row containing missing values or values outside the scale range
    (`geom_bar()`).”



    
![png](qc2_files/qc2_14_2.png)
    



    
![png](qc2_files/qc2_14_3.png)
    


Low-quality cells being filtered out can be identified by the PCA plot.


```R
filter_by_total_counts <- sce$sum > umi
filter_by_mt_percent <- sce$subsets_MT_percent < mtpct
filter_by_nfeature <- sce$detected > detect

sce <- runColDataPCA(sce, variables = list("sum", "detected", "subsets_MT_percent", "percent.top_500"))

sce$use <- (
  filter_by_total_counts &
    filter_by_mt_percent &
    filter_by_nfeature
)

plotReducedDim(sce, dimred="PCA_coldata", colour_by="use") +
  theme(plot.margin = margin(10, 10, 10, 10))

```


    
![png](qc2_files/qc2_17_0.png)
    



```R
sce <- sce[,sce$use]
sce
```


    class: SingleCellExperiment 
    dim: 36604 1985 
    metadata(1): Samples
    assays(1): counts
    rownames(36604): MIR1302-2HG FAM138A ... Htag2 Htag3
    rowData names(4): ID Symbol Type subsets_MT
    colnames(1985): AAACCCATCTGAACGT-1 AAACGAACAAGAGAGA-1 ...
      TTTGTTGAGTATAGAC-1 TTTGTTGGTGGGATTG-1
    colData names(15): Sample Barcode ... log10_detected use
    reducedDimNames(1): PCA_coldata
    mainExpName: NULL
    altExpNames(0):


### Process other samples
Other samples can be processed in the same way, and it's baked into a function `preprocess()` below.


```R
preprocess <- function(sce){
  rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
  
  mtgenes <- rowData(sce)[grep("^MT-", rowData(sce)$Symbol),]$Symbol
  is.mito <- rownames(sce) %in% mtgenes
  
  sce <- addPerCellQC(
    sce,
    subsets = list(MT=mtgenes),
    percent_top = c(50, 100, 200, 500), 
    detection_limit = 5
  )
  
  sce$log10_sum <- log10(sce$sum + 1)
  sce$log10_detected <- log10(sce$detected + 1)
  
  sce <- sce[, sce$sum!=0]
  return(sce)
}

filtering <- function(sce, umi, mtpct, detect){
  filter_by_total_counts <- sce$sum > umi
  filter_by_mt_percent <- sce$subsets_MT_percent < mtpct
  filter_by_nfeature <- sce$detected > detect
  
  sce <- runColDataPCA(sce, variables = list("sum", "detected", "subsets_MT_percent", "percent.top_500"))
  
  sce$use <- (
    filter_by_total_counts &
      filter_by_mt_percent &
      filter_by_nfeature
  )
  
  sce <- sce[, sce$use]
  return(sce)
}
```


```R
umi <- 500
mtpct <- 15
detect <- 100

preprocessed_sce_list <- lapply(sce_list, preprocess)
filtered_sce_list <- lapply(preprocessed_sce_list, filtering, umi, mtpct, detect)
```

### Change cellbarcode name
Each cell has their own cellbarcode. To prevent errors due to the presence of the duplicate cellbarcode across samples, we will append the sample information to the cellbarcode. 

Currently, cellbarcode from the first sample looks like below.


```R
print(head(colnames(filtered_sce_list$sce_1)))
```

    [1] "AAACCCAAGGGTGAAA-1" "AAACCCATCAGTCTTT-1" "AAACGAAGTGCGGATA-1"
    [4] "AAACGAATCGTAACTG-1" "AAACGCTGTCTCACGG-1" "AAAGGATTCGTGTCAA-1"



```R
for (i in seq_along(sample_numbers)) {
    new_names <- paste0(substring(colnames(filtered_sce_list[[i]]), 1, 16), "-L", sample_numbers[i])
    colnames(filtered_sce_list[[i]]) <- new_names
}

```

After appending the sample information, the cellbarcode will look like below.


```R
print(head(colnames(filtered_sce_list$sce_1)))
```

    [1] "AAACCCAAGGGTGAAA-L1" "AAACCCATCAGTCTTT-L1" "AAACGAAGTGCGGATA-L1"
    [4] "AAACGAATCGTAACTG-L1" "AAACGCTGTCTCACGG-L1" "AAAGGATTCGTGTCAA-L1"


### Adding sample information as metadata

For downstream analysis, we will combine all SCE objects add metadata about their conditions (`CONTROL`, `ASYMPTOMATIC` or `SYMPTOMATIC`)




```R
for (i in seq_along(sample_numbers)) {
    filtered_sce_list[[i]]$library <- paste0("L", sample_numbers[i])
}

sce <- do.call(cbind, filtered_sce_list)

table(sce$library)
```


    
      L1  L12   L2   L3   L4   L5   L6   L7   L8   L9 
    1823 3293 2166 1835  833 1350  324  419 1345  678 



```R
sce$Condition <- "SYMPTOMATIC"
sce$Condition[sce$library == "L12"] <- "CONTROL"
sce$Condition[sce$library %in% c("L4","L5","L6")] <- "ASYMPTOMATIC"

table(sce$Condition)
```


    
    ASYMPTOMATIC      CONTROL  SYMPTOMATIC 
            2507         3293         8266 


## References

Butler, A., Hoffman, P., Smibert, P., Papalexi, E. & Satija, R. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat. Biotechnol. 36, 411-420 (2018).

McCarthy, D. J., Campbell, K. R., Lun, A. T. & Wills, Q. F. Scater: pre-processing, quality control, normalization and visualization of single-cell RNA-seq data in R. Bioinformatics 33, 1179-1186 (2017).

Pekayvaz, K., Leunig, A., Kaiser, R. et al. Protective immune trajectories in early viral containment of non-pneumonic SARS-CoV-2 infection. Nat Commun 13, 1018 (2022).


    R version 4.5.3 (2026-03-11)
    Platform: x86_64-conda-linux-gnu
    Running under: Ubuntu 24.04.4 LTS
    
    Matrix products: default
    BLAS/LAPACK: /opt/conda/lib/libopenblasp-r0.3.33.so;  LAPACK version 3.12.0
    
    locale:
     [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
     [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
     [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    
    time zone: Etc/UTC
    tzcode source: system (glibc)
    
    attached base packages:
    [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    [8] base     
    
    other attached packages:
     [1] scater_1.38.1               ggplot2_4.0.3              
     [3] scuttle_1.20.0              dplyr_1.2.1                
     [5] DropletUtils_1.30.0         SingleCellExperiment_1.32.0
     [7] SummarizedExperiment_1.40.0 Biobase_2.70.0             
     [9] GenomicRanges_1.62.1        Seqinfo_1.0.0              
    [11] IRanges_2.44.0              S4Vectors_0.48.1           
    [13] BiocGenerics_0.56.0         generics_0.1.4             
    [15] MatrixGenerics_1.22.0       matrixStats_1.5.0          
    
    loaded via a namespace (and not attached):
     [1] tidyselect_1.2.1          viridisLite_0.4.3        
     [3] IRdisplay_1.1             vipor_0.4.7              
     [5] farver_2.1.2              viridis_0.6.5            
     [7] R.utils_2.13.0            S7_0.2.2                 
     [9] fastmap_1.2.0             digest_0.6.39            
    [11] rsvd_1.0.5                lifecycle_1.0.5          
    [13] Cairo_1.7-0               statmod_1.5.2            
    [15] magrittr_2.0.5            compiler_4.5.3           
    [17] rlang_1.2.0               tools_4.5.3              
    [19] labeling_0.4.3            S4Arrays_1.10.1          
    [21] dqrng_0.4.1               DelayedArray_0.36.1      
    [23] repr_1.1.7                RColorBrewer_1.1-3       
    [25] abind_1.4-8               BiocParallel_1.44.0      
    [27] HDF5Array_1.38.0          pbdZMQ_0.3-14            
    [29] withr_3.0.3               R.oo_1.27.1              
    [31] grid_4.5.3                beachmat_2.26.0          
    [33] Rhdf5lib_1.32.0           edgeR_4.8.2              
    [35] scales_1.4.0              cli_3.6.6                
    [37] crayon_1.5.3              ragg_1.5.2               
    [39] DelayedMatrixStats_1.32.0 ggbeeswarm_0.7.3         
    [41] rhdf5_2.54.1              parallel_4.5.3           
    [43] XVector_0.50.0            base64enc_0.1-6          
    [45] vctrs_0.7.3               Matrix_1.7-5             
    [47] jsonlite_2.0.0            BiocSingular_1.26.1      
    [49] BiocNeighbors_2.4.0       ggrepel_0.9.8            
    [51] irlba_2.3.7               beeswarm_0.4.0           
    [53] systemfonts_1.3.2         h5mread_1.2.1            
    [55] locfit_1.5-9.12           limma_3.66.0             
    [57] glue_1.8.1                codetools_0.2-20         
    [59] cowplot_1.2.0             gtable_0.3.6             
    [61] ScaledMatrix_1.18.0       tibble_3.3.1             
    [63] pillar_1.11.1             htmltools_0.5.9          
    [65] rhdf5filters_1.22.0       IRkernel_1.3.2           
    [67] R6_2.6.1                  textshaping_1.0.5        
    [69] sparseMatrixStats_1.22.0  evaluate_1.0.5           
    [71] lattice_0.22-9            R.methodsS3_1.8.2        
    [73] Rcpp_1.1.1-1.1            uuid_1.2-2               
    [75] gridExtra_2.3             SparseArray_1.10.10      
    [77] pkgconfig_2.0.3          

