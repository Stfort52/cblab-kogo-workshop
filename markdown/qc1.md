## QC1: DropletUtils

2026-07-21


## Remove empty droplets

Many droplets are empty, mainly consisted of ambient RNA.

We have to filter out these droplets.


1. Construct an ambient pool by summing up UMIs from low-count droplet
2. Test the possibility that ambient RNA accidentally generated this droplet profile
These process are wrapped up by package DropletUtils


```R
library(Seurat)
library(scater)
library(DropletUtils)
library(scran)
library(ggplot2)

set.seed(42) # for reproducibility


data_path <- if (dir.exists("data")) file.path("data", "qc1", "") else file.path("..", "data", "qc1", "")
save_path <- file.path("qc1/")

if (!dir.exists(save_path)) {
  dir.create(save_path)
}

```

### Load data


```R
raw_sce_list <- list()
sample_numbers <- c(1:9, 12)

for (i in sample_numbers) {
    file_path <- paste0(data_path, sprintf("20094_%04d_A_B/raw_feature_bc_matrix", i))
    raw_sce_list[[paste0("rawsce_", i)]] <- read10xCounts(file_path, type = "sparse", compressed = TRUE)
}
```

### Examine and check the data


```R
raw_sce_list$rawsce_1
```


    class: SingleCellExperiment 
    dim: 36604 6794880 
    metadata(1): Samples
    assays(1): counts
    rownames(36604): ENSG00000243485 ENSG00000237613 ... Htag2 Htag3
    rowData names(3): ID Symbol Type
    colnames: NULL
    colData names(2): Sample Barcode
    reducedDimNames(0):
    mainExpName: NULL
    altExpNames(0):



```R
filtered_sce <- list()

for (i in seq_along(sample_numbers)) {
  rawsce <- raw_sce_list[[i]]
  br.out <- barcodeRanks(counts(rawsce))
  
  e.out <- emptyDrops(counts(rawsce))  # Cells that have UMI counts lower than 100 (by default) are empty cells.
  print(table(Sig = e.out$FDR <= 0.05, Limited = e.out$Limited))

  is.cell <- e.out$FDR <= 0.05
  print(sum(is.cell, na.rm = TRUE))

  p <- ggplot(data.frame(br.out), aes(x = rank, y = total)) + 
    geom_point() + 
    scale_x_continuous(trans = "log10") +
    scale_y_continuous(trans = "log10") +
    ggtitle(paste0("raw_sce_", sample_numbers[i])) +
    theme_classic()
  p <- p + 
    geom_hline(aes(yintercept = min(br.out$fitted, na.rm = TRUE),  color = "FDR_0.05"), linetype = "dashed") +
    geom_hline(aes(yintercept = as.numeric(metadata(br.out)$knee), color = "knee"), linetype = "dashed") +
    geom_hline(aes(yintercept = as.numeric(metadata(br.out)$inflection), color = "inflection"), linetype = "dashed") +
    scale_color_manual(
        name = "Cutoffs", 
        values = c(
            "FDR_0.05" = "red",
            "knee" = "dodgerblue",
            "inflection"= "forestgreen"
        )
    ) +
    theme(legend.position = "bottom")
  ggsave(filename = paste0(save_path, "/DropletUtils_", sample_numbers[i], ".png"), plot = p, width = 6, height = 6)
  
  colnames(rawsce) <- colData(rawsce)$Barcode
  rawsce <- rawsce[,which(e.out$FDR <= 0.05)]
  
  filtered_sce[[paste0("DropletUtils_sce_", sample_numbers[i])]] <- rawsce
}
```

           Limited
    Sig     FALSE  TRUE
      FALSE 67355     0
      TRUE   3693  6467
    [1] 10160


    Warning message in min(br.out$fitted, na.rm = TRUE):
    “no non-missing arguments to min; returning Inf”
    Warning message in scale_y_continuous(trans = "log10"):
    “[1m[22m[32mlog-10[39m transformation introduced infinite values.”


           Limited
    Sig     FALSE  TRUE
      FALSE 74890     0
      TRUE    992  6574
    [1] 7566


    Warning message in min(br.out$fitted, na.rm = TRUE):
    “no non-missing arguments to min; returning Inf”
    Warning message in scale_y_continuous(trans = "log10"):
    “[1m[22m[32mlog-10[39m transformation introduced infinite values.”


           Limited
    Sig     FALSE  TRUE
      FALSE 68766     0
      TRUE    469  5881
    [1] 6350


    Warning message in min(br.out$fitted, na.rm = TRUE):
    “no non-missing arguments to min; returning Inf”
    Warning message in scale_y_continuous(trans = "log10"):
    “[1m[22m[32mlog-10[39m transformation introduced infinite values.”


           Limited
    Sig     FALSE  TRUE
      FALSE 30276     0
      TRUE    500  4650
    [1] 5150


    Warning message in min(br.out$fitted, na.rm = TRUE):
    “no non-missing arguments to min; returning Inf”
    Warning message in scale_y_continuous(trans = "log10"):
    “[1m[22m[32mlog-10[39m transformation introduced infinite values.”


           Limited
    Sig     FALSE TRUE
      FALSE  9764    0
      TRUE    536 4488
    [1] 5024


    Warning message in min(br.out$fitted, na.rm = TRUE):
    “no non-missing arguments to min; returning Inf”
    Warning message in scale_y_continuous(trans = "log10"):
    “[1m[22m[32mlog-10[39m transformation introduced infinite values.”


           Limited
    Sig     FALSE  TRUE
      FALSE 24018     0
      TRUE    408  3305
    [1] 3713


    Warning message in min(br.out$fitted, na.rm = TRUE):
    “no non-missing arguments to min; returning Inf”
    Warning message in scale_y_continuous(trans = "log10"):
    “[1m[22m[32mlog-10[39m transformation introduced infinite values.”


           Limited
    Sig     FALSE  TRUE
      FALSE 13523     0
      TRUE    263  2322
    [1] 2585


    Warning message in min(br.out$fitted, na.rm = TRUE):
    “no non-missing arguments to min; returning Inf”
    Warning message in scale_y_continuous(trans = "log10"):
    “[1m[22m[32mlog-10[39m transformation introduced infinite values.”


           Limited
    Sig     FALSE  TRUE
      FALSE 68688     0
      TRUE    373  4995
    [1] 5368


    Warning message in min(br.out$fitted, na.rm = TRUE):
    “no non-missing arguments to min; returning Inf”
    Warning message in scale_y_continuous(trans = "log10"):
    “[1m[22m[32mlog-10[39m transformation introduced infinite values.”


           Limited
    Sig     FALSE  TRUE
      FALSE 66869     0
      TRUE    515  2510
    [1] 3025


    Warning message in min(br.out$fitted, na.rm = TRUE):
    “no non-missing arguments to min; returning Inf”
    Warning message in scale_y_continuous(trans = "log10"):
    “[1m[22m[32mlog-10[39m transformation introduced infinite values.”


           Limited
    Sig     FALSE  TRUE
      FALSE 60423     0
      TRUE    948  6294
    [1] 7242


    Warning message in min(br.out$fitted, na.rm = TRUE):
    “no non-missing arguments to min; returning Inf”
    Warning message in scale_y_continuous(trans = "log10"):
    “[1m[22m[32mlog-10[39m transformation introduced infinite values.”


### Save the result data


```R
for (i in sample_numbers) {
    sce <- filtered_sce[[paste0("DropletUtils_sce_", i)]]
    save(sce, file = paste0(save_path, "/Dropletutils_filtered_sce_", i, ".RData"))
}
```

Reference
Lun, A. T., Riesenfeld, S., Andrews, T., Gomes, T., & Marioni, J. C. (2019). EmptyDrops: distinguishing cells from empty droplets in droplet-based single-cell RNA sequencing data. Genome biology, 20(1), 1-9.

Pekayvaz, K., Leunig, A., Kaiser, R., Joppich, M., Brambs, S., Janjic, A., ... & Nicolai, L. (2022). Protective immune trajectories in early viral containment of non-pneumonic SARS-CoV-2 infection. Nature communications, 13(1), 1-21.


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
     [1] scran_1.38.1                DropletUtils_1.30.0        
     [3] scater_1.38.1               ggplot2_4.0.3              
     [5] scuttle_1.20.0              SingleCellExperiment_1.32.0
     [7] SummarizedExperiment_1.40.0 Biobase_2.70.0             
     [9] GenomicRanges_1.62.1        Seqinfo_1.0.0              
    [11] IRanges_2.44.0              S4Vectors_0.48.1           
    [13] BiocGenerics_0.56.0         generics_0.1.4             
    [15] MatrixGenerics_1.22.0       matrixStats_1.5.0          
    [17] Seurat_5.5.0                SeuratObject_5.4.0         
    [19] sp_2.2-1                   
    
    loaded via a namespace (and not attached):
      [1] RcppAnnoy_0.0.23          splines_4.5.3            
      [3] later_1.4.8               pbdZMQ_0.3-14            
      [5] tibble_3.3.1              R.oo_1.27.1              
      [7] polyclip_1.10-7           fastDummies_1.7.6        
      [9] lifecycle_1.0.5           edgeR_4.8.2              
     [11] globals_0.19.1            lattice_0.22-9           
     [13] MASS_7.3-65               magrittr_2.0.5           
     [15] limma_3.66.0              plotly_4.12.0            
     [17] metapod_1.18.0            httpuv_1.6.17            
     [19] otel_0.2.0                sctransform_0.4.3        
     [21] spam_2.11-4               spatstat.sparse_3.2-0    
     [23] reticulate_1.46.0         cowplot_1.2.0            
     [25] pbapply_1.7-4             RColorBrewer_1.1-3       
     [27] abind_1.4-8               Rtsne_0.17               
     [29] purrr_1.2.2               R.utils_2.13.0           
     [31] ggrepel_0.9.8             irlba_2.3.7              
     [33] listenv_1.0.0             spatstat.utils_3.2-3     
     [35] goftest_1.2-3             RSpectra_0.16-2          
     [37] spatstat.random_3.4-5     dqrng_0.4.1              
     [39] fitdistrplus_1.2-6        parallelly_1.47.0        
     [41] DelayedMatrixStats_1.32.0 codetools_0.2-20         
     [43] DelayedArray_0.36.1       tidyselect_1.2.1         
     [45] farver_2.1.2              ScaledMatrix_1.18.0      
     [47] viridis_0.6.5             base64enc_0.1-6          
     [49] spatstat.explore_3.8-0    jsonlite_2.0.0           
     [51] BiocNeighbors_2.4.0       progressr_0.19.0         
     [53] ggridges_0.5.7            survival_3.8-6           
     [55] systemfonts_1.3.2         tools_4.5.3              
     [57] ragg_1.5.2                ica_1.0-3                
     [59] Rcpp_1.1.1-1.1            glue_1.8.1               
     [61] gridExtra_2.3             SparseArray_1.10.10      
     [63] IRdisplay_1.1             dplyr_1.2.1              
     [65] HDF5Array_1.38.0          withr_3.0.3              
     [67] fastmap_1.2.0             bluster_1.20.0           
     [69] rhdf5filters_1.22.0       digest_0.6.39            
     [71] rsvd_1.0.5                R6_2.6.1                 
     [73] mime_0.13                 textshaping_1.0.5        
     [75] scattermore_1.2           tensor_1.5.1             
     [77] spatstat.data_3.1-9       R.methodsS3_1.8.2        
     [79] h5mread_1.2.1             tidyr_1.3.2              
     [81] data.table_1.18.4         httr_1.4.8               
     [83] htmlwidgets_1.6.4         S4Arrays_1.10.1          
     [85] uwot_0.2.4                pkgconfig_2.0.3          
     [87] gtable_0.3.6              lmtest_0.9-40            
     [89] S7_0.2.2                  XVector_0.50.0           
     [91] htmltools_0.5.9           dotCall64_1.2            
     [93] scales_1.4.0              png_0.1-9                
     [95] spatstat.univar_3.2-0     reshape2_1.4.5           
     [97] uuid_1.2-2                nlme_3.1-169             
     [99] repr_1.1.7                zoo_1.8-15               
    [101] rhdf5_2.54.1              stringr_1.6.0            
    [103] KernSmooth_2.23-26        parallel_4.5.3           
    [105] miniUI_0.1.2              vipor_0.4.7              
    [107] pillar_1.11.1             grid_4.5.3               
    [109] vctrs_0.7.3               RANN_2.6.2               
    [111] promises_1.5.0            BiocSingular_1.26.1      
    [113] beachmat_2.26.0           xtable_1.8-8             
    [115] cluster_2.1.8.2           beeswarm_0.4.0           
    [117] evaluate_1.0.5            locfit_1.5-9.12          
    [119] cli_3.6.6                 compiler_4.5.3           
    [121] rlang_1.2.0               crayon_1.5.3             
    [123] future.apply_1.20.2       plyr_1.8.9               
    [125] ggbeeswarm_0.7.3          stringi_1.8.7            
    [127] viridisLite_0.4.3         deldir_2.0-4             
    [129] BiocParallel_1.44.0       lazyeval_0.2.3           
    [131] spatstat.geom_3.8-1       Matrix_1.7-5             
    [133] IRkernel_1.3.2            RcppHNSW_0.7.0           
    [135] patchwork_1.3.2           sparseMatrixStats_1.22.0 
    [137] future_1.70.0             Rhdf5lib_1.32.0          
    [139] statmod_1.5.2             shiny_1.14.0             
    [141] ROCR_1.0-12               igraph_2.3.2             

