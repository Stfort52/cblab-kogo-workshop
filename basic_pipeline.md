1. Normalization
2. Batch-aware HVG selection
3. Dimensional Reduction
4. Batch-Correciton
5. Annotation
6. DEG Analysis


```R
library(Seurat)
library(harmony)
library(ggplot2)
set.seed(42) # for reproducibility

save_path = "basic_pipeline"
setwd(save_path)
```

    Warning message:
    “package ‘Seurat’ was built under R version 4.3.2”
    Loading required package: SeuratObject
    
    Loading required package: sp
    
    Warning message:
    “package ‘sp’ was built under R version 4.3.3”
    
    Attaching package: ‘SeuratObject’
    
    
    The following object is masked from ‘package:base’:
    
        intersect
    
    
    Warning message:
    “package ‘harmony’ was built under R version 4.3.3”
    Loading required package: Rcpp
    
    Warning message:
    “package ‘Rcpp’ was built under R version 4.3.3”
    Warning message:
    “package ‘ggplot2’ was built under R version 4.3.3”


### 0. load data
data from https://www.nature.com/articles/s41591-023-02327-2<br>
You can access the data via CellxGene (https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293)<br>
In this example, we will use the sampled data


```R
count_matrix <- read.csv("/BiO/data/process/basic_pipeline_data/HLCA_pulmonary_fibrosis_immune_raw.csv", row.names = 1)
meta.data <- read.csv("/BiO/data/process/basic_pipeline_data/HLCA_pulmonary_fibrosis_immune_meta.csv", row.names = 1)

so <- CreateSeuratObject(counts = count_matrix, meta.data = meta.data, assay = 'RNA', min.cells = 0, min.features = 0, project = 'HLCA_Pulmonary_Fibrosis_immune')
# genes are in rows, cells are in columns

# so stand for 's'eurat 'o'bject
# In this example, we use filtered data, so set min.cells and min.features to 0 (no filtering)

# > head(so, n = 3)
#                                                      orig.ident nCount_RNA nFeature_RNA            disease                 study
# F01173_GCTGGGTTCCTGTAGA_haberman HLCA_Pulmonary_Fibrosis_immune       5525         1877 pulmonary fibrosis Banovich_Kropski_2020
# F00431_CTAGAGTCATGCCACG_haberman HLCA_Pulmonary_Fibrosis_immune       2784         1017 pulmonary fibrosis Banovich_Kropski_2020
# F01172_AGTAGTCGTCCGACGT_haberman HLCA_Pulmonary_Fibrosis_immune       1617         1012 pulmonary fibrosis Banovich_Kropski_2020
```

    Warning message:
    “Data is of class data.frame. Coercing to dgCMatrix.”


### 1. Normalization


```R
# seurat style
so <- NormalizeData(so)
```

    Normalizing layer: counts
    



```R
# scran style
# To remove cell-specific biases, cells are clustered using quickCluster() and cell-specific size factors are calculated using computeSumFactors() of scran R package. 
# Raw counts of each cell are divided by cell-specific size factor and log2-transformed with a pseudocount of 1.
library(scran)
library(scater)

sce <- as.SingleCellExperiment(so)

clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
sce.norm <- logNormCounts(sce,  pseudo_count = 1)
so <- as.Seurat(sce.norm,
                    counts = "counts",
                    data = "logcounts")
```

    Warning message:
    “package ‘scran’ was built under R version 4.3.3”
    Loading required package: SingleCellExperiment
    
    Warning message:
    “package ‘SingleCellExperiment’ was built under R version 4.3.2”
    Loading required package: SummarizedExperiment
    
    Warning message:
    “package ‘SummarizedExperiment’ was built under R version 4.3.2”
    Loading required package: MatrixGenerics
    
    Warning message:
    “package ‘MatrixGenerics’ was built under R version 4.3.3”
    Loading required package: matrixStats
    
    Warning message:
    “package ‘matrixStats’ was built under R version 4.3.3”
    
    Attaching package: ‘MatrixGenerics’
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
        colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
        colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
        colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
        colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
        colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
        colWeightedMeans, colWeightedMedians, colWeightedSds,
        colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
        rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
        rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
        rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
        rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
        rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
        rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
        rowWeightedSds, rowWeightedVars
    
    
    Loading required package: GenomicRanges
    
    Warning message:
    “package ‘GenomicRanges’ was built under R version 4.3.3”
    Loading required package: stats4
    
    Loading required package: BiocGenerics
    
    Warning message:
    “package ‘BiocGenerics’ was built under R version 4.3.2”
    
    Attaching package: ‘BiocGenerics’
    
    
    The following object is masked from ‘package:SeuratObject’:
    
        intersect
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
        match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
        Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
        table, tapply, union, unique, unsplit, which.max, which.min
    
    
    Loading required package: S4Vectors
    
    Warning message:
    “package ‘S4Vectors’ was built under R version 4.3.3”
    
    Attaching package: ‘S4Vectors’
    
    
    The following object is masked from ‘package:utils’:
    
        findMatches
    
    
    The following objects are masked from ‘package:base’:
    
        expand.grid, I, unname
    
    
    Loading required package: IRanges
    
    Warning message:
    “package ‘IRanges’ was built under R version 4.3.3”
    
    Attaching package: ‘IRanges’
    
    
    The following object is masked from ‘package:sp’:
    
        %over%
    
    
    Loading required package: GenomeInfoDb
    
    Warning message:
    “package ‘GenomeInfoDb’ was built under R version 4.3.2”
    Loading required package: Biobase
    
    Warning message:
    “package ‘Biobase’ was built under R version 4.3.3”
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    
    Attaching package: ‘Biobase’
    
    
    The following object is masked from ‘package:MatrixGenerics’:
    
        rowMedians
    
    
    The following objects are masked from ‘package:matrixStats’:
    
        anyMissing, rowMedians
    
    
    
    Attaching package: ‘SummarizedExperiment’
    
    
    The following object is masked from ‘package:Seurat’:
    
        Assays
    
    
    The following object is masked from ‘package:SeuratObject’:
    
        Assays
    
    
    Loading required package: scuttle
    
    Warning message:
    “package ‘scuttle’ was built under R version 4.3.3”
    Warning message:
    “package ‘scater’ was built under R version 4.3.2”
    Warning message:
    “Layer ‘scale.data’ is empty”


### 2. Batch-aware Feature selection
In this data, there are multiple studies. <br>
Therefore, we will select HVG with consideration of the source of data (study). <br>
batch-aware HVG selection is implemented in Scanpy python library <br>


```R
table(so$study)
```


    
     Banovich_Kropski_2020          Kaminski_2020 Misharin_Budinger_2018 
                      2338                   2328                    810 
             Sheppard_2020 
                      2306 



```R
# > table(so$study)
#   Budinger_2020 Lambrechts_2021  Wunderink_2021      Zhang_2021 
#            1041             833            2657            3439 
```


```R
source('/BiO/data/batch_aware_in_seurat.R')
```

    Loading required package: reticulate
    
    Warning message:
    “package ‘reticulate’ was built under R version 4.3.3”



```R
batch_key = 'study'
nHVG = 2000
batch_aware_in_seurat(so, batch_key = batch_key, nHVG = nHVG, conda_env = "/BiO/prog/miniforge3/envs/QC", save_path = save_path)
```

    Warning message:
    “No layers found matching search pattern provided”
    Warning message:
    “No layers found matching search pattern provided”
    Warning message:
    “Layer ‘scale.data’ is empty”
    Warning message:
    “Assay RNA changing from Assay5 to Assay”



```R
HVG = read.csv(paste0(save_path, "/hvg_", nHVG, "_", batch_key, ".csv"))
# > head(HVG)
#         gene
# 1        A2M
# 2       AACS
# 3 AARS.AARS1
# 4      ABCA1
# 5     ABCA13
# 6      ABCA2
```

### 3. Dimensionality Reduction
Normalized data is scaled and principal components (PCs) are calculated by a gene-by-cell matrix with HVGs.


```R
all.genes <- rownames(so)
so <- ScaleData(so, features = all.genes)

so <- RunPCA(so, features = HVG$gene)

p <- ElbowPlot(so, ndims = 50)
ggsave(p, filename = paste0(save_path, "/elbow_plot.png"), width = 4, height = 4)

PCs <- 8
so <- FindNeighbors(so, dims = 1:PCs)
so <- FindClusters(so, resolution = 0.5)

so <- RunUMAP(so, dims = 1:PCs)
# so <- RunTSNE(so, dims = 1:PCs)

DimPlot(so, group.by = 'study')

# save as png
p <- DimPlot(so, group.by = 'study')
ggsave(p, filename = paste0(save_path, "/umap_study.png"), width = 8, height = 6)
```

    Centering and scaling data matrix
    
    PC_ 1 
    Positive:  CST3, PSAP, CTSB, PLAUR, LYZ, C5AR1, S100A9, LST1, OLR1, CTSL 
    	   CTSZ, CD14, TIMP1, MS4A7, MRC1.MRC1L1, CFD, CD163, MNDA, C1QA, MARCO 
    	   VSIG4, LGALS3, MS4A6A, ETS2, CSTB, C1QB, NAMPT, S100A8, FBP1, SOD2.1.SOD2.1.SOD2.LOC100129518 
    Negative:  CD69, CD247, IL32, CST7, NKG7, GZMA, GZMB, SYNE2, PRF1, CD7 
    	   SKAP1, FYN, GNG2, KLRD1, CD3E, CD2, GNLY, CD96, CTSW, ETS1 
    	   ITK, PARP8, SYTL3, LCK, GZMH, CD3D, PRKCH, SLC38A1, AC243829.1.CCL4, STAT4 
    PC_ 2 
    Positive:  BTG1, DDIT4, CTSD, TUBA1B, HSP90AA1, SFTPC, CCL3, DNAJB1, CTSW, MS4A7 
    	   DUSP2, NKG7, GZMA, GIMAP7, LAT, CST7, CHMP1B, FKBP11, AC243829.1.CCL4, IL32 
    	   LYZ, CD3E, HOPX, GNLY, SH2D1A, GZMB, CD69, CD7, TRAF3IP3, GZMM 
    Negative:  SIPA1L1, ELL2, DENND4A, FOXO1, KDM6A, MT.ND4L, PITPNC1, SLC16A10, DENND5A, MGAT5 
    	   FMN1, NBPF19, KIF13B, FNDC3A, SH3PXD2B, B4GALT1, SPAG9, RGPD5, XYLT1, ATP13A3 
    	   PELI1, AC006978.2.GS1.114I9.3.ZNRF2.ENSG00000281593, SLC8A1, CTNNB1, SMYD3, USP36, NFKB1, AKT3, SYTL3, ABCA1 
    PC_ 3 
    Positive:  APOC1, ACP5, APOE, FABP4, CD9, GPNMB, C1QB, LPL, FABP5, C1QC 
    	   CYP27A1, C1QA, MARCO, NUPR1, FBP1, FN1, TREM2, SERPING1, LGALS3, CSTB 
    	   GLDN, MRC1.MRC1L1, SCD, MSR1, AC005027.3.INHBA, VSIG4, SLCO2B1, GSN, ABCG1, HSPB1 
    Negative:  FCN1, NAMPT, VCAN, APOBEC3A.B.APOBEC3A, FGL2, RGS2, CD300E, EREG, G0S2, TIMP1 
    	   S100A12, SOD2.1.SOD2.1.SOD2.LOC100129518, S100A8, TNFRSF1B, SLC25A37, NLRP3, IL1B, SLC2A3, AREGB.AREG, CD93 
    	   PLAUR, SERPINB9, THBS1, S100A9, C5AR1, LRRK2, NR4A1, ARL5B, CEBPD, SULF2 
    PC_ 4 
    Positive:  GNLY, KLRD1, NKG7, FGFBP2, KLRF1, PRF1, GZMB, LOC100130872.SPON2, CTSW, CLIC3 
    	   CMC1, FCGR3A, SH2D1B, CST7, S1PR5, ADGRG1.GPR56, C1ORF21.C1orf21, TTC38, PLAC8, AC243829.1.CCL4 
    	   CD247, TXK, CCL3, PRSS23, MCTP2, OSBPL5, FCRL6, PTGDR, YES1, HOPX 
    Negative:  IL7R, CD3D, LTB, ICOS, CTLA4, INPP4B, BATF, ADAM19, SPOCK2, PBX4 
    	   DUSP4, THEMIS, CD40LG, CD28, CD5, CD3G, CD6, CXCR6, SLAMF1, BIRC3 
    	   CD2, SIT1, RGS1, CXCR3, HNRNPLL, TRAT1, CD3E, DPP4, SUSD3, FAM129A.NIBAN1 
    PC_ 5 
    Positive:  HBB, NR4A1, HLA.DRB5, GBP5, FOS, ATF3, MT.ND4L, ITGA4, SCGB3A1, JUN 
    	   SLC9A3R1, GBP1, GLA, BHLHE40, WARS1.WARS, FOSB, HSPH1, CRIP1.AL928654.7.CRIP1.1.CRIP1.1.ENSG00000257341.AL928654.3, RAB24, PSAP 
    	   SLC20A1, ENSG00000234290.AC116366.6.IRF1, CDC42EP3, DDX3Y, CX3CR1, RABGEF1, APOL6, DNAJB1, TNFSF10, SCGB1A1 
    Negative:  FYN, SYTL3, PPP1R16B, PITPNC1, S100A12, FOXO1, XYLT1, PBX4, KLRC4.KLRK1, ANK3 
    	   THBS1, BACH2, ZNF831.RP3.492J12.2, ZEB1, BICDL1.CCDC64, DENND4A, AIM1.CRYBG1, GRK5, KLRB1, CASK 
    	   BCL11B, C1orf56.C1ORF56, MPP7, RNASE2, NCK2, CXCR4, ZBTB16, LUZP1, DUSP16, ENSG00000272980.RP11.517H2.6.Z94721.2 
    
    Computing nearest neighbor graph
    
    Warning message:
    “package ‘future’ was built under R version 4.3.3”
    Computing SNN
    


    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 7782
    Number of edges: 252571
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.9204
    Number of communities: 12
    Elapsed time: 0 seconds


    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”
    16:07:37 UMAP embedding parameters a = 0.9922 b = 1.112
    
    16:07:37 Read 7782 rows and found 8 numeric columns
    
    16:07:37 Using Annoy for neighbor search, n_neighbors = 30
    
    16:07:37 Building Annoy index with metric = cosine, n_trees = 50
    
    0%   10   20   30   40   50   60   70   80   90   100%
    
    [----|----|----|----|----|----|----|----|----|----|
    
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    |
    
    16:07:38 Writing NN index file to temp file /tmp/Rtmp13zLhx/filefd85032f6badc
    
    16:07:38 Searching Annoy index using 1 thread, search_k = 3000
    
    16:07:40 Annoy recall = 100%
    
    16:07:40 Commencing smooth kNN distance calibration using 1 thread
     with target n_neighbors = 30
    
    16:07:40 Initializing from normalized Laplacian + noise (using RSpectra)
    
    16:07:41 Commencing optimization for 500 epochs, with 307036 positive edges
    
    16:07:41 Using rng type: pcg
    
    16:07:49 Optimization finished
    



    
![png](output_14_3.png)
    



```R
ElbowPlot(so, ndims = 50)
```


    
![png](output_15_0.png)
    



```R
DimPlot(so, group.by = 'study')
```


    
![png](output_16_0.png)
    


### 4. Batch Correction by Harmony
ref) Korsunsky, I., Millard, N., Fan, J. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16, 1289–1296 (2019). https://doi.org/10.1038/s41592-019-0619-0 <br>
https://portals.broadinstitute.org/harmony/articles/quickstart.html


```R
so <- RunHarmony(so, 'study')
so <- FindNeighbors(so, reduction = "harmony")
so <- FindClusters(so, resolution = 0.5) 
so <- RunUMAP(so, dims = 1:PCs, reduction = 'harmony', reduction.name = 'umap.harmony') # use same dimension number as before
# so <- RunTSNE(so, dims = 1:PCs, reduction = 'harmony', reduction.name = 'tsne.harmony')

# save as png
p <- DimPlot(so, group.by = 'study', reduction = 'umap.harmony')
ggsave(p, filename = paste0(save_path, "/umap_harmony_study.png"), width = 8, height = 6)
p <- DimPlot(so, group.by = 'seurat_clusters', reduction = 'umap.harmony', label = TRUE)
ggsave(p, filename = paste0(save_path, "/umap_harmony_seurat_clusters.png"), width = 8, height = 6)
```

    Transposing data matrix
    
    Initializing state using k-means centroids initialization
    
    Harmony 1/10
    
    Harmony 2/10
    
    Harmony 3/10
    
    Harmony 4/10
    
    Harmony 5/10
    
    Harmony converged after 5 iterations
    
    Computing nearest neighbor graph
    
    Computing SNN
    


    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 7782
    Number of edges: 268158
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8761
    Number of communities: 9
    Elapsed time: 0 seconds


    16:07:57 UMAP embedding parameters a = 0.9922 b = 1.112
    
    16:07:57 Read 7782 rows and found 8 numeric columns
    
    16:07:57 Using Annoy for neighbor search, n_neighbors = 30
    
    16:07:57 Building Annoy index with metric = cosine, n_trees = 50
    
    0%   10   20   30   40   50   60   70   80   90   100%
    
    [----|----|----|----|----|----|----|----|----|----|
    
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    |
    
    16:07:58 Writing NN index file to temp file /tmp/Rtmp13zLhx/filefd8501ca0cce1
    
    16:07:58 Searching Annoy index using 1 thread, search_k = 3000
    
    16:08:00 Annoy recall = 100%
    
    16:08:00 Commencing smooth kNN distance calibration using 1 thread
     with target n_neighbors = 30
    
    16:08:01 Initializing from normalized Laplacian + noise (using RSpectra)
    
    16:08:01 Commencing optimization for 500 epochs, with 316290 positive edges
    
    16:08:01 Using rng type: pcg
    
    16:08:10 Optimization finished
    



```R
DimPlot(so, group.by = 'study', reduction = 'umap.harmony')
DimPlot(so, group.by = 'seurat_clusters', reduction = 'umap.harmony', label = TRUE)
```


    
![png](output_19_0.png)
    



    
![png](output_19_1.png)
    


### 5. Celltypist Prediction


```R
source("/BiO/data/celltypist_in_seurat.R")
celltypist_in_seurat(so, conda_env = "/BiO/prog/miniforge3/envs/QC", save_path = save_path, model_path = '/BiO/data/Immune_All_High.pkl')
```

    Warning message:
    “Assay RNA changing from Assay5 to Assay”


### 6. Marker gene expression visualization


```R
marker.genes<- list(T.cell = c('CD3D', 'CD3E'),
                    NK.cell = c('KLRD1', 'NKG7', 'GNLY'),
                    Macrophage = c('FABP4', 'APOE', 'MARCO'),
                    Monocyte = c('S100A12', 'FCN1'))
p <- FeaturePlot(so, features = as.vector(unlist(marker.genes)), ncol = 4, reduction = 'umap.harmony')
ggsave(p, filename = paste0(save_path, "/feature_plot_markers.png"), width = 13, height = 9)

p <- DotPlot(so, features = unlist(marker.genes), group.by = 'seurat_clusters')
p <- p + theme(axis.text.x = element_text(angle = 90))
ggsave(p, filename = paste0(save_path, "/dot_plot_markers.png"), width = 10, height = 6)

so$celltype = as.character(so$seurat_clusters)
so$celltype[so$celltype %in% c(1, 6, 8)] <- 'T cell'
so$celltype[so$celltype %in% c(3)] <- 'NK cell'
so$celltype[so$celltype %in% c(0, 4, 7)] <- 'Macrophage'
so$celltype[so$celltype %in% c(2, 5)] <- 'Monocyte'
so$celltype = factor(so$celltype, levels = c('T cell', 'NK cell', 'Macrophage', 'Monocyte'))

p <- DimPlot(so, group.by = 'celltype', reduction = 'umap.harmony', label = TRUE)
ggsave(p, filename = paste0(save_path, "/umap_harmony_celltype.png"), width = 8, height = 6)

p <- DotPlot(so, features = as.vector(unlist(marker.genes)), group.by = 'celltype')
p <- p + theme(axis.text.x = element_text(angle = 90))
ggsave(p, filename = paste0(save_path, "/dot_plot_celltype_markers.png"), width = 10, height = 6)

# saveRDS(so, file = paste0(save_path, "/HLCA_pulmonary_fibrosis_immune.rds"))
```

    Warning message:
    “[1m[22mThe `facets` argument of `facet_grid()` is deprecated as of ggplot2 2.2.0.
    [36mℹ[39m Please use the `rows` argument instead.
    [36mℹ[39m The deprecated feature was likely used in the [34mSeurat[39m package.
      Please report the issue at [3m[34m<https://github.com/satijalab/seurat/issues>[39m[23m.”
    Warning message:
    “Scaling data with a low number of groups may produce misleading results”



```R
FeaturePlot(so, features = as.vector(unlist(marker.genes)), ncol = 4, reduction = 'umap.harmony')
DotPlot(so, features = unlist(marker.genes), group.by = 'seurat_clusters') + theme(axis.text.x = element_text(angle = 90))
DimPlot(so, group.by = 'celltype', reduction = 'umap.harmony', label = TRUE)
DotPlot(so, features = as.vector(unlist(marker.genes)), group.by = 'celltype') + theme(axis.text.x = element_text(angle = 90))
```


    
![png](output_24_0.png)
    



    
![png](output_24_1.png)
    


    Warning message:
    “Scaling data with a low number of groups may produce misleading results”



    
![png](output_24_3.png)
    



    
![png](output_24_4.png)
    


### 7. Patient-aware DEG selection


```R
donor_id <- read.csv('/BiO/data/HLCA_pulmonary_fibrosis_donor_subset.csv', row.names = 1)
so$donor.id <- donor_id[Cells(so), ]
```


```R
library(EnhancedVolcano)

Tcells.control = Cells(so)[so$celltype == 'T cell' & so$disease == 'normal']
Tcells.pulmonary = Cells(so)[so$celltype == 'T cell' & so$disease == 'pulmonary fibrosis']

DEG.vanilla <- FindMarkers(so, ident.1 = Tcells.pulmonary, ident.2 = Tcells.control)
DEG.patient.correction <- FindMarkers(so, ident.1 = Tcells.pulmonary, ident.2 = Tcells.control, test.use = 'MAST', latent.vars = 'donor.id')

p <- EnhancedVolcano(DEG.vanilla,
    lab = rownames(DEG.vanilla),
    x = 'avg_log2FC',
    y = 'p_val_adj', 
    FCcutoff = log2(1.5),
    pCutoff = 0.05)
ggsave(p, filename = paste0(save_path, "/volcano_Tcells_vanilla.png"), width = 8, height = 6)

p <- EnhancedVolcano(DEG.patient.correction,
    lab = rownames(DEG.patient.correction),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    FCcutoff = log2(1.5),
    pCutoff = 0.05)
ggsave(p, filename = paste0(save_path, "/volcano_Tcells_patient_correction.png"), width = 8, height = 6)
```

    Warning message in .nextMethod(object = object, value = value):
    “Coefficients donor.idhomosapiens_None_2023_None_sikkemalisa_002_d10_1101_2022_03_10_483747VUILD65 are never estimible and will be dropped.”
    
     Completed [--------------------------------------------]   0% with 0 failures
    
     Completed [--------------------------------------------]   1% with 0 failures
    
     Completed [>-------------------------------------------]   1% with 0 failures
    
     Completed [>-------------------------------------------]   2% with 0 failures
    
     Completed [>-------------------------------------------]   3% with 0 failures
    
     Completed [=>------------------------------------------]   3% with 0 failures
    
     Completed [=>------------------------------------------]   4% with 0 failures
    
     Completed [=>------------------------------------------]   5% with 0 failures
    
     Completed [=>------------------------------------------]   6% with 0 failures
    
     Completed [==>-----------------------------------------]   6% with 0 failures
    
     Completed [==>-----------------------------------------]   7% with 0 failures
    
     Completed [==>-----------------------------------------]   8% with 0 failures
    
     Completed [===>----------------------------------------]   8% with 0 failures
    
     Completed [===>----------------------------------------]   9% with 0 failures
    
     Completed [===>----------------------------------------]  10% with 0 failures
    
     Completed [====>---------------------------------------]  10% with 0 failures
    
     Completed [====>---------------------------------------]  11% with 0 failures
    
     Completed [====>---------------------------------------]  12% with 0 failures
    
     Completed [=====>--------------------------------------]  13% with 0 failures
    
     Completed [=====>--------------------------------------]  14% with 0 failures
    
     Completed [=====>--------------------------------------]  15% with 0 failures
    
     Completed [======>-------------------------------------]  15% with 0 failures
    
     Completed [======>-------------------------------------]  16% with 0 failures
    
     Completed [======>-------------------------------------]  17% with 0 failures
    
     Completed [=======>------------------------------------]  17% with 0 failures
    
     Completed [=======>------------------------------------]  18% with 0 failures
    
     Completed [=======>------------------------------------]  19% with 0 failures
    
     Completed [========>-----------------------------------]  19% with 0 failures
    
     Completed [========>-----------------------------------]  20% with 0 failures
    
     Completed [========>-----------------------------------]  21% with 0 failures
    
     Completed [========>-----------------------------------]  22% with 0 failures
    
     Completed [=========>----------------------------------]  22% with 0 failures
    
     Completed [=========>----------------------------------]  23% with 0 failures
    
     Completed [=========>----------------------------------]  24% with 0 failures
    
     Completed [==========>---------------------------------]  24% with 0 failures
    
     Completed [==========>---------------------------------]  25% with 0 failures
    
     Completed [==========>---------------------------------]  26% with 0 failures
    
     Completed [===========>--------------------------------]  26% with 0 failures
    
     Completed [===========>--------------------------------]  27% with 0 failures
    
     Completed [===========>--------------------------------]  28% with 0 failures
    
     Completed [============>-------------------------------]  28% with 0 failures
    
     Completed [============>-------------------------------]  29% with 0 failures
    
     Completed [============>-------------------------------]  30% with 0 failures
    
     Completed [============>-------------------------------]  31% with 0 failures
    
     Completed [=============>------------------------------]  31% with 0 failures
    
     Completed [=============>------------------------------]  32% with 0 failures
    
     Completed [=============>------------------------------]  33% with 0 failures
    
     Completed [==============>-----------------------------]  33% with 0 failures
    
     Completed [==============>-----------------------------]  34% with 0 failures
    
     Completed [==============>-----------------------------]  35% with 0 failures
    
     Completed [===============>----------------------------]  35% with 0 failures
    
     Completed [===============>----------------------------]  36% with 0 failures
    
     Completed [===============>----------------------------]  37% with 0 failures
    
     Completed [================>---------------------------]  38% with 0 failures
    
     Completed [================>---------------------------]  39% with 0 failures
    
     Completed [================>---------------------------]  40% with 0 failures
    
     Completed [=================>--------------------------]  40% with 0 failures
    
     Completed [=================>--------------------------]  41% with 0 failures
    
     Completed [=================>--------------------------]  42% with 0 failures
    
     Completed [==================>-------------------------]  42% with 0 failures
    
     Completed [==================>-------------------------]  43% with 0 failures
    
     Completed [==================>-------------------------]  44% with 0 failures
    
     Completed [===================>------------------------]  44% with 0 failures
    
     Completed [===================>------------------------]  45% with 0 failures
    
     Completed [===================>------------------------]  46% with 0 failures
    
     Completed [===================>------------------------]  47% with 0 failures
    
     Completed [====================>-----------------------]  47% with 0 failures
    
     Completed [====================>-----------------------]  48% with 0 failures
    
     Completed [====================>-----------------------]  49% with 0 failures
    
     Completed [=====================>----------------------]  49% with 0 failures
    
     Completed [=====================>----------------------]  50% with 0 failures
    
     Completed [=====================>----------------------]  51% with 0 failures
    
     Completed [======================>---------------------]  51% with 0 failures
    
     Completed [======================>---------------------]  52% with 0 failures
    
     Completed [======================>---------------------]  53% with 0 failures
    
     Completed [=======================>--------------------]  53% with 0 failures
    
     Completed [=======================>--------------------]  54% with 0 failures
    
     Completed [=======================>--------------------]  55% with 0 failures
    
     Completed [=======================>--------------------]  56% with 0 failures
    
     Completed [========================>-------------------]  56% with 0 failures
    
     Completed [========================>-------------------]  57% with 0 failures
    
     Completed [========================>-------------------]  58% with 0 failures
    
     Completed [=========================>------------------]  58% with 0 failures
    
     Completed [=========================>------------------]  59% with 0 failures
    
     Completed [=========================>------------------]  60% with 0 failures
    
     Completed [==========================>-----------------]  60% with 0 failures
    
     Completed [==========================>-----------------]  61% with 0 failures
    
     Completed [==========================>-----------------]  62% with 0 failures
    
     Completed [===========================>----------------]  63% with 0 failures
    
     Completed [===========================>----------------]  64% with 0 failures
    
     Completed [===========================>----------------]  65% with 0 failures
    
     Completed [============================>---------------]  65% with 0 failures
    
     Completed [============================>---------------]  66% with 0 failures
    
     Completed [============================>---------------]  67% with 0 failures
    
     Completed [=============================>--------------]  67% with 0 failures
    
     Completed [=============================>--------------]  68% with 0 failures
    
     Completed [=============================>--------------]  69% with 0 failures
    
     Completed [==============================>-------------]  69% with 0 failures
    
     Completed [==============================>-------------]  70% with 0 failures
    
     Completed [==============================>-------------]  71% with 0 failures
    
     Completed [==============================>-------------]  72% with 0 failures
    
     Completed [===============================>------------]  72% with 0 failures
    
     Completed [===============================>------------]  73% with 0 failures
    
     Completed [===============================>------------]  74% with 0 failures
    
     Completed [================================>-----------]  74% with 0 failures
    
     Completed [================================>-----------]  75% with 0 failures
    
     Completed [================================>-----------]  76% with 0 failures
    
     Completed [=================================>----------]  76% with 0 failures
    
     Completed [=================================>----------]  77% with 0 failures
    
     Completed [=================================>----------]  78% with 0 failures
    
     Completed [==================================>---------]  78% with 0 failures
    
     Completed [==================================>---------]  79% with 0 failures
    
     Completed [==================================>---------]  80% with 0 failures
    
     Completed [==================================>---------]  81% with 0 failures
    
     Completed [===================================>--------]  81% with 0 failures
    
     Completed [===================================>--------]  82% with 0 failures
    
     Completed [===================================>--------]  83% with 0 failures
    
     Completed [====================================>-------]  83% with 0 failures
    
     Completed [====================================>-------]  84% with 0 failures
    
     Completed [====================================>-------]  85% with 0 failures
    
     Completed [=====================================>------]  85% with 0 failures
    
     Completed [=====================================>------]  86% with 0 failures
    
     Completed [=====================================>------]  87% with 0 failures
    
     Completed [======================================>-----]  88% with 0 failures
    
     Completed [======================================>-----]  89% with 0 failures
    
     Completed [======================================>-----]  90% with 0 failures
    
     Completed [=======================================>----]  90% with 0 failures
    
     Completed [=======================================>----]  91% with 0 failures
    
     Completed [=======================================>----]  92% with 0 failures
    
     Completed [========================================>---]  92% with 0 failures
    
     Completed [========================================>---]  93% with 0 failures
    
     Completed [========================================>---]  94% with 0 failures
    
     Completed [=========================================>--]  94% with 0 failures
    
     Completed [=========================================>--]  95% with 0 failures
    
     Completed [=========================================>--]  96% with 0 failures
    
     Completed [=========================================>--]  97% with 0 failures
    
     Completed [==========================================>-]  97% with 0 failures
    
     Completed [==========================================>-]  98% with 0 failures
    
     Completed [==========================================>-]  99% with 0 failures
    
     Completed [===========================================>]  99% with 0 failures
    
     Completed [===========================================>] 100% with 0 failures
    
     Completed [============================================] 100% with 0 failures
                                                                                  
    
    
    Done!
    
    Combining coefficients and standard errors
    
    Calculating log-fold changes
    
    Calculating likelihood ratio tests
    
    Refitting on reduced model...
    
    
     Completed [--------------------------------------------]   0% with 0 failures
    
     Completed [--------------------------------------------]   1% with 0 failures
    
     Completed [>-------------------------------------------]   1% with 0 failures
    
     Completed [>-------------------------------------------]   2% with 0 failures
    
     Completed [>-------------------------------------------]   3% with 0 failures
    
     Completed [=>------------------------------------------]   3% with 0 failures
    
     Completed [=>------------------------------------------]   4% with 0 failures
    
     Completed [=>------------------------------------------]   5% with 0 failures
    
     Completed [=>------------------------------------------]   6% with 0 failures
    
     Completed [==>-----------------------------------------]   6% with 0 failures
    
     Completed [==>-----------------------------------------]   7% with 0 failures
    
     Completed [==>-----------------------------------------]   8% with 0 failures
    
     Completed [===>----------------------------------------]   8% with 0 failures
    
     Completed [===>----------------------------------------]   9% with 0 failures
    
     Completed [===>----------------------------------------]  10% with 0 failures
    
     Completed [====>---------------------------------------]  10% with 0 failures
    
     Completed [====>---------------------------------------]  11% with 0 failures
    
     Completed [====>---------------------------------------]  12% with 0 failures
    
     Completed [=====>--------------------------------------]  13% with 0 failures
    
     Completed [=====>--------------------------------------]  14% with 0 failures
    
     Completed [=====>--------------------------------------]  15% with 0 failures
    
     Completed [======>-------------------------------------]  15% with 0 failures
    
     Completed [======>-------------------------------------]  16% with 0 failures
    
     Completed [======>-------------------------------------]  17% with 0 failures
    
     Completed [=======>------------------------------------]  17% with 0 failures
    
     Completed [=======>------------------------------------]  18% with 0 failures
    
     Completed [=======>------------------------------------]  19% with 0 failures
    
     Completed [========>-----------------------------------]  19% with 0 failures
    
     Completed [========>-----------------------------------]  20% with 0 failures
    
     Completed [========>-----------------------------------]  21% with 0 failures
    
     Completed [========>-----------------------------------]  22% with 0 failures
    
     Completed [=========>----------------------------------]  22% with 0 failures
    
     Completed [=========>----------------------------------]  23% with 0 failures
    
     Completed [=========>----------------------------------]  24% with 0 failures
    
     Completed [==========>---------------------------------]  24% with 0 failures
    
     Completed [==========>---------------------------------]  25% with 0 failures
    
     Completed [==========>---------------------------------]  26% with 0 failures
    
     Completed [===========>--------------------------------]  26% with 0 failures
    
     Completed [===========>--------------------------------]  27% with 0 failures
    
     Completed [===========>--------------------------------]  28% with 0 failures
    
     Completed [============>-------------------------------]  28% with 0 failures
    
     Completed [============>-------------------------------]  29% with 0 failures
    
     Completed [============>-------------------------------]  30% with 0 failures
    
     Completed [============>-------------------------------]  31% with 0 failures
    
     Completed [=============>------------------------------]  31% with 0 failures
    
     Completed [=============>------------------------------]  32% with 0 failures
    
     Completed [=============>------------------------------]  33% with 0 failures
    
     Completed [==============>-----------------------------]  33% with 0 failures
    
     Completed [==============>-----------------------------]  34% with 0 failures
    
     Completed [==============>-----------------------------]  35% with 0 failures
    
     Completed [===============>----------------------------]  35% with 0 failures
    
     Completed [===============>----------------------------]  36% with 0 failures
    
     Completed [===============>----------------------------]  37% with 0 failures
    
     Completed [================>---------------------------]  38% with 0 failures
    
     Completed [================>---------------------------]  39% with 0 failures
    
     Completed [================>---------------------------]  40% with 0 failures
    
     Completed [=================>--------------------------]  40% with 0 failures
    
     Completed [=================>--------------------------]  41% with 0 failures
    
     Completed [=================>--------------------------]  42% with 0 failures
    
     Completed [==================>-------------------------]  42% with 0 failures
    
     Completed [==================>-------------------------]  43% with 0 failures
    
     Completed [==================>-------------------------]  44% with 0 failures
    
     Completed [===================>------------------------]  44% with 0 failures
    
     Completed [===================>------------------------]  45% with 0 failures
    
     Completed [===================>------------------------]  46% with 0 failures
    
     Completed [===================>------------------------]  47% with 0 failures
    
     Completed [====================>-----------------------]  47% with 0 failures
    
     Completed [====================>-----------------------]  48% with 0 failures
    
     Completed [====================>-----------------------]  49% with 0 failures
    
     Completed [=====================>----------------------]  49% with 0 failures
    
     Completed [=====================>----------------------]  50% with 0 failures
    
     Completed [=====================>----------------------]  51% with 0 failures
    
     Completed [======================>---------------------]  51% with 0 failures
    
     Completed [======================>---------------------]  52% with 0 failures
    
     Completed [======================>---------------------]  53% with 0 failures
    
     Completed [=======================>--------------------]  53% with 0 failures
    
     Completed [=======================>--------------------]  54% with 0 failures
    
     Completed [=======================>--------------------]  55% with 0 failures
    
     Completed [=======================>--------------------]  56% with 0 failures
    
     Completed [========================>-------------------]  56% with 0 failures
    
     Completed [========================>-------------------]  57% with 0 failures
    
     Completed [========================>-------------------]  58% with 0 failures
    
     Completed [=========================>------------------]  58% with 0 failures
    
     Completed [=========================>------------------]  59% with 0 failures
    
     Completed [=========================>------------------]  60% with 0 failures
    
     Completed [==========================>-----------------]  60% with 0 failures
    
     Completed [==========================>-----------------]  61% with 0 failures
    
     Completed [==========================>-----------------]  62% with 0 failures
    
     Completed [===========================>----------------]  63% with 0 failures
    
     Completed [===========================>----------------]  64% with 0 failures
    
     Completed [===========================>----------------]  65% with 0 failures
    
     Completed [============================>---------------]  65% with 0 failures
    
     Completed [============================>---------------]  66% with 0 failures
    
     Completed [============================>---------------]  67% with 0 failures
    
     Completed [=============================>--------------]  67% with 0 failures
    
     Completed [=============================>--------------]  68% with 0 failures
    
     Completed [=============================>--------------]  69% with 0 failures
    
     Completed [==============================>-------------]  69% with 0 failures
    
     Completed [==============================>-------------]  70% with 0 failures
    
     Completed [==============================>-------------]  71% with 0 failures
    
     Completed [==============================>-------------]  72% with 0 failures
    
     Completed [===============================>------------]  72% with 0 failures
    
     Completed [===============================>------------]  73% with 0 failures
    
     Completed [===============================>------------]  74% with 0 failures
    
     Completed [================================>-----------]  74% with 0 failures
    
     Completed [================================>-----------]  75% with 0 failures
    
     Completed [================================>-----------]  76% with 0 failures
    
     Completed [=================================>----------]  76% with 0 failures
    
     Completed [=================================>----------]  77% with 0 failures
    
     Completed [=================================>----------]  78% with 0 failures
    
     Completed [==================================>---------]  78% with 0 failures
    
     Completed [==================================>---------]  79% with 0 failures
    
     Completed [==================================>---------]  80% with 0 failures
    
     Completed [==================================>---------]  81% with 0 failures
    
     Completed [===================================>--------]  81% with 0 failures
    
     Completed [===================================>--------]  82% with 0 failures
    
     Completed [===================================>--------]  83% with 0 failures
    
     Completed [====================================>-------]  83% with 0 failures
    
     Completed [====================================>-------]  84% with 0 failures
    
     Completed [====================================>-------]  85% with 0 failures
    
     Completed [=====================================>------]  85% with 0 failures
    
     Completed [=====================================>------]  86% with 0 failures
    
     Completed [=====================================>------]  87% with 0 failures
    
     Completed [======================================>-----]  88% with 0 failures
    
     Completed [======================================>-----]  89% with 0 failures
    
     Completed [======================================>-----]  90% with 0 failures
    
     Completed [=======================================>----]  90% with 0 failures
    
     Completed [=======================================>----]  91% with 0 failures
    
     Completed [=======================================>----]  92% with 0 failures
    
     Completed [========================================>---]  92% with 0 failures
    
     Completed [========================================>---]  93% with 0 failures
    
     Completed [========================================>---]  94% with 0 failures
    
     Completed [=========================================>--]  94% with 0 failures
    
     Completed [=========================================>--]  95% with 0 failures
    
     Completed [=========================================>--]  96% with 0 failures
    
     Completed [=========================================>--]  97% with 0 failures
    
     Completed [==========================================>-]  97% with 0 failures
    
     Completed [==========================================>-]  98% with 0 failures
    
     Completed [==========================================>-]  99% with 0 failures
    
     Completed [===========================================>]  99% with 0 failures
    
     Completed [===========================================>] 100% with 0 failures
    
     Completed [============================================] 100% with 0 failures
                                                                                  
    
    
    Done!
    

