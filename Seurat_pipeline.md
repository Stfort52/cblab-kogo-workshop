Seurat_pipeline
================
2023-07-07

### **Load pacakage and dataset**

SingleCellExperiment object which quality control (QC) was done should
be loaded.

``` r
library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
```

``` r
load('/home/data/sce.RData')
sce
```

    ## class: SingleCellExperiment 
    ## dim: 36604 14066 
    ## metadata(10): Samples Samples ... Samples Samples
    ## assays(1): counts
    ## rownames(36604): MIR1302-2HG FAM138A ... Htag2 Htag3
    ## rowData names(3): ID Symbol Type
    ## colnames(14066): AAACCCAAGGGTGAAA-L1 AAACCCATCAGTCTTT-L1 ...
    ##   TTTGTTGGTGGTTTAC-L12 TTTGTTGTCACTCACC-L12
    ## colData names(17): Sample Barcode ... library Condition
    ## reducedDimNames(1): PCA_coldata
    ## mainExpName: NULL
    ## altExpNames(0):

### **Create Seurat object**

``` r
seurat <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
seurat
```

    ## An object of class Seurat 
    ## 36604 features across 24389 samples within 1 assay 
    ## Active assay: RNA (36604 features, 0 variable features)

### **Normalization**

``` r
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize")
```

### **Feature selection (calculate highly variable genes)**

To find genes contain useful information about the biology of the data,
highly variable genes (HVGs) are defined by selecting the most variable
genes based on their expression across cells.

``` r
seurat <- FindVariableFeatures(seurat, nfeatures = 1000)
```

### **Dimensionality reduction**

Normalized data is scaled and principal components (PCs) are calculated
by a gene-by-cell matrix with HVGs.

``` r
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)

seurat <- RunPCA(seurat, features = VariableFeatures(seurat))
```

``` r
p <- ElbowPlot(seurat, ndims = 50)
ggsave(filename = "KOGO_Seurat/elbow_plot.png", plot = p)
```

We set to 10 PCs for clustering and visualization. After clustering and
visualization, cells are plotted in the two-dimensional TSNE or UMAP
plot and cell information can be also shown.

``` r
PCs <- 10

seurat <- FindNeighbors(seurat,
                        dims = 1:PCs)
seurat <- FindClusters(seurat,
                       resolution = 0.5)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 14066
    ## Number of edges: 777229
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.9098
    ## Number of communities: 17
    ## Elapsed time: 1 seconds

``` r
seurat <- RunUMAP(seurat, dims = 1:PCs)
seurat <- RunTSNE(seurat, dims = 1:PCs)
```

``` r
p <- DimPlot(seurat, label = TRUE)
ggsave(filename = "KOGO_Seurat/umap.png", plot = p)
```

### **Canonical cell type marker gene expression**

Based on ‘seurat_clusters’ after clustering, cells are grouped according
to their cell types as annotated based on known cell lineage-specific
marker genes

``` r
marker.genes<- list(T.cell = c('CD3D','CD3E','TRAC'),
                    Monocyte = c('CD14','S100A8','FCGR3A'),
                    NK.cell = c('NCAM1','NKG7','TRDC'),
                    B.cell = c('CD79A','CD83','MS4A1'),
                    Classical.Dendritic = c('CLEC9A','FCER1A','CLEC10A'),
                    Plasmacytoid.Dendritic = c('CLEC4C','IL3RA','LILRA4'),
                    Plasma.cell = c('JCHAIN','TNFRSF17','SDC1'))
```

## T.cell

``` r
p <- FeaturePlot(seurat, features = marker.genes$T.cell, order = T, ncol = 3)
ggsave(filename = "KOGO_Seurat/marker_T.png", plot = p, width = 10, height = 3)
```

## Monocyte

``` r
p <- FeaturePlot(seurat, features = marker.genes$Monocyte, order = T, ncol = 3)
ggsave(filename = "KOGO_Seurat/marker_monocyte.png", plot = p, width = 10, height = 3)
```

## NK.cell

``` r
p <- FeaturePlot(seurat, features = marker.genes$NK.cell, order = T, ncol = 3)
ggsave(filename = "KOGO_Seurat/marker_NK.png", plot = p, width = 10, height = 3)
```

## B.cell

``` r
p <- FeaturePlot(seurat, features = marker.genes$B.cell, order = T, ncol = 3)
ggsave(filename = "KOGO_Seurat/marker_B.png", plot = p, width = 10, height = 3)
```

## Classical.Dendritic

``` r
p <- FeaturePlot(seurat, features = marker.genes$Classical.Dendritic, order = T, ncol = 3)
ggsave(filename = "KOGO_Seurat/marker_classical-dendritic.png", plot = p, width = 10, height = 3)
```

## Plasmacytoid.Dendritic

``` r
p <- FeaturePlot(seurat, features = marker.genes$Plasmacytoid.Dendritic, order = T, ncol = 3)
ggsave(filename = "KOGO_Seurat/marker_plasmacytoid-dendritic.png", plot = p, width = 10, height = 3)
```

## Plasma.cell

``` r
p <- FeaturePlot(seurat, features = marker.genes$Plasma.cell, order = T, ncol = 3)
ggsave(filename = "KOGO_Seurat/marker_plasma.png", plot = p, width = 10, height = 3)
```

### **Marker gene expression on Heatmap**

``` r
library(pheatmap)
library(RColorBrewer)

avgExprs <- AverageExpression(seurat,
                              features = unlist(marker.genes),
                              assays = "RNA", slot = "data")

scaledExprs <- t(scale(t(avgExprs$RNA)))
scaledExprs[scaledExprs > -min(scaledExprs)] <- -min(scaledExprs)

palette_length = 100
my_color = my_color <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(palette_length)

my_breaks <- c(seq(min(scaledExprs), 0,
                   length.out=ceiling(palette_length/2) + 1),
               seq(max(scaledExprs)/palette_length,
                   max(scaledExprs),
                   length.out=floor(palette_length/2)))
```

![](KOGO_Seurat_pipeline_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
pdf(file = "KOGO_Seurat/heatmap.pdf")
pheatmap(scaledExprs,
              cluster_cols = T, cluster_rows = F, clustering_method = "ward.D2",
              treeheight_col = 0,
              breaks = my_breaks, color=my_color,
              labels_row = as.expression(lapply(rownames(scaledExprs), function(a) bquote(italic(.(a))))),
              angle_col = 315)
dev.off()
```

    ## pdf 
    ##   3

### **Cell type annotation**

``` r
seurat$celltype = as.character(seurat$seurat_clusters)

seurat$celltype[seurat$celltype %in% c('0')] <- 'NK.cell'
seurat$celltype[seurat$celltype %in% c('1','5','7','8')] <- 'T.cell'
seurat$celltype[seurat$celltype %in% c('2','6')] <- 'B.cell'
seurat$celltype[seurat$celltype %in% c('13')] <- 'Plasma.cell'
seurat$celltype[seurat$celltype %in% c('3','4','9','11')] <- 'Monocyte'
seurat$celltype[seurat$celltype %in% c('10')] <- 'Classical.Dendritic'
seurat$celltype[seurat$celltype %in% c('12')] <- 'Plasmacytoid.Dendritic'
```

``` r
p <- DimPlot(seurat, group.by = 'celltype', label=FALSE)
ggsave(filename = "KOGO_Seurat/annotation.png", plot = p, width = 9, height = 7)
```

### **Reference**

Pekayvaz, K., Leunig, A., Kaiser, R. et al. Protective immune
trajectories in early viral containment of non-pneumonic SARS-CoV-2
infection. Nat Commun 13, 1018 (2022).

Lun, A. T., McCarthy, D. J. & Marioni, J. C. A step-by-step workflow for
low-level analysis of single-cell RNA-seq data with Bioconductor.
F1000Res 5, 2122 (2016).

Yuhan Hao, Stephanie Hao, Erica Andersen-Nissen, William M. Mauck,
Shiwei Zheng, Andrew Butler, Maddie J. Lee, Aaron J. Wilk, Charlotte
Darby, Michael Zager, Paul Hoffman, Marlon Stoeckius, Efthymia Papalexi,
Eleni P. Mimitou, Jaison Jain, Avi Srivastava, Tim Stuart, Lamar M.
Fleming, Bertrand Yeung, Angela J. Rogers, Juliana M. McElrath,
Catherine A. Blish, Raphael Gottardo, Peter Smibert, Rahul Satija,
Integrated analysis of multimodal single-cell data,
Cell,184(13):3573-3587(2021).
