KOGO_CellChat
================
2023-07-07

### **Library CellChat**

``` r
library(devtools)
library(Seurat)
library(scater)
library(harmony)
library(DropletUtils)
library(scran)
library(CellChat)
library(ggplot2)
library(patchwork)
library(igraph)
library(ComplexHeatmap)
```

### **Load the Data**

``` r
lung_seurat <- readRDS("/home/data/lung_seurat.rds")
lung_seurat
```

    ## An object of class Seurat 
    ## 25916 features across 66452 samples within 2 assays 
    ## Active assay: RNA (23916 features, 2000 variable features)
    ##  1 other assay present: integrated
    ##  3 dimensional reductions calculated: pca, tsne, umap

``` r
p <- UMAPPlot(lung_seurat, group.by = 'celltype', label = TRUE)
ggsave(filename = "KOGO_CellChat/umap.png", plot = p)
```

``` r
p <- UMAPPlot(lung_seurat, group.by = 'celltype', split.by = 'group', label = TRUE)
ggsave(filename = "KOGO_CellChat/umap-group.png", plot = p, width = 15, height = 6)
```

### **Preparing the data for CellChat**

``` r
CellChatDB <- CellChatDB.human
p <- showDatabaseCategory(CellChatDB) + theme_classic(base_line_size = 0)
ggsave(filename = "KOGO_CellChat/DB-category.png", plot = p)
```

``` r
expr <- lung_seurat@assays$RNA@data
meta <- lung_seurat@meta.data
```

### **Run CellChat**

1.  Create CellChat object

``` r
table(lung_seurat$disease)
```

    ## 
    ##     N     Y 
    ## 21939 44513

``` r
covid <- rownames(meta)[meta$group == 'S/C']
covid_input = expr[, covid]
covid_meta = meta[covid, ]
unique(covid_meta$celltype)
```

    ##  [1] M1         Neutrophil Ambient    T          Epithelial NK        
    ##  [7] Plasma     M2         mDC        B          pDC        Mast      
    ## Levels: B Plasma T NK pDC mDC Mast Neutrophil M1 M2 Epithelial Ambient

``` r
cellchat_covid <- createCellChat(object = covid_input, meta = covid_meta, group.by = 'celltype')
```

    ## [1] "Create a CellChat object from a data matrix"
    ## Set cell identities for the new CellChat object 
    ## The cell groups used for CellChat analysis are  B Plasma T NK pDC mDC Mast Neutrophil M1 M2 Epithelial Ambient

``` r
cellchat_covid <- addMeta(cellchat_covid, meta = covid_meta)
cellchat_covid <- setIdent(cellchat_covid, ident.use = 'celltype')

levels(cellchat_covid@meta$celltype)
```

    ##  [1] "B"          "Plasma"     "T"          "NK"         "pDC"       
    ##  [6] "mDC"        "Mast"       "Neutrophil" "M1"         "M2"        
    ## [11] "Epithelial" "Ambient"

``` r
groupSize_covid <- as.numeric(table(cellchat_covid@idents))
cellchat_covid@DB <- CellChatDB
```

2.  Run CellChat

``` r
# Preprocessing the expression data for cell-cell communication analysis
cellchat_covid <- subsetData(cellchat_covid)
cellchat_covid <- identifyOverExpressedGenes(cellchat_covid)
cellchat_covid <- identifyOverExpressedInteractions(cellchat_covid)
cellchat_covid <- projectData(cellchat_covid, PPI.human)

# Compute the communication probability and infer cellular communication network
cellchat_covid <- computeCommunProb(cellchat_covid)
```

    ## triMean is used for calculating the average gene expression per cell group. 
    ## [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2023-07-04 14:54:58]"
    ## [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2023-07-04 14:58:18]"

``` r
cellchat_covid <- filterCommunication(cellchat_covid, min.cells = 0)

# Extract the inferred cellular communication network as a data frame
df.net_covid <- subsetCommunication(cellchat_covid)

# Infer the cell-cell communication at a signaling pathway level
cellchat_covid <- computeCommunProbPathway(cellchat_covid)

# Calculate the aggregated cell-cell communication network
cellchat_covid <- aggregateNet(cellchat_covid)

# Compute centrality
cellchat_covid <- netAnalysis_computeCentrality(cellchat_covid, slot.name = "netP")
```

### **Visualization**

``` r
pdf(file = "KOGO_CellChat/netVisual_circle.pdf")
netVisual_circle(cellchat_covid@net$weight, vertex.weight = groupSize_covid,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
```

    ## png 
    ##   2

### **Reference**

Jin, S., Guerrero-Juarez, C. F., Zhang, L., Chang, I., Ramos, R., Kuan,
C. H., … & Nie, Q. (2021). Inference and analysis of cell-cell
communication using CellChat. Nature communications, 12(1), 1-20. Zhang,
Z., Cui, F., Cao, C., Wang, Q., & Zou, Q. (2022). Single-cell RNA
analysis reveals the potential risk of organ-specific cell types
vulnerable to SARS-CoV-2 infections. Computers in biology and medicine,
140, 105092.
