CellChat
================
2024-07-08

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
lung_seurat <- readRDS("/BiO/home/data/CellChat/lung_seurat.rds")
lung_seurat
```

    ## An object of class Seurat 
    ## 25916 features across 66452 samples within 2 assays 
    ## Active assay: RNA (23916 features, 2000 variable features)
    ##  1 other assay present: integrated
    ##  3 dimensional reductions calculated: pca, tsne, umap

``` r
UMAPPlot(lung_seurat, group.by = 'celltype', label = TRUE)
```
![image](https://github.com/CB-postech/Workshop-hands-on-materials/assets/98519284/64667fa6-b693-4d0e-a13e-54efafd5c80c)

``` r
UMAPPlot(lung_seurat, group.by = 'celltype', split.by = 'group', label = TRUE)
```
![image](https://github.com/CB-postech/Workshop-hands-on-materials/assets/98519284/6f4d2ca0-ea6c-4bab-b405-8cbe3c0c785d)

### **Preparing the data for CellChat**

``` r
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB) + theme_classic(base_line_size = 0)
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
# remove neutrophil and ambient cluster for analysis and select patient group (Severe covid: S/C)
covid <- rownames(meta)[meta$group == 'S/C' & meta$celltype %in% c('M1','M2','mDC','T','NK','B','Epithelial','pDC','Plasma','Mast')]
covid_input = expr[, covid]
covid_meta = meta[covid, ]
covid_meta$celltype = as.character(covid_meta$celltype)
unique(covid_meta$celltype)
```
    ## 'M1' 'T' Epithelial' 'NK' 'Plasma' 'M2' 'mDC' 'B' 'pDC' 'Mast'

``` r
cellchat_covid <- createCellChat(object = covid_input, meta = covid_meta, group.by = 'celltype')
```

    ## [1] "Create a CellChat object from a data matrix"
    ## Set cell identities for the new CellChat object 
    ## The cell groups used for CellChat analysis are  B Epithelial M1 M2 Mast mDC NK pDC Plasma T 

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
netVisual_circle(cellchat_covid@net$weight, vertex.weight = groupSize_covid,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
```
![image](https://github.com/CB-postech/Workshop-hands-on-materials/assets/98519284/e9526333-a65f-42e7-a30e-82d577433ed3)



### Process same process for healthy control(HC) group

```r
healthy <- rownames(meta)[meta$group == 'HC' & meta$celltype %in% c('M1','M2','mDC','T','NK','B','Epithelial','pDC','Plasma','Mast')]
healthy_input = expr[, healthy]
healthy_meta = meta[healthy, ]
healthy_meta$celltype = as.character(healthy_meta$celltype)
```

```r
cellchat_HC <- createCellChat(object = healthy_input, meta = healthy_meta, group.by = 'celltype')
```

    ## [1] "Create a CellChat object from a data matrix"
    ## Set cell identities for the new CellChat object 
    ## The cell groups used for CellChat analysis are  B Epithelial M1 M2 Mast mDC NK pDC Plasma T 

```r
cellchat_HC <- addMeta(cellchat_HC, meta = healthy_meta)
cellchat_HC <- setIdent(cellchat_HC, ident.use = 'celltype')
groupSize_HC <- as.numeric(table(cellchat_HC@idents))
cellchat_HC@DB <- CellChatDB
```

```r
# Preprocessing the expression data for cell-cell communication analysis
cellchat_HC <- subsetData(cellchat_HC)
cellchat_HC <- identifyOverExpressedGenes(cellchat_HC)
cellchat_HC <- identifyOverExpressedInteractions(cellchat_HC)
cellchat_HC <- projectData(cellchat_HC, PPI.human)

# Compute the communication probability and infer cellular communication network
cellchat_HC <- computeCommunProb(cellchat_HC)
cellchat_HC <- filterCommunication(cellchat_HC, min.cells = 0)

# Extract the inferred cellular communication network as a data frame
df.net_HC <- subsetCommunication(cellchat_HC)

# Infer the cell-cell communication at a signaling pathway level
cellchat_HC <- computeCommunProbPathway(cellchat_HC)

# Calculate the aggregated cell-cell communication network
cellchat_HC <- aggregateNet(cellchat_HC)

# Compute centrality
cellchat_HC <- netAnalysis_computeCentrality(cellchat_HC, slot.name = "netP")

```
### **Visualization**
```
netVisual_circle(cellchat_HC@net$weight, vertex.weight = groupSize_HC, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
```
![image](https://github.com/CB-postech/Workshop-hands-on-materials/assets/98519284/77257d17-3f71-46d4-b54b-aae3897542c0)





### **Reference**

Jin, S., Guerrero-Juarez, C. F., Zhang, L., Chang, I., Ramos, R., Kuan,
C. H., … & Nie, Q. (2021). Inference and analysis of cell-cell
communication using CellChat. Nature communications, 12(1), 1-20. 

Zhang, Z., Cui, F., Cao, C., Wang, Q., & Zou, Q. (2022). Single-cell RNA
analysis reveals the potential risk of organ-specific cell types
vulnerable to SARS-CoV-2 infections. Computers in biology and medicine,
140, 105092.
