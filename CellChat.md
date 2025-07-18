Cell-Cell Interaction Analysis with CellChat
================
2025-07-22

---
### **Load library**

``` r
library(Seurat)
library(CellChat)
```

### **Load the Data**

``` r
seurat <- readRDS("/BiO/data/HLCA_pulmonary_fibrosis_immune.rds")
seurat
```

    ## An object of class Seurat 
    ## 19354 features across 7782 samples within 1 assay 
    ## Active assay: RNA (19354 features, 0 variable features)
    ##  3 layers present: counts, data, scale.data
    ##  4 dimensional reductions calculated: pca, umap, harmony, umap.harmony

### **UMAP visualization**

``` r
DimPlot(seurat, group.by = 'celltype', reduction = "umap.harmony")
```

<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/122061dd-b366-4e14-938e-8a305e3dcaa6" />

``` r
DimPlot(seurat, group.by = 'celltype', reduction = "umap.harmony", split.by = "disease")
```

<img width="600" height="320" alt="image" src="https://github.com/user-attachments/assets/217f66a5-7c18-49d1-8f44-0197c1be3b28" />

---
### **Input data processing**
**1-1. Create CellChat object from Seurat object**
``` r
seurat_PF <- subset(seurat, subset = disease == "pulmonary fibrosis")
cellchat_PF <- createCellChat(object = seurat_PF, group.by = "celltype", assay = "RNA")
```

**1-2. Create CellChat object from expression matrix and metadata**
``` r
expr <- seurat[["RNA"]]$data
meta <- seurat@meta.data

cells_PF <- rownames(meta)[meta$disease == "pulmonary fibrosis"]

expr_PF <- expr[, cells_PF]
meta_PF <- meta[cells_PF,]
meta_PF$celltype <- as.character(meta_PF$celltype)

cellchat_PF <- createCellChat(object = expr_PF, meta = meta_PF, group.by = "celltype")
```

**2. Set the ligand-receptor interaction database**
``` r
CellChatDB <- CellChatDB.human # # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
```

<img width="1770" height="312" alt="image" src="https://github.com/user-attachments/assets/a52532d5-fec7-4073-b64d-7435b1975968" />


**3. Subset the expression data using CellChatDB genes**
``` r
cellchat_PF@DB <- CellChatDB
cellchat_PF <- subsetData(cellchat_PF)
```

**4. Identify over-expressed ligands/receptors and L-R interactions in each cell group**
``` r
cellchat_PF <- identifyOverExpressedGenes(cellchat_PF)
cellchat_PF <- identifyOverExpressedInteractions(cellchat_PF)
```

    ## The number of highly variable ligand-receptor pairs used for signaling inference is 872

**5. (Optional) Smooth the gene expression because of shallow sequencing depth**
``` r
cellchat_PF <- smoothData(cellchat_PF, adj = PPI.human)
```

---
### **Inference of cell-cell communication networks**

**1. Compute the communication probability**
``` r
cellchat_PF <- computeCommunProb(cellchat_PF, raw.use = FALSE) # Set raw.use = FALSE to use the smoothed data
```
    
    ## triMean is used for calculating the average gene expression per cell group. 
    ## [1] ">>> Run CellChat on sc/snRNA-seq data <<< [2025-07-09 16:42:49.714313]"
    ## [1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2025-07-09 16:44:45.242436]"

**2. Filter the cell-cell interaction, based on the number of cells in each group**
``` r
cellchat_PF <- filterCommunication(cellchat_PF, min.cells = 0)
```

**3. Extract the inferred cellular communication network as a data frame**
``` r
df_net_PF <- subsetCommunication(cellchat_PF)
```
<img width="1770" height="299" alt="image" src="https://github.com/user-attachments/assets/55c45fb2-c047-44b2-bf08-67e27712b427" />


**4. Infer the cell-cell communication at a signaling pathway level**
``` r
cellchat_PF <- computeCommunProbPathway(cellchat_PF)
```

**5. Calculate the aggregated cell-cell communication network**
``` r
cellchat_PF <- aggregateNet(cellchat_PF)
```

**6. Save the CellChat object**
``` r
saveRDS(cellchat_PF, file = "/BiO/home/edu03/cellchat_lung_PF.rds")
```

---
### **Visualization**
**Circle plot**
``` r
groupSize_PF <- as.numeric(table(cellchat_PF@idents))
netVisual_circle(cellchat_PF@net$weight, vertex.weight = groupSize_PF, weight.scale = T, label.edge= F)
```
<img width="300" height="300" alt="image" src="https://github.com/user-attachments/assets/4deb6b30-54f5-4f46-a51b-a47da4354cca" />

---
### **Identify signaling roles and major contributing signaling**
**1. Compute the network centrality scores**
``` r
cellchat_PF <- netAnalysis_computeCentrality(cellchat_PF, slot.name = "netP")
```

**2. Visualization**
``` r
pathways.show <- "TGFb"
netAnalysis_signalingRole_network(cellchat_PF, signaling = pathways.show, width = 6, height = 2, font.size = 10)
```
<img width="400" height="200" alt="image" src="https://github.com/user-attachments/assets/c08a1806-da05-4332-9d89-8664c369c749" />

``` r
netVisual_bubble(cellchat_PF, sources.use = c(3:4), targets.use = c(1:2), signaling = c("TGFb"), remove.isolate = TRUE)
```
<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/1d6c1356-573c-4baf-88cd-deb657c47892" />

``` r
netVisual_chord_gene(cellchat_PF, sources.use = c(3:4), targets.use = c(1:2), signaling = c("TGFb"))
```
<img width="300" height="320" alt="image" src="https://github.com/user-attachments/assets/b58c211f-ba2d-48ac-9ee9-d789e717282b" />

---
### Process same process for normal group
**1. Input data processing & Running CellChat**
```r
seurat_NM <- subset(seurat, subset = disease == "normal")
cellchat_NM <- createCellChat(object = seurat_NM, group.by = "celltype", assay = "RNA")

cellchat_NM@DB <- CellChatDB

cellchat_NM <- subsetData(cellchat_NM)

cellchat_NM <- identifyOverExpressedGenes(cellchat_NM)
cellchat_NM <- identifyOverExpressedInteractions(cellchat_NM)

cellchat_NM <- smoothData(cellchat_NM, adj = PPI.human)

cellchat_NM <- computeCommunProb(cellchat_NM, raw.use = FALSE)
cellchat_NM <- filterCommunication(cellchat_NM, min.cells = 0)
cellchat_NM <- computeCommunProbPathway(cellchat_NM)
cellchat_NM <- aggregateNet(cellchat_NM)

```
**2. Visualization**
``` r
groupSize_NM <- as.numeric(table(cellchat_NM@idents))
netVisual_circle(cellchat_NM@net$weight, vertex.weight = groupSize_NM,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
```
<img width="300" height="300" alt="image" src="https://github.com/user-attachments/assets/b6d96df9-d78f-43f8-be31-376cf26baf6c" />

``` r
# Compute centrality
cellchat_NM <- netAnalysis_computeCentrality(cellchat_NM, slot.name = "netP")

pathways.show <- "TGFb"
netAnalysis_signalingRole_network(cellchat_NM, signaling = pathways.show, width = 6, height = 2, font.size = 10)
```
<img width="400" height="200" alt="image" src="https://github.com/user-attachments/assets/099a98fc-f133-4b74-842d-78cb731a1891" />

---
### **Reference**

Jin, S., Guerrero-Juarez, C. F., Zhang, L., Chang, I., Ramos, R., Kuan, C. H., … & Nie, Q. (2021). Inference and analysis of cell-cell communication using CellChat. Nature communications, 12(1), 1-20. 

Jin, S., Plikus, M.V. & Nie, Q. (2025). CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics. Nat Protoc, 20, 180–219.

Sikkema, L., Ramírez-Suástegui, C., Strobl, D.C., … & Theis, F.J. (2023). An integrated cell atlas of the lung in health and disease. Nat Med, 29, 1563–1577.
