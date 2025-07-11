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
seurat <- readRDS("/BiO/data/HLCA_pulmonary_fibrosis_immune_sampled.rds")
seurat
```

    ## An object of class Seurat 
    ## 19354 features across 7408 samples within 1 assay 
    ## Active assay: RNA (19354 features, 0 variable features)
    ##  3 layers present: counts, data, scale.data
    ##  4 dimensional reductions calculated: pca, umap, harmony, umap.harmony

### **UMAP visualization**

``` r
DimPlot(seurat, group.by = 'celltype', reduction = "umap.harmony")
```
<img width="400" height="300" alt="image" src="https://github.com/user-attachments/assets/bbbf95c3-7d8e-449b-90b7-7825a3994a35" />

``` r
DimPlot(seurat, group.by = 'celltype', reduction = "umap.harmony", split.by = "disease")
```
<img width="600" height="320" alt="image" src="https://github.com/user-attachments/assets/86a5a685-67ab-46ca-868d-c8820abf3c48" />

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
<img width="500" height="400" alt="image" src="https://github.com/user-attachments/assets/792c759b-f82d-41d0-9edc-b7000450b409" />

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
<img width="300" height="300" alt="image" src="https://github.com/user-attachments/assets/e5f9a27e-832a-4010-aea8-f77723f74cf4" />

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
<img width="400" height="200" alt="image" src="https://github.com/user-attachments/assets/87c8bee5-e98f-43a1-b930-c22317d7a2ff" />

``` r
netVisual_bubble(cellchat_PF, sources.use = c(3:4), targets.use = c(1:2), signaling = c("TGFb"), remove.isolate = TRUE)
```
<img width="300" height="300" alt="image" src="https://github.com/user-attachments/assets/7d567df2-e49d-4de0-a5f1-5cef6d6dc719" />

``` r
netVisual_chord_gene(cellchat_PF, sources.use = c(3:4), targets.use = c(1:2), signaling = c("TGFb"))
```
<img width="300" height="320" alt="image" src="https://github.com/user-attachments/assets/9c8a0258-e66e-4249-9512-c3e47a9bb3e3" />

---
### Process same process for normal group (to be updated)

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

---
### **Reference**

Jin, S., Guerrero-Juarez, C. F., Zhang, L., Chang, I., Ramos, R., Kuan, C. H., … & Nie, Q. (2021). Inference and analysis of cell-cell communication using CellChat. Nature communications, 12(1), 1-20. 

Jin, S., Plikus, M.V. & Nie, Q. (2025). CellChat for systematic analysis of cell–cell communication from single-cell transcriptomics. Nat Protoc, 20, 180–219.

Sikkema, L., Ramírez-Suástegui, C., Strobl, D.C., … & Theis, F.J. (2023). An integrated cell atlas of the lung in health and disease. Nat Med, 29, 1563–1577.
