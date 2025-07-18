Monocle3
================
2025-07-22

## Import required libraries

```r
library(Seurat)
library(tidyverse)
library(monocle3)
```

## Configure output directory
```r
OUT = "./monocle-outputs"
if (!dir.exists(OUT)) { dir.create(OUT, recursive = TRUE) }
```


## Trajectory analysis using monocle3

Here, we describe a brief trajectory analysis of T cell subset using monocle3.
The dataset has various celltypes including T cell.

```r
seurat <- readRDS("/BiO/data/HLCA_pulmonary_fibrosis_immune.rds")
seurat
```

```r
DimPlot(seurat, reduction = "umap.harmony", group.by = "celltype",
    label = TRUE, label.size = 10, repel = TRUE) +
    NoLegend() + ggtitle(NULL)
ggsave(paste0(OUT, "/umap-total.png"), width = 7, height = 7, dpi = 300)
```

<img src="/images/monocle/umap-total.png" width="400" height="400">

Since we want to draw a trajectory graph of T cells, we will subset only the T cells from the whole dataset and re-normalize.

```r
seurat_t <- seurat[, seurat$celltype == "T cell"]
# normalization
seurat_t <- NormalizeData(seurat_t)
```

```r
hvg <- SelectIntegrationFeatures(
    object.list = SplitObject(seurat_t, split.by = "study"),
    assay = c("RNA", "RNA", "RNA", "RNA"),
    nfeatures = 2000, fvf.nfeatures = 2000
)
VariableFeatures(seurat_t) <- hvg
```

### Generate CDS object

Monocle3 package uses differently structured object named cell_data_set (cds)
We recombine the normalized expressions, metadata for cells, and metadata for genes to create the cds object.

```r
cds <- new_cell_data_set(
    expression_data = seurat_t@assays$RNA$data,
    cell_metadata = seurat_t@meta.data,
    gene_metadata = data.frame(
        gene_short_name = row.names(seurat_t),
        row.names = row.names(seurat_t)
        )
    )
cds
```

### Dimension reduction for CDS object

Monocle3 allows dimension reduction using hvgs. As we import normalized count in cds object, we preprocess the object without additional normalization.

``` r
cds <- preprocess_cds(cds, "PCA", num_dim = 30, norm_method = "none", use_genes = hvg)
```

We than see the explained variance of each component to select the optimal number of components to use

```r
plot_pc_variance_explained(cds) +
    theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)
        )
ggsave(paste0(OUT, "/elbow-PCA.png"), width = 7, height = 7, dpi = 300)
```

<img src="/images/monocle/elbow-PCA.png" width="400" height="400">

```r
cds <- preprocess_cds(cds, "PCA", num_dim = 10, norm_method = "none", use_genes = hvg)
```

### Correcting Batch effects

Since the dataset has various sample and batch effects, we perform Mutual Nearest Neighbor (MNN) batch effect correction implemented batchelor, which is included in monocle3 package.
The sample ID information is in 'study' metadata.

```r
# before batch correction
cds <- reduce_dimension(cds, preprocess_method = "PCA")

plot_cells(cds, color_cells_by = "study",
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           group_label_size = 10,
           cell_size = 1, cell_stroke = 2) +
    theme(
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)
        )
ggsave(paste0(OUT, "/monocle_umap-study-bfCorrection.png"), width = 9, height = 7, dpi = 300)
```

<img src="/images/monocle/monocle_umap-study-bfCorrection.png" width="450" height="350">

```r
cds <- align_cds(cds, alignment_group = "study")
cds <- reduce_dimension(cds, preprocess_method = 'Aligned')
```

```r
# after batch correction
plot_cells(cds, color_cells_by = "study",
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           group_label_size = 10,
           cell_size = 1, cell_stroke = 2) +
    theme(
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)
        )
ggsave(paste0(OUT, "/monocle_umap-study-afCorrection.png"), width = 9, height = 7, dpi = 300)
```

<img src="/images/monocle/monocle_umap-study-afCorrection.png" width="450" height="350">

### Cluster cells and learn the trajectory graph

After clustering, we will fit a principal graph within each partition using the learn_graph() function.

```r
# clustering
cds <- cluster_cells(cds, resolution = 0.0015)

plot_cells(cds, color_cells_by = "cluster",
           show_trajectory_graph = FALSE,
           group_label_size = 10,
           cell_size = 1, cell_stroke = 2)
ggsave(paste0(OUT, "/monocle_umap-clusters.png"), width = 7, height = 7, dpi = 300)
```

<img src="/images/monocle/monocle_umap-clusters.png" width="350" height="350">

```r
# learn graph
cds = learn_graph(cds)
```

```r
plot_cells(cds, color_cells_by = "cluster",
           show_trajectory_graph = TRUE,
           label_principal_points = TRUE,
           label_cell_groups = FALSE,
           graph_label_size = 3,
           trajectory_graph_color = "firebrick1", trajectory_graph_segment_size = 2,
           cell_size = 1, cell_stroke = 1
          ) +
    theme(
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)
    )
ggsave(paste0(OUT, "/monocle_trajectory-clusters.png"), width = 8, height = 7, dpi = 300)
```

<img src="/images/monocle/monocle_trajectory-clusters.png" width="400" height="350">

### Order cells in pseudotime

CCR7 and LEF1 are known as naive T cell marker, so we set the cluster where expression of these genes is high as the root state.

```r
plot_cells(cds, genes = c("CCR7", "LEF1"),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           graph_label_size = 3,
           cell_size = 1, cell_stroke = 2
          ) +
    labs(color = "Expression") +
    theme(
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)
    ) +
    scale_color_viridis_c(option = "magma")
ggsave(paste0(OUT, "/monocle_root-expression.png"), width = 10.5, height = 5, dpi = 300)
```

<img src="/images/monocle/monocle_root-expression.png" width="525" height="250">

```r
# order cells while setting root principal node
cds <- order_cells(cds, root_pr_nodes = 'Y_72')
```

Plotting the cells and coloring them by pseudotime shows how they were ordered.

```r
plot_cells(cds, color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE,
           label_principal_points = TRUE, label_leaves = FALSE, label_branch_points = FALSE,
           label_cell_groups = FALSE,
           graph_label_size = 3,
           trajectory_graph_color = "firebrick1", trajectory_graph_segment_size = 2,
           cell_size = 1, cell_stroke = 2
          ) +
    labs(color = "Pseudotime") +
    theme(
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)
    ) +
    scale_color_viridis_c(option = "viridis")
ggsave(paste0(OUT, "/monocle_pseudotime.png"), width = 8, height = 7, dpi = 300)
```

<img src="/images/monocle/monocle_pseudotime.png" width="400" height="350">

Check the expression of genes related to t cell function.

```r
# CD4+ / CD8+ T cells
plot_cells(cds, genes = c("CD4", "CD8A"),
           show_trajectory_graph = FALSE,
           label_principal_points = TRUE,
           label_cell_groups = FALSE,
           graph_label_size = 3,
           cell_size = 1, cell_stroke = 2
          ) +
    labs(color = "Expression") +
    theme(
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)
    ) +
    scale_color_viridis_c(option = "magma")
ggsave(paste0(OUT, "/monocle_expression1.png"), width = 10.5, height = 5, dpi = 300)
```

<img src="/images/monocle/monocle_expression1.png" width="525" height="250">

```r
# Cytotoxic CD8 T cells
plot_cells(cds, genes = c("CCL5", "GZMK", "GNLY", "NKG7"),
           show_trajectory_graph = FALSE,
           label_principal_points = TRUE,
           label_cell_groups = FALSE,
           graph_label_size = 3,
           cell_size = 1, cell_stroke = 1
          ) +
    labs(color = "Expression") +
    theme(
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15)
    ) +
    scale_color_viridis_c(option = "magma")
ggsave(paste0(OUT, "/monocle_expression2.png"), width = 10.5, height = 10, dpi = 300)
```

<img src="/images/monocle/monocle_expression2.png" width="525" height="500">

## Reference

Butler, A., Hoffman, P., Smibert, P., Papalexi, E. & Satija, R. Integrating single-cell transcriptomic data across different conditions, technologies, and species. Nat. Biotechnol. 36, 411–420 (2018).

Cao, J. et al. The single-cell transcriptional landscape of mammalian organogenesis. Nature 566, 496–502 (2019).

Haghverdi L, Lun ATL, Morgan MD, Marioni JC (2018). 'Batch effects in single-cell RNA-sequencing data are corrected by matching mutual nearest neighbors.' Nat. Biotechnol., 36(5), 421-427. doi: 10.1038/nbt.4091

L. Ma, M.O. Hernandez, Y. Zhao, M. Mehta, B. Tran, M. Kelly, Z. Rae, J.M. Hernandez, J.L. Davis, S.P. Martin, D.E. Kleiner, S.M. Hewitt, K. Ylaya, B.J. Wood, T.F. Greten, X.W. Wang. Tumor cell biodiversity drives microenvironmental reprogramming in liver cancer. Canc. Cell, 36 (4): 418-430 (2019)

Sikkema, L., Ramírez-Suástegui, C., Strobl, D.C. et al. An integrated cell atlas of the lung in health and disease. Nat Med 29, 1563–1577 (2023)
