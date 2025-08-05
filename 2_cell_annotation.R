#initial cell annotation
DimPlot(sce, reduction = "umap", label = TRUE, group.by = "seurat_clusters")

#differential marker genes identification
markers <- FindAllMarkers(sce,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 1)
top10.markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

#Hepatocytic cancer cell (with CNV analysis)
DotPlot(sce, group.by = "seurat_clusters", features = c("AFP", "ALB", "EPCAM"), cols = c("gray", "red"))
#Monocyte/Macrophage
DotPlot(sce, group.by = "seurat_clusters", features = c("PTPRC", "CD14", "FCGR3A", "CD68", "CSF1R", "S100A9"), cols = c("white", "red"))
#Progenitor cell (with celltypist annotation transfer)
DotPlot(sce, group.by = "seurat_clusters", features = c("SPI1", "RUNX1", "CEBPA", "IRF8", "MYB"), cols = c("gray", "red"))
#NK cell
DotPlot(sce, group.by = "seurat_clusters", features = c("NCAM1", "FCGR3A", "CD7", "KLRD1", "KLRF1", "NKG7", "GNLY", "GZMB"), cols = c("gray", "red"))
#T cell

#B cell


#For cell subtype identification, the cluster is first screened out and annotated with pathway and gene expression


