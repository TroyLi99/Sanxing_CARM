DimPlot(sce, reduction = "umap", label = TRUE, group.by = "seurat_clusters")

#differential marker genes identification
markers <- FindAllMarkers(sce,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 1)
top10.markers <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

#initial cell annotation
#Hepatocytic cancer cell (with CNV analysis)
DotPlot(sce, group.by = "seurat_clusters", features = c("AFP", "ALB", "EPCAM"), cols = c("gray", "red"))
#Monocyte/Macrophage
DotPlot(sce, group.by = "seurat_clusters", features = c("PTPRC", "CD14", "FCGR3A", "CD68", "CSF1R", "S100A9"), cols = c("white", "red"))
#Progenitor cell (with celltypist annotation transfer)
DotPlot(sce, group.by = "seurat_clusters", features = c("SPI1", "RUNX1", "CEBPA", "IRF8", "MYB"), cols = c("gray", "red"))
#NK cell
DotPlot(sce, group.by = "seurat_clusters", features = c("NCAM1", "FCGR3A", "CD7", "KLRD1", "KLRF1", "NKG7", "GNLY", "GZMB"), cols = c("gray", "red"))
#T cell
DotPlot(sce, group.by = "seurat_clusters", features = c("CD8A", "CD8B", "CD3D", "GNLY", "NKG7", "CD4", "TRGC1", "TRGC2"), cols = c("gray", "red"))
#B cell
DotPlot(sce, group.by = "seurat_clusters", features = c("CD19", "CD79A", "MS4A1"), cols = c("gray", "red"))

#For cell subtype identification, each cluster is first screened out and annotated with pathway and gene expression
#GSVA
library(GSVA)
library(GSEABase)
library(tidyverse)
library(org.Hs.eg.db)
library(pheatmap)

expr <- AverageExpression(sce, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,] 
expr <- as.matrix(expr)
head(expr)

url <- "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset" #each cell's specific download_geneset can be found in supplementary materials
file_name <- "xxx.gmt"
download.file(url, file_name)
GeneSet <- getGmt("./xxx.gmt")

result <- gsvaParam(expr, GeneSet, maxDiff = T)
gsva_es2 <- gsva(result)
result_c <- rbind(result_c, gsva_es2)
rownames(result_c) <- c("xxx", "xxx", "xxx", "xxx")
colnames(result_c) <- c("xxx", "xxx", "xxx", "xxx")
#plot
pheatmap(result_c, show_colnames = T, 
         scale = "row",angle_col = "315",
         cluster_row = F,cluster_col = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

#GO and KEGG pathway enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tibble)

run_enrichment <- function(cluster_id) {
  cluster_genes <- markers %>%
    filter(cluster == cluster_id & p_val_adj < 0.05) %>%
    pull(gene)
  
  gene_ids <- bitr(cluster_genes, 
                   fromType = "SYMBOL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)

  go_res <- enrichGO(gene = gene_ids$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "BP", 
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)
  
  kegg_res <- enrichKEGG(gene = gene_ids$ENTREZID,
                         organism = "hsa",
                         pvalueCutoff = 0.05)
  
  return(list(go = go_res, kegg = kegg_res))
}


enrich_results <- lapply(levels(sce), run_enrichment)
names(sce) <- levels(sce)

# GO
go_terms <- lapply(names(enrich_results), function(cl){
  res <- enrich_results[[cl]]$go
  if(nrow(res) > 0){
    res@result %>%
      mutate(Cluster = cl) %>%
      slice_min(pvalue, n = 10) 
  }
}) %>% bind_rows()

ggplot(go_terms, aes(x = -log10(pvalue), y = Description, fill = Cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Cluster, scales = "free_y") +
  labs(title = "Top GO Biological Processes in Each Subtype",
       x = "-log10(p-value)",
       y = "") +
  theme_minimal()

# KEGG
combined_kegg <- merge_result(lapply(enrich_results, function(x) x$kegg))
combined_kegg@compareClusterResult$Description <- gsub(" - Hs \\(human\\)", "", 
                                                       combined_kegg@compareClusterResult$Description)
combined_kegg@compareClusterResult$Description

heatmap_data <- combined_kegg@compareClusterResult %>%
  group_by(Cluster) %>%
  slice_min(pvalue, n = 10) %>%
  acast(Description ~ Cluster, value.var = "pvalue", fill = 1)

pheatmap(-log10(heatmap_data),
         color = colorRampPalette(c("white", "red"))(100),
         main = "KEGG Pathway Enrichment Heatmap",
         angle_col = 45)









