library(Seurat)
library(tidyverse)
library(dbplyr)
library(dplyr)
set.seed(123)

tmp = list.dirs('~/SCP/sanxing_project/1/')[-1]
tmp
ct = Read10X('~/SCP/sanxing_project/1/aGPC3-CAR-DR18-T-SIRPa//')
sce.all=CreateSeuratObject(counts = ct,
                           min.cells = 3,
                           min.features = 200)
#rename seurat_obj into orig.ident
new.cluster.ids <-  "seurat_obj"
names(new.cluster.ids) <- levels(sce.all$orig.ident)
sce.all <- RenameIdents(sce.all, new.cluster.ids)
gplots::balloonplot(table(sce.all$orig.ident,Idents(sce.all)))
sce.all$orig.ident = Idents(sce.all)

#identify human and mouse cells in samples
tail(rownames(sce.all), 3)
head(rownames(sce.all), 3)
sce.all[["prop_mouse"]] <- PercentageFeatureSet(sce.all, pattern = "^GRCm39-")
sce.all[["prop_human"]] <- PercentageFeatureSet(sce.all, pattern = "^GRCh38-")
head(sce.all@meta.data)
VlnPlot(sce.all,features =c("prop_human","prop_mouse"))
cell_species <- sce.all@meta.data %>% 
  mutate(species = case_when(
    prop_mouse > 85 ~ "mouse",
    prop_human > 85 ~ "human",
    TRUE ~ "mixed"
  ))
sce.all <- AddMetaData(sce.all, metadata = cell_species$species, col.name = "species")
head(sce.all@meta.data)
table(sce.all@meta.data$species)

#observe the distribution of prop_human, pro_mouse in each sample
library(ggplot2)
meta <- sce.all@meta.data
prop <- data.frame(
  proportion = c(sce.all$prop_mouse, sce.all$prop_human),
  species = rep(c("mouse", "human"), each = nrow(meta))
)

ggplot(prop, aes(x = proportion, fill = species)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  scale_fill_manual(values = c("mouse" = "blue", "human" = "red")) +
  labs(title = "Human_Mouse Proportion Distribution in SIRPa",
       x = "Proportion",
       y = "Count") +
  theme_minimal()

#only retain human cells and clean gene names
sce.all <- subset(sce.all, species != "mixed")
sce.all <- subset(sce.all, species == "human")
head(rownames(sce.all))
new_rownames <- gsub("^GRCh38-", "", rownames(sce.all))
rownames(sce.all) <- new_rownames
rownames(sce.all)
genes_to_remove <- grep("GRCm39-", rownames(sce.all), value = TRUE)
print(genes_to_remove)
sce.all <- sce.all[!rownames(sce.all) %in% genes_to_remove, ]
tail(rownames(sce.all))

#before filter
Idents(sce.all) <- sce.all@meta.data$orig.ident
sce.all[["percent.mt"]] <- PercentageFeatureSet(sce.all, pattern = "^MT-")
VlnPlot(sce.all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

nFeature_RNA_data <- data.frame(cell = rownames(sce.all@meta.data), 
                                nFeature_RNA = sce.all@meta.data$nFeature_RNA)
ggplot(nFeature_RNA_data, aes(x = nFeature_RNA)) +
  geom_bar(stat = "bin", fill = "black") +
  labs(x = "Number of Features", y = "Count", title = "Distribution of nFeature_RNA across cells") +
  theme_minimal()

plot1 <- FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce.all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
#after filter
sce.all_f <- subset(sce.all, subset = nFeature_RNA > xx & nFeature_RNA < xx & percent.mt < xx) #filtering threshold can be adjusted and found in supplementary materials
VlnPlot(sce.all_f, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(sce.all_f, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce.all_f, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#normalization, scale, PCA, UMAP
sce.all_f <- NormalizeData(sce.all_f, normalization.method = "LogNormalize", scale.factor = 10000)
sce.all_f <- FindVariableFeatures(sce.all_f, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(sce.all_f), 10)
plot1 <- VariableFeaturePlot(sce.all_f)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(sce.all_f)
sce.all_f <- ScaleData(sce.all_f, features = all.genes)

sce.all_f <- RunPCA(sce.all_f, features = VariableFeatures(object = sce.all_f), dims= 1:30)
DimPlot(sce.all_f, reduction = "pca")
ElbowPlot(sce.all_f)

sce.all_f <- FindNeighbors(sce.all_f, dims = 1:20)
sce.all_f <- FindClusters(sce.all_f, resolution = 0.5)
sce.all_f <- RunUMAP(sce.all_f, dims = 1:30)
DimPlot(sce.all_f, reduction = "umap")
DimPlot(sce.all_f, reduction = "umap", pt.size = 0.5, group.by = "orig.ident")

#doublet detection
library(DoubletFinder)
sweep.res.list <- paramSweep(sce, PCs = 1:20, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
sweep.stats[order(sweep.stats$BCreal),]
bcmvn <- find.pK(sweep.stats) 
pK_bcmvn <- as.numeric(bcmvn$pK[which.max(bcmvn$BCmetric)]) 
homotypic.prop <- modelHomotypic(sce.all_f$seurat_clusters) 
nExp_poi <- round(x *nrow(sce.all_f@meta.data))   # x value can be adjusted referring to the instruction of DoubletFinder package
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop)) 

sce.all_f <- doubletFinder(sce.all_f, PCs = 1:20, pN = 0.25, pK = pK_bcmvn,
                           nExp = nExp_poi.adj, reuse.pANN = F, sct = F)
DimPlot(sce, reduction = "umap", group.by = "DF.classifications_0.25_31_56", raster = FALSE)
colnames(sce.all_f@meta.data)
sce <- subset(sce.all_f, subset = DF.classifications_0.25_31_56 == "Singlet")
saveRDS(sce, file = "X.rds")