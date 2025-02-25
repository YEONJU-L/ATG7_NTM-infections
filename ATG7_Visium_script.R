
#### ATG7_Visum analyses script.R ####

# load library ------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(harmony)
library(pheatmap)
library(ape)

# Load data ---------------------------------------------------------------
wt <- Load10X_Spatial(
  "Spatial_matrix_wt",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper=FALSE,
  image=NULL)

ko <- Load10X_Spatial(
  "Spatial_matrix_ko/",
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE,
  to.upper=FALSE,
  image=NULL)

wt$orig.ident <- "ATG7_WT"
ko$orig.ident <- "ATG7_KO"

wt$barcodes <- rownames(wt@meta.data)
ko$barcodes <- rownames(ko@meta.data)

# QC -----------------------------------------------------------------

summary(wt$nCount_Spatial)
summary(ko$nCount_Spatial)

summary(wt$nFeature_Spatial)
summary(ko$nFeature_Spatial)

atg <- merge(wt, ko, add.cell.ids = c("ATG7_WT", "ATG7_KO"))

VlnPlot(atg, c("nCount_Spatial","nFeature_Spatial"), cols = c("#66c2a5", "#fc8d62"), group.by = "orig.ident")

SpatialFeaturePlot(atg, features = "nCount_Spatial") + theme(legend.position = "right")
SpatialFeaturePlot(atg, features = "nFeature_Spatial") + theme(legend.position = "right")

wt$lowQC <- ifelse(wt$nFeature_Spatial >300, "PASS", "FAIL")
ko$lowQC <- ifelse(ko$nFeature_Spatial >300, "PASS", "FAIL")

plot1 <- FeatureScatter(wt, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", group.by = "lowQC", cols = c("PASS" = "black", "FAIL" = "red")) + theme(legend.position = "none")

plot2 <- FeatureScatter(ko, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", group.by = "lowQC", cols = c("PASS" = "black", "FAIL" = "red"))

plot2 <- plot2 + 
  scale_x_continuous(breaks = seq(0, 60000, by = 20000), limits = c(0, 60000)) +
  scale_y_continuous(breaks = seq(0, 10000, by = 2500), limits = c(0, 10000))

wt <- subset(wt, subset = lowQC == "PASS") 
ko <- subset(ko, subset = lowQC == "PASS") 

atg <- merge(wt, ko, add.cell.ids = c("ATG7_WT", "ATG7_KO"))


# Processing --------------------------------------------------------------

atg <- JoinLayers(atg)

atg <- NormalizeData(atg, normalization.method = "LogNormalize")
atg <- FindVariableFeatures(atg, assay = "Spatial", selection.method = "vst", nfeatures = 2000)
atg <- ScaleData(atg, features = rownames(atg))

atg <- RunPCA(atg, features = VariableFeatures(object = atg))
ElbowPlot(atg, ndims =50)

atg <- atg %>% RunUMAP(dims = 1:10) %>%
  FindNeighbors(dims = 1:10) %>% FindClusters(verbose = FALSE, resolution = 1.0)

DimPlot(atg, split.by = "orig.ident", label = T)

atg <- RunHarmony(atg, group.by.vars = "orig.ident")

atg <- atg %>% RunUMAP(reduction = "harmony", dims = 1:10) %>%
  FindNeighbors(reduction = "harmony", dims = 1:10) %>% 
  FindClusters(verbose = FALSE, resolution = 1.0)


# Nodal clusters group (cluster tree) ------------------------------------------------------------

atg <- BuildClusterTree(
  atg,
  assay = NULL,
  features = NULL,
  dims = 1:10,
  reduction = "pca",
  graph = NULL,
  slot = "data",
  reorder = FALSE,
  reorder.numeric = FALSE,
  verbose = TRUE
)

a <- Tool(object = atg, slot = "BuildClusterTree")

ape::plot.phylo(x = a, direction = "rightwards", label.offset=TRUE)

# group annotation
Idents(atg) <- atg$seurat_clusters

new.cluster.ids <- c("0" = "2", "1" ="2", "2" ="5", "3" ="2", "4" ="6", "5" ="5",
                     "6" ="3", "7" ="2", "8"="5", "9"="5", "10"="6",
                     "11"="1", "12"="4", "13"="1", "14"="3", "15"="6", "16"="7", "17"="8")

atg <- RenameIdents(atg, new.cluster.ids)


atg$clusters <- atg@active.ident

atg@meta.data$clusters <- factor(atg@meta.data$clusters, 
                                 levels = c("1","2","3","4","5", "6", "7", "8"))

atg@meta.data$clusters <- factor(atg@meta.data$clusters, 
                                 levels = c("1","2","3","4","5", "6", "7", "8", "9", "10", "11"))

Idents(atg) <- atg$clusters

cols = c("#7bc07d", "#3467a6", "#7BD3EA", "#f4bb84", "#FFF78A", "#de1779", "#bbabcf", "#6c68a7")

DimPlot(atg, label = T, cols = cols)
SpatialDimPlot(atg, pt.size.factor = 2, stroke =NA, label =F, images = "slice1") + 
  scale_fill_manual(values = cols) + theme(legend.position = "none")
SpatialDimPlot(atg, pt.size.factor = 2, stroke =NA, label =F, images = "slice1.2") + 
  scale_fill_manual(values = cols) + theme(legend.position = "none")

ggplot(atg@meta.data, aes(atg@meta.data$orig.ident, fill = atg@meta.data$clusters)) +
  geom_bar(color = "black", width = 0.5)+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.title = element_blank())+
  theme(axis.text.y = element_text(size= 9, color = "black"), 
        axis.text.x = element_text(size= 9,  vjust =0.5),
        legend.text = element_text(size =11, lineheight =5)) +scale_fill_manual(values = cols) +
  ggplot(atg@meta.data, aes(atg@meta.data$orig.ident, fill = atg@meta.data$clusters)) +
  geom_bar(position = "fill", color = "black", width = 0.5)+
  theme_bw()+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.title = element_blank())+
  theme(axis.text.y = element_text(size= 9, color = "black"), 
        axis.text.x = element_text(size= 9,  vjust =0.5),
        legend.text = element_text(size =11, lineheight =5)) +scale_fill_manual(values = cols)


# Deconvolution Heatmap  -----------------------------------------------------------
# deconvolution result from Cell2location.py

cellscore <- read.csv("cell2location_adata_result.csv")
rownames(cellscore) <- cellscore$X
new_row_names <- gsub("-ATG7_(WT|KO)$", "", rownames(cellscore))
rownames(cellscore) <- new_row_names

# merge with atg data
df <- atg@meta.data
data <- merge(df, cellscore, by = "row.names", all.x = TRUE)
rownames(data) <- data$Row.names
head(data)

atg2 <- AddMetaData(atg, data)

df <- atg2@meta.data

cellscore_data <- df %>% dplyr::select("clusters", "Epithelial_cells", "Endothelial_cells", "Stromal_cells",
                                       "Alv_Mf", "Int_Mf",  "Neutrophil",
                                       "Monocytes", "DC", "Basophil",
                                       "B_cells",  "T_cells", "NK_cells", "Mesothelial_cells")
cluster_means <- aggregate(. ~ clusters, cellscore_data, mean)
rownames(cluster_means) <- cluster_means$clusters
cluster_means$clusters <- NULL

cluster_means <- t(cluster_means)
cluster_means

# pheatmap
library(RColorBrewer)

display.brewer.pal(7,"RdYlBu")

pheatmap(cluster_means, cluster_rows = FALSE,
         cluster_cols = FALSE, scale = "row",
         display_numbers = TRUE, fontsize = 13)


# Pathway enriched analysis -----------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(msigdbr)
library(fgsea)
library(enrichR)

# DEG 
atg_markers <- FindMarkers(atg, ident.1 = "ATG7_WT", 
                           group.by ="orig.ident", test.use = "MAST")

gene_list <- list()
gene_list$Atg7_WT <- atg_markers %>% filter(avg_log2FC >= 1.0 &  avg_log2FC <= 1.0 & p_val_adj < 0.05)

gene_list$Atg7_CKO <- atg_markers %>% filter(avg_log2FC < -1.0 &  avg_log2FC >= -1.0 &  p_val_adj < 0.05)

id_list <- list()
for (i in seq_along(gene_list)){
  gene_symbols <- as.character(rownames(gene_list[[i]]))
  tmp_id <- mapIds(org.Mm.eg.db, gene_symbols, 'ENTREZID', 'SYMBOL') 
  tmp_id <- tmp_id[!is.na(tmp_id)]
  id_list <- append(id_list, list(tmp_id))
}

names(id_list) <- names(gene_list)
pathway_GO <- compareCluster(geneClusters = id_list, fun='enrichGO', OrgDb='org.Mm.eg.db')

dotplot(pathway_GO, showCategory=8)



# Addmodule scores --------------------------------------------------------
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringi)
library(data.table)
library(msigdbr)

# Load DB (Example: Inflammatory Pathway; Apply the Same Code for Other Pathways)
inpGS <- data.table(msigdbr(species = "Mus musculus", category = "H"))

inpGS$gs_name <- gsub("HALLMARK_", "", inpGS$gs_name)
inpGS$gs_name <- gsub("_SIGNALING", "", inpGS$gs_name)
inpGS <- split(inpGS$gene_symbol, inpGS$gs_name)

## Run addmodulescore
p <- AddModuleScore(other, features = infla_list, name = "HALLMARK")
colnames(p@meta.data)[grep("HALLMARK", colnames(p@meta.data))] <- names(infla_list)

p5 <- subset(p, subset = clusters == "5")

## Plot
p1 <- FeaturePlot(p, reduction = "umap", pt.size = 0.1,
                  features = names(inpGS), order = TRUE, 
                  split.by ="sample") &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()

p2 <- SpatialFeaturePlot(p, features = names(infla_list), pt.size.factor = 2, stroke =NA)
p3 <- SpatialFeaturePlot(p5, features = names(infla_list), pt.size.factor = 2, stroke =NA)

## cluster 5
p5 <- subset(p, subset = clusters == "5")

d1 <- p5@meta.data
comparisons_list <- list(
  c("ATG7_KO", "ATG7_WT"))

ggplot(d1, aes(x=orig.ident, y=avg, fill = orig.ident)) + 
  geom_boxplot()+ theme_bw() + theme(legend.position = "none") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size =14)) +
  theme(axis.text.x = element_text(size =12, vjust = 0.5, color = "black")) +
  labs(x = "condtion", y = "Inflammatory response") +
  scale_fill_manual(values = c("ATG7_WT" = "#CCCCCC", "ATG7_KO" = "#FFC0CB"))+
  ggpubr::stat_compare_means(method = "t.test", label = "p.format", hide.ns = T, 
                             comparisons = list(c("Atg7 WT", "ATG7_KO")))















