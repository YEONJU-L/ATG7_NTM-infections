
#### ATG7_Single-cell RNA-sequencing analyses script.R ####

# load library ------------------------------------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(harmony)
library(pheatmap)
library(ape)

wt.data <- Read10X(data.dir = "CR_24_15074_TS_C_SS3_1/Cell_matrix/")
wt <- CreateSeuratObject(counts = wt.data, project = "ATG7_WT")
dim(wt)

ko.data <- Read10X(data.dir = "CR_24_15075_TS_C_SS3_1/Cell_matrix/")
ko <- CreateSeuratObject(counts = ko.data, project = "ATG7_cKO")
dim(ko)

wt$barcodes <- rownames(wt@meta.data)
ko$barcodes <- rownames(ko@meta.data)

wt$orig.ident <- "ATG7_WT"
ko$orig.ident <- "ATG7_KO"

merge <- merge(wt, ko, add.cell.ids = c("WT", "cKO"))
merge@meta.data$orig.ident <- factor(merge@meta.data$orig.ident, levels = c("ATG7_WT", "ATG7_cKO"))


# QC ----------------------------------------------------------------------

VlnPlot(merge, c("nCount_RNA","nFeature_RNA"), cols = c("#66c2a5", "#fc8d62"), group.by = "orig.ident", pt.size =0)

ggplot(merge@meta.data, aes(x = orig.ident, y = `nFeature_RNA`)) +
  geom_boxplot() +  theme_minimal()+
  labs(title = "nFeature_RNA", y = "nCount") +
  theme(legend.position = "none", plot.title = element_text(size =15, hjust = 0.5),
        axis.text = element_text(size=11, color ="black"),
        axis.title = element_text(size=13))

summary(ko$nFeature_RNA)
summary(wt$nFeature_RNA)
summary(merge$nFeature_RNA)

# filtering low QC
merge$lowQC <- ifelse(merge$nFeature_RNA >200, "PASS", "FAIL")
table(merge$orig.ident, merge$lowQC)

FeatureScatter(merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
               group.by = "lowQC", cols = c("PASS" = "black", "FAIL" = "red")) + theme(legend.position = "none")

merge <- subset(merge, subset = lowQC == "PASS") # 32285 25040


# Processing ----------------------------------------------------------

merge <- JoinLayers(merge)
merge <- NormalizeData(merge)
merge <- FindVariableFeatures(merge, selection.method = "vst", 
                              nfeatures = 2000)
merge <- ScaleData (merge, features = rownames(merge))
merge <- RunPCA(merge, features = VariableFeatures(object = merge))

ElbowPlot(merge, ndims = 50)

merge <- merge %>% RunUMAP(dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.5)

DimPlot(merge, split.by = "orig.ident")

# Batch Correction
merge <- RunHarmony(merge, group.by.vars = "orig.ident")

merge <- merge %>% RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.5)

merge <- merge %>% FindClusters(resolution = 1.0)

DimPlot(merge,label =T)


# Annotation --------------------------------------------------------------

# myeloid markers
# https://azimuth.hubmapconsortium.org/references/#Mouse%20-%20Pansci
mast_cells <- c("Cpa3", "Mcpt1", "Tpsb2", "Mcpt2", "Mcpt9", "Olfr920", "Mrgprb2", "Cd200r3", "Mrgprb1", "Tpsab1")
baso <- c("Rasa3", "Dhx40", "Hgf", "Arhgap10", "Il6", "Nrgn", "Hdc", "Nbeal2", "Slc6a4", "Pim1") 
monocytes <- c("Hc", "Vcan", "Fgr", "Igsf6", "Pot1b", "RP23-349H13.2", "Sirpb1c", "Klra2", "Treml4", "Gpr141")
neutrophils <- c("Irg1", "S100a9", "Ncf2", "Clec4d", "Stfa2l1", "R3hdm4", "Trem1", "A530064D06Rik", "RP24-127D13.2", "Retnlg")
AM <-  c("Olr1", "Ncf1", "Cxcl2", "Fpr2", "Clec7a", "Cyb561a3", "RP24-443M8.2", "Abcg1", "Clec4n", "Fpr1") 
IM <- c("Sirpa", "Stab1", "Msr1", "Ccr5", "Nrros", "Marco", "Aoah", "Sdc3", "Apobec1", "Mrc1")

# annotation
Idents(merge) <- merge$seurat_clusters 

new.cluster.ids <- c("0" = "B_cell", "1" ="CD4+Tcell", "2" ="CD4+/8+Tcell",
                     "3" ="Cap", "4" ="Int_Macrophage", "5" ="Monocyte",
                     "6" ="CD8+Tcell", "7" ="Col13a1+fibroblast", "8"="Int_Macrophage", "9"="Cap", 
                     "10"="B_cell", "11"="Neutrophil", "12"="Alv_Macrophage", "13"="Cap",
                     "14"="Basophil", "16"="DC", "17"="Cap", "18"="DC", 
                     "19"='Monocyte', "20"="NK_cell", "21"="Col13a1+fibroblast", "22"="B_cell", "23"="Plasma_cell", "24"="AT2", 
                     "25"="AT1", "26"="Pericyte2", 
                     "27"="Col14a1+fibroblast", "28"="Pericyte1", "29"="Vein", "30"="Cap", 
                     "31"="Int_Macrophage", "33"="Treg_cell", "34"="Lymph", "35"="Alv_Macrophage", "36"="Ciliated")

merge <- RenameIdents(merge, new.cluster.ids)
merge$celltype <- merge@active.ident

merge@meta.data$celltype <- factor(merge@meta.data$celltype, 
                                   levels = c("AT1", "AT2","Ciliated",
                                              "Cap", "Vein", "Lymph",
                                              "Col13a1+fibroblast", "Col14a1+fibroblast", "Pericyte1", "Pericyte2",
                                              "Int_Macrophage", "Alv_Macrophage", "Neutrophil", "Basophil","Monocyte", "DC",
                                              "CD4+Tcell","CD8+Tcell","CD4+/8+Tcell", "Treg_cell","NK_cell" ,
                                              "B_cell", "Plasma_cell"))

Idents(merge) <- merge$celltype

cols <- c(
  "#EAC7C7", "#CA8A8B", "#CDC2AE", 
  "#BBD6B8", "#A0C3D2", "#FF8AAE", 
  "#537188", "#DAEAF1", "#FFA6C1", "#D0C9C0", 
  "#3467a6","#de1779","#FFF78A", "#bbabcf", "#7bc07d",
  "#FDC5F5", "#E5D4FF", "#D0A2F7",
  "#96DED1", "#FFC4C4", "#A7BED3", "#F4B183"
)

DimPlot(merge, label =T, cols = cols)


# Addmodule scores ---------------------------------------------------------

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringi)
library(data.table)
library(msigdbr)
library(ggpubr)

## Load DB 
inpGS <- data.table(msigdbr(species = "Mus musculus",  category = "C2", subcategory = "CP:REACTOME"))
inpGS_1 <- data.table(msigdbr(species = "Mus musculus", category = "H"))
inpGS2 <- data.table(msigdbr(species = "Mus musculus",  category = "C5", subcategory = "GO:BP"))
inpGS3 <- data.table(msigdbr(species = "Mus musculus",  category = "C2", subcategory = "CP:KEGG"))

inpGS <- split(inpGS$gene_symbol, inpGS$gs_name)
inpGS_1 <- split(inpGS_1$gene_symbol, inpGS_1$gs_name)
inpGS2 <- split(inpGS2$gene_symbol, inpGS2$gs_name)
inpGS3 <- split(inpGS3$gene_symbol, inpGS3$gs_name)

path_list <- list()
path_list$INFLAMMATORY_RESPONSE <- inpGS_1$HALLMARK_INFLAMMATORY_RESPONSE
path_list$REACTIVE_OXYGEN_SPECIES_PATHWAY  <- inpGS_1$HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY

path_list$REACTOME_PYROPTOSIS  <- inpGS$REACTOME_PYROPTOSIS
path_list$REACTOME_REGULATED_NECROSIS  <- inpGS$REACTOME_REGULATED_NECROSIS
path_list$APOPTOSIS  <- inpGS3$KEGG_APOPTOSIS

path_list$GSDME <- "Gsdme"


## Run addmodulescore
p <- AddModuleScore(merge, features = path_list, name = "HALLMARK")
colnames(p@meta.data)[grep("HALLMARK", colnames(p@meta.data))] <- names(path_list)

p2 <- subset(merge, subset = celltype %in% c("Int_Macrophage", "Alv_Macrophage", "Neutrophil", "Monocyte", "DC")))

d1 <- p2@meta.data
df <- d1 %>% dplyr::mutate(avg = rowMeans(dplyr::select(.,c(names(path_list)))))
head(df)

my_comparisons <- list(c("ATG7_WT", "ATG7_cKO"))

ggplot(df, aes(x = orig.ident, y = avg, fill = orig.ident, color = orig.ident)) +   
  geom_violin(aes(fill = orig.ident), trim = T) +
  stat_summary(fun = median, geom = "point", shape = 95, size = 10, color = "#555555") + # meidan line
  theme_classic() + 
  theme(legend.position = "none",  plot.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 30, vjust = 0.5, hjust = 0.5, size = 11, color = "black"),
        axis.text.y = element_text(size = 11, color = "black"),
        axis.title.y = element_text(size = 13),
        axis.title.x = element_blank()) +  
  labs(y = "[Module Score]", title = "ROS") +  # necrosis
  scale_fill_manual(values = c("ATG7_WT" = "#CCCCCC", "ATG7_cKO" = "#FFC0CB")) +
  scale_color_manual(values = c("ATG7_WT" = "#BEBEBE", "ATG7_cKO" = "#FFC0CB")) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox", label = "p.format", size = 4.5)

FeaturePlot(p, reduction = "umap", pt.size = 0.1,
            features = "GSDME", order = TRUE, 
            split.by ="orig.ident") &
  scale_color_distiller(palette = "RdYlBu") & coord_fixed()

df_k <- subset(df, subset = orig.ident == "ATG7_cKO")

# Statistical Test
# Shapiro-Wilk Test
by(df_k$avg, df_k$celltype, shapiro.test)

# Levene's Test
library(car)
leveneTest(avg ~ celltype, data = df_k)

# ANOVA 
anova_result <- aov(avg ~ celltype, data = df_k)
summary(anova_result)

# Tukey HSD
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

