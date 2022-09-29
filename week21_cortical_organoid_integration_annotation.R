library(Seurat)
library(ggplot2)
library(sctransform)
library(presto)
require(data.table)
library(dplyr)
library(scCustomize)

# Loading Cellranger processed data
chap.w21.data <- Read10X_h5('chap_week21_filtered_feature_bc_matrix.h5')
chap_P.w21.data <- Read10X_h5('chap_p_week21_filtered_feature_bc_matrix.h5')
C_KO.w21.data <- Read10X_h5('C_KO_week21_filtered_feature_bc_matrix.h5')
A_PC.w21.data <- Read10X_h5('A_PC_week21_filtered_feature_bc_matrix.h5')
Apex.w21.data <- Read10X_h5('Apex_week21_filtered_feature_bc_matrix.h5')
A_KO.w21.data <- Read10X_h5('A_KO_week21_filtered_feature_bc_matrix.h5')

Apex.w21 <- CreateSeuratObject(counts = Apex.w21.data, project = "Apex WT/I135L.w21", min.cells = 3, min.features = 200)
Apex.w21[["percent.mt"]] <- PercentageFeatureSet(Apex.w21, pattern = "^MT[-\\.]")
VlnPlot(Apex.w10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
Apex.w21 <- subset(Apex.w21, subset = nFeature_RNA > 200 & percent.mt < 30)

correction.w21 <- CreateSeuratObject(counts = A_PC.w21.data, project = "Apex WT/WT.w21", min.cells = 3, min.features = 200)
correction.w21[["percent.mt"]] <- PercentageFeatureSet(correction.w21, pattern = "^MT[-\\.]")
VlnPlot(correction.w21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
correction.w21 <- subset(correction.w21, subset = nFeature_RNA > 200 & percent.mt < 30)

chap.w21 <- CreateSeuratObject(counts = chap.w21.data, project = "Chap WT/WT.w21", min.cells = 3, min.features = 200)
chap.w21[["percent.mt"]] <- PercentageFeatureSet(chap.w21, pattern = "^MT[-\\.]")
VlnPlot(chap.w21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
chap.w21 <- subset(chap.w21, subset = nFeature_RNA > 200 & percent.mt < 30)

induction.w21 <- CreateSeuratObject(counts = chap_P.w21.data, project = "Chap WT/I135L.w21", min.cells = 3, min.features = 200)
induction.w21[["percent.mt"]] <- PercentageFeatureSet(induction.w21, pattern = "^MT[-\\.]")
VlnPlot(induction.w21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
induction.w21 <- subset(induction.w21, subset = nFeature_RNA > 200 & percent.mt < 30)

C_KO.w21 <- CreateSeuratObject(counts = C_KO.w21.data, project = "Chap KO/KO.w21", min.cells = 3, min.features = 200)
C_KO.w21[["percent.mt"]] <- PercentageFeatureSet(C_KO.w21, pattern = "^MT[-\\.]")
VlnPlot(C_KO.w21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
C_KO.w21 <- subset(C_KO.w21, subset = nFeature_RNA > 200 & percent.mt < 30)

A_KO.w21 <- CreateSeuratObject(counts = A_KO.w21.data, project = "Apex KO/KO.w21", min.cells = 3, min.features = 200)
A_KO.w21[["percent.mt"]] <- PercentageFeatureSet(A_KO.w21, pattern = "^MT[-\\.]")
VlnPlot(A_KO.w21, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
A_KO.w21 <- subset(A_KO.w21, subset = nFeature_RNA > 200 & percent.mt < 30)

#########
correction.w21 <- SCTransform(correction.w21, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)

Apex.w21 <- SCTransform(Apex.w21, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)

chap.w21 <- SCTransform(chap.w21, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)

induction.w21 <- SCTransform(induction.w21, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)

C_KO.w21 <- SCTransform(C_KO.w21, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)

A_KO.w21 <- SCTransform(A_KO.w21, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)


brain_2017 <- readRDS("brain_2017_rename_SCT_v2.rds")


s_list = list(brain_2017, chap.w21, induction.w21, C_KO.w21, correction.w21, Apex.w21, A_KO.w21 )

features <- SelectIntegrationFeatures(object.list = s_list, nfeatures = 3000)
s_list <- PrepSCTIntegration(object.list = s_list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = s_list, normalization.method = "SCT", 
                                  anchor.features = features)
organoid.combined.cca.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
saveRDS(organoid.combined.cca.sct, "brain_2017_w21_all.sct.cca_200_30_5.8.22.rds")

organoid.combined.cca.sct <- readRDS("brain_2017_w21_all.sct.cca_200_30_5.8.22.rds")
organoid.combined.cca.sct <- RunPCA(organoid.combined.cca.sct, verbose = T)

# Plot PCA
PCAPlot(organoid.combined.cca.sct,
        split.by = "orig.ident")
organoid.combined.cca.sct <- RunUMAP(organoid.combined.cca.sct, reduction = "pca", dims = 1:30)
#organoid.combined.cca.sct <- RunTSNE(organoid.combined.cca.sct, reduction = "pca", dims = 1:30)

organoid.combined.cca.sct <- FindNeighbors(organoid.combined.cca.sct, dims = 1:30)
organoid.combined.cca.sct_1 <- FindClusters(organoid.combined.cca.sct, resolution = 1.5)

DimPlot(organoid.combined.cca.sct, reduction = 'umap', label = F, repel = T, split.by = "orig.ident") + NoLegend()

DimPlot(organoid.combined.cca.sct_1, reduction = 'umap', label = T, repel = T)
organoid.combined.cca.sct_1$orig.ident <- factor(organoid.combined.cca.sct_1$orig.ident, levels = c("Chap WT/WT.w21", "Chap WT/I135L.w21", "Chap KO/KO.w21", "Apex WT/WT.w21", "Apex WT/I135L.w21", "Apex KO/KO.w21", "brain_2017"))
DimPlot(organoid.combined.cca.sct_1, reduction = 'umap', label = T, repel = T)

Idents(organoid.combined.cca.sct_1) <- "seurat_clusters"

DefaultAssay(object = organoid.combined.cca.sct_1) <- "SCT"

organoid.combined.cca.sct_1<- PrepSCTFindMarkers(organoid.combined.cca.sct_1)
# Example for extracting markers genes in cluster 0, 
# this needs to be done for all clusters separately.
markers_0 <- FindMarkers(organoid.combined.cca.sct_1, assay = "SCT", 
                         only.pos = TRUE, 
                         ident.1 = "0", ident.2 = NULL)
saveRDS(markers_0, "markers_0_w21_prepSCT_findmarkers.rds")
library(dplyr)
top200_ <- 
  markers_0 %>%
  top_n(n = 200, wt = avg_log2FC)  #
write.csv(top200_,"top200_cellranger_week21_200_30_cca.sct_findmarkers_0_8.7.22.csv", row.names = T)
new_ident_1 <- setNames(c("U1", #0 
                          "oRG",  #1  
                          "EN_deep",#2 
                          "oRG", #3 
                          "U2", #4 
                          "EN_upper",  #5 
                          "tRG", #6 
                          "EN_upper", #7
                          "Astrocyte", #8
                          "IN", #9
                          "Cajal-Retzius", #10 
                          "IPC", #11, 
                          "EN_upper", #12
                          "Glyc", #13 
                          "U3",#14 
                          "Glyc",  #15 
                          "U4", #16
                          "Endothelial", #17 
                          "S/M_phase", #18 
                          "U5", #19 
                          "IN", #20 
                          "U6", #21
                          "IN", #22
                          "EN_deep", #23
                          "S_phase", #24
                          "U7", #25
                          "Cortical_hem", #26
                          "M_phase", #27
                          "U8", #28
                          "OPC", #29
                          "U9", #30
                          "U10", #31
                          "U11", #32
                          "Microglia", #33 
                          "S/M_phase", #34
                          "IN", #35
                          "Choroid", #36
                          "Oligodendrocyte" #37
),

levels(organoid.combined.cca.sct_1_rename))
organoid.combined.cca.sct_1_rename <- RenameIdents(organoid.combined.cca.sct_1_rename, new_ident_1)

organoid.combined.cca.sct_1_rename$CellType <- Idents(organoid.combined.cca.sct_1_rename)

saveRDS(organoid.combined.cca.sct_1_rename, "w21_brain_2017_sct_cca_200_30_rename_findmarkers_8.7.22.rds")
