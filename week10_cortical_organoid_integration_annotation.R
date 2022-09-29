library(Seurat)
library(ggplot2)
library(sctransform)
library(presto)
require(data.table)
library(dplyr)
library(scCustomize)
# Loading Cellranger processed data

chap.w10.data <- Read10X_h5('chap_week10_R2_filtered_feature_bc_matrix.h5')
chap_P.w10.data <- Read10X_h5('chap_p_week10_R2_filtered_feature_bc_matrix.h5')
C_KO.w10.data <- Read10X_h5('c_ko_week10_R2_filtered_feature_bc_matrix.h5')
A_PC.w10.data <- Read10X_h5('A_PC_week10_R2_filtered_feature_bc_matrix.h5')
Apex.w10.data <- Read10X_h5('Apex_week10_r2_filtered_feature_bc_matrix.h5')
A_KO.w10.data <- Read10X_h5('A_KO_week10_r2_filtered_feature_bc_matrix.h5')


Apex.w10 <- CreateSeuratObject(counts = Apex.w10.data, project = "Apex WT/I135L.w10", min.cells = 3, min.features = 200)
Apex.w10[["percent.mt"]] <- PercentageFeatureSet(Apex.w10, pattern = "^MT[-\\.]")
Apex.w10 <- subset(Apex.w10, subset = nFeature_RNA > 200 & percent.mt < 30)
VlnPlot(Apex.w10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

correction.w10 <- CreateSeuratObject(counts = A_PC.w10.data, project = "Apex WT/WT.w10", min.cells = 3, min.features = 200)
correction.w10[["percent.mt"]] <- PercentageFeatureSet(correction.w10, pattern = "^MT[-\\.]")
VlnPlot(correction.w10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
correction.w10 <- subset(correction.w10, subset = nFeature_RNA > 200 & percent.mt < 30)

chap.w10 <- CreateSeuratObject(counts = chap.w10.data, project = "Chap WT/WT.w10", min.cells = 3, min.features = 200)
chap.w10[["percent.mt"]] <- PercentageFeatureSet(chap.w10, pattern = "^MT[-\\.]")
VlnPlot(chap.w10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
chap.w10 <- subset(chap.w10, subset = nFeature_RNA > 200 & percent.mt < 30)

induction.w10 <- CreateSeuratObject(counts = chap_P.w10.data, project = "Chap WT/I135L.w10", min.cells = 3, min.features = 200)
induction.w10[["percent.mt"]] <- PercentageFeatureSet(induction.w10, pattern = "^MT[-\\.]")
VlnPlot(induction.w10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
induction.w10 <- subset(induction.w10, subset = nFeature_RNA > 200 & percent.mt < 30)

C_KO.w10 <- CreateSeuratObject(counts = C_KO.w10.data, project = "Chap KO/KO.w10", min.cells = 3, min.features = 200)
C_KO.w10[["percent.mt"]] <- PercentageFeatureSet(C_KO.w10, pattern = "^MT[-\\.]")
VlnPlot(C_KO.w10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
C_KO.w10 <- subset(C_KO.w10, subset = nFeature_RNA > 200 & percent.mt < 30)

A_KO.w10 <- CreateSeuratObject(counts = A_KO.w10.data, project = "Apex KO/KO.w10", min.cells = 3, min.features = 200)
A_KO.w10[["percent.mt"]] <- PercentageFeatureSet(A_KO.w10, pattern = "^MT[-\\.]")
VlnPlot(A_KO.w10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
A_KO.w10 <- subset(A_KO.w10, subset = nFeature_RNA > 200 & percent.mt < 30)

#########
correction.w10 <- SCTransform(correction.w10, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)

Apex.w10 <- SCTransform(Apex.w10, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)

chap.w10 <- SCTransform(chap.w10, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)

induction.w10 <- SCTransform(induction.w10, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)

C_KO.w10 <- SCTransform(C_KO.w10, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)

A_KO.w10 <- SCTransform(A_KO.w10, vst.flavor = "v2", verbose = T) %>%
  RunPCA(npcs = 30, verbose = T)


brain_2017 <- readRDS("brain_2017_rename_SCT_v2.rds")

s_list = list(brain_2017, chap.w10, induction.w10, C_KO.w10, correction.w10, Apex.w10, A_KO.w10 )

features <- SelectIntegrationFeatures(object.list = s_list, nfeatures = 3000)
s_list <- PrepSCTIntegration(object.list = s_list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = s_list, normalization.method = "SCT", 
                                  anchor.features = features)
organoid.combined.cca.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
saveRDS(organoid.combined.cca.sct, "brain_2017_w10_all.sct.cca_200_30_5.8.22.rds")

organoid.combined.cca.sct <- readRDS("brain_2017_w10_all.sct.cca_200_30_5.8.22.rds")
organoid.combined.cca.sct <- RunPCA(organoid.combined.cca.sct, verbose = T)

# Plot PCA
PCAPlot(organoid.combined.cca.sct,
        split.by = "orig.ident")
organoid.combined.cca.sct <- RunUMAP(organoid.combined.cca.sct, reduction = "pca", dims = 1:30)
#organoid.combined.cca.sct <- RunTSNE(organoid.combined.cca.sct, reduction = "pca", dims = 1:30)

organoid.combined.cca.sct <- FindNeighbors(organoid.combined.cca.sct, dims = 1:30)
organoid.combined.cca.sct_1 <- FindClusters(organoid.combined.cca.sct, resolution = 1.5)

DimPlot(organoid.combined.cca.sct, reduction = 'umap', label = F, repel = T, split.by = "orig.ident", ncol = 3) 
DimPlot(organoid.combined.cca.sct, reduction = 'umap', label = F, repel = T, split.by = "orig.ident", ncol =3) + NoLegend()
DimPlot(organoid.combined.cca.sct_1, reduction = 'umap', label = T, repel = T)

organoid.combined.cca.sct_1$orig.ident <- factor(organoid.combined.cca.sct_1$orig.ident, levels = c("Chap WT/WT.w10", "Chap WT/I135L.w10", "Chap KO/KO.w10", "Apex WT/WT.w10", "Apex WT/I135L.w10", "Apex KO/KO.w10", "brain_2017"))

Idents(organoid.combined.cca.sct_1) <- "seurat_clusters"

DefaultAssay(object = organoid.combined.cca.sct_1) <- "SCT"
organoid.combined.cca.sct_1<- PrepSCTFindMarkers(organoid.combined.cca.sct_1)
# Example for extracting markers genes in cluster 0, 
# this needs to be done for all clusters separately.
markers_0 <- FindMarkers(organoid.combined.cca.sct_1, assay = "SCT", 
                          only.pos = TRUE, 
                          ident.1 = "0", ident.2 = NULL)
saveRDS(markers_0, "markers_0_w10_prepSCT_findmarkers.rds")
library(dplyr)
top200_ <- 
  markers_0 %>%
  top_n(n = 200, wt = avg_log2FC)  #
write.csv(top200_,"top200_cellranger_week10_200_30_cca.sct_findmarkers_0_8.7.22.csv", row.names = T)
# annotating the cell types for each cluster based on marker gene expression. 
new_ident_1 <- setNames(c("EN_upper", #0 
                          "Glyc",  #1  
                          "U1",#2
                          "Cortical_hem", #3
                          "IN", #4
                          "tRG",  #5  
                          "IPC_EN_newborn", #6
                          "EN_upper", #7
                          "EN_upper", #8
                          "vRG", #9  
                          "IN", #10 , 
                          "oRG", #11, 
                          "Cajal-Retzius", #12 
                          "U2", #13
                          "Microglia",#14 
                          "tRG",  #15 
                          "IPC", #16  
                          "U3", #17
                          "U4", #18
                          "U5", #19
                          "EN_deep", #20 
                          "IN", #21
                          "U6", #22
                          "U7", #23
                          "M_phase", #24
                          "U8", #25
                          "U9", #26
                          "S/M_phase", #27 
                          "S_phase", #28
                          "IN", #29
                          "M_phase", #30
                          "vRG", #31
                          "Choroid", #32 
                          "U10", #33 
                          "IN" #34
                          
),
levels(organoid.combined.cca.sct_1_rename))
organoid.combined.cca.sct_1_rename <- RenameIdents(organoid.combined.cca.sct_1_rename, new_ident_1)
organoid.combined.cca.sct_1_rename$CellType <- Idents(organoid.combined.cca.sct_1_rename)

Idents(organoid.combined.cca.sct_1_rename) <- "orig.ident"
DimPlot(organoid.combined.cca.sct_1_rename, reduction = "umap", split.by = "orig.ident", ncol=3, label = F, cols = c("gray","gray", "gray", "gray", "gray","gray","gray")) + NoLegend()

DimPlot(organoid.combined.cca.sct_1_rename, reduction = 'umap', label = T, repel = T, label.size  = 4.5)


saveRDS(organoid.combined.cca.sct_1_rename, "w10_brain_2017_sct_cca_200_30_rename_8.6.22.rds")
