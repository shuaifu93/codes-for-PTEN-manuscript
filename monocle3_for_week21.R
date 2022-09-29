library(Seurat)
library(monocle3)
library(patchwork)
library(cowplot)
library(SeuratWrappers)
organoid.combined.cca.sct_1_rename <- readRDS("w21_brain_2017_sct_cca_200_30_rename.rds")
DimPlot(organoid.combined.cca.sct_1_rename)

seurat_cleaned <- subset(organoid.combined.cca.sct_1_rename, idents = c("Neuron", "oRG", "EN_deep",
                                                                        "EN_upper", "tRG", "IPC", "M_phase",
                                                                        "S_phase"))
Idents(seurat_cleaned) <- "CellType"
Idents(seurat_cleaned) <- "orig.ident"
DimPlot(seurat_cleaned)
cds <- as.cell_data_set(seurat_cleaned, assay = "RNA")
cds <- cluster_cells(cds)

rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name


cds <- estimate_size_factors(cds)
cds <- learn_graph(cds, verbose = T, use_partition = F)
p7<- plot_cells(cds,
                color_cells_by = 'CellType',
                label_groups_by_cluster=F,
                label_leaves=F,
                label_branch_points=F,
                graph_label_size=4)
cds <- order_cells(cds)
p8<- plot_cells(cds,
                color_cells_by = "pseudotime",
                label_cell_groups=F,
                label_leaves=FALSE,
                label_branch_points=F,
                graph_label_size=4)
plot_grid(p7/p8)

saveRDS(cds, "cds_for_pseudotime_tRG_cycling_sox2_root_w21.rds" )
a<- as.Seurat(cds, assay = NULL)  #this works great now.  https://github.com/satijalab/seurat/issues/3746

Idents(a) <- "CellType"
DimPlot(a, reduction = "UMAP", ncol=3, label = T, order = T)
RidgePlot(a, features = c("monocle3_pseudotime"),  sort = T)#Ridg
Idents(a) <- "orig.ident"
cd<- subset(x = a, idents = c( 'correction.w21', 'Apex.w21'))
ab<- subset(x = a, idents = c( 'chap.w21', 'induction.w21'))
abc<- subset(ab, downsample = 3221)

cde<- subset(cd, downsample = 3407)
Idents(abc)<- "CellType"
p11<- RidgePlot(abc, features = c("monocle3_pseudotime"),  idents = "EN_upper", group.by = "orig.ident", cols = c('blue', 'red')
               , combine = F)
p12<- RidgePlot(abc, features = c("monocle3_pseudotime"),  idents = "EN_deep", group.by = "orig.ident", cols = c('blue', 'red'))
plot_grid(p10+p12)
plot_grid(p9+p11)

Idents(cde)<- "CellType"
p9<- RidgePlot(cde, features = c("monocle3_pseudotime"),  idents = "EN_upper", group.by = "orig.ident", cols = c('blue', 'red')
)#RidgePlot(integrated, features = c("monocle3_pseudotime"), group.by = "monocle3_clusters",  sort = T)
p10<- RidgePlot(cde, features = c("monocle3_pseudotime"),  idents = "EN_deep", group.by = "orig.ident", cols = c('blue', 'red'))#RidgePlot(integrated, features = c("monocle3_pseudotime"), group.by = "monocle3_clusters",  sort = T)
DimPlot(cde, group.by = "orig.ident", label = T)
plot_grid(p10+p9)


##########
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
plot_cells(cds, 
           genes="SOX2",
           show_trajectory_graph=T,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           reduction_method = "UMAP") 


