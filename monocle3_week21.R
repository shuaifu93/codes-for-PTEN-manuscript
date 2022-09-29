library(Seurat)
library(monocle3)
library(patchwork)
library(cowplot)
library(SeuratWrappers)
organoid.combined.cca.sct_1_rename <- readRDS("w21_brain_2017_sct_cca_200_30_rename_findmarkers_8.7.22.rds")
DimPlot(organoid.combined.cca.sct_1_rename, label = T) 

seurat_cleaned <- subset(organoid.combined.cca.sct_1_rename, 
                         idents = c( "oRG", "IPC","tRG", "M_phase", "S_phase", 
                       "S/M_phase","Cajal-Retzius", "EN_deep", "EN_upper"))
                                                                         
                                                                        

Idents(seurat_cleaned) <- "CellType"

DimPlot(seurat_cleaned, label = T)
cds <- as.cell_data_set(seurat_cleaned, assay = "RNA")
cds <- cluster_cells(cds, resolution=1e-5)

rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
cds <- estimate_size_factors(cds)

plot_cells(cds,
           genes="SOX2",
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)
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



saveRDS(cds, "cds_for_pseudotime_sox2_oRG_root_w21_8.9.22.rds" )
cds <- readRDS("cds_for_pseudotime_sox2_oRG_root_w21_8.9.22.rds")
a<- as.Seurat(cds, assay = NULL)  #this works great now.  https://github.com/satijalab/seurat/issues/3746

Idents(a) <- "CellType"
DimPlot(a, reduction = "UMAP", ncol=3, label = T, order = T)
RidgePlot(a, features = c("monocle3_pseudotime"),  sort = T)#Ridg
Idents(a) <- "orig.ident"
cd<- subset(x = a, idents = c( 'correction.w21', 'Apex.w21'))
ab<- subset(x = a, idents = c( 'chap.w21', 'induction.w21'))
abc<- subset(ab, downsample = 2808)

cde<- subset(cd, downsample = 3047)
Idents(abc)<- "CellType"
p11<- RidgePlot(abc, features = c("monocle3_pseudotime"),  idents = "EN_upper", group.by = "orig.ident", cols = c('blue', 'red')
               , combine = F)
p12<- RidgePlot(abc, features = c("monocle3_pseudotime"),  idents = "EN_deep", group.by = "orig.ident", cols = c('blue', 'red'))


Idents(cde)<- "CellType"
p9<- RidgePlot(cde, features = c("monocle3_pseudotime"),  idents = "EN_upper", group.by = "orig.ident", cols = c('blue', 'red')
)
p10<- RidgePlot(cde, features = c("monocle3_pseudotime"),  idents = "EN_deep", group.by = "orig.ident", cols = c('blue', 'red'))
DimPlot(cde, group.by = "orig.ident", label = T)
plot_grid(p9+p11)
plot_grid(p10+p12)


pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=5)

saveRDS(pr_graph_test_res, "pr_graph_test_res_w10_8.10.22.rds")
pr_graph_test_res <- readRDS("pr_graph_test_res_w10.rds")

pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
gene_list_1 <- pr_deg_ids[1:10]
hgenes <- plot_cells(cds, 
                     genes = gene_list_1,
                     show_trajectory_graph=T,
                     label_cell_groups=FALSE,
                     label_leaves=FALSE)



##########
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
plot_cells(cds, 
           genes="DCX",
           show_trajectory_graph=T,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           reduction_method = "UMAP") 

tail(cds@colData$Size_Factor)

plot_cells(cds_subset1,
           label_cell_groups=T,
           color_cells_by = "pseudotime",
           show_trajectory_graph=F)
plot_cells(cds,
           label_cell_groups=T,
           genes=c("PTEN","TIGIT"),
           show_trajectory_graph=F)
# This solved "no" expression issue. 
cds <- estimate_size_factors(cds)


rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(so[["RNA"]])
cds@rowRanges


cds@preprocess_aux@listData$gene_loadings
row.names(cds@preprocess_aux$gene_loadings)
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
gene_module_df = monocle3::find_gene_modules(cds[var.genes,], # cds[var.genes, ]
                                             # resolution=1e-2, 
                                             reduction_method = "UMAP",
                                             max_components = 3,
                                             cores=nCores,
                                             umap.fast_sgd = F, # Slower but deterministic results
                                             random_seed = 2019)





rownames(cds@preprocess_aux$gene_loadings)
cds@preprocess_aux
cds@preprocess_aux
rownames(cds)
rowData(cds)$gene_short_name <- rownames(cds)





cds <- readRDS("cds_for_pseudotime_vRG_cycling_root_w10.rds")
cds_subset <- cds[,colData(cds)$orig.ident %in% c("correction.w10", "Apex.w10")] 

my_genes<- c("SOX2", "EOMES", "TBR1", "HOPX", "SATB2")
my_genes<- c("EOMES")

cds_subset1 <- cds[my_genes,]

plot_genes_in_pseudotime(cds_subset1,
                         color_cells_by="pseudotime", 
                         min_expr=0.5)

plot_in_pseudotime(cds_subset1, nrow = NULL, ncol = NULL, use_short_names = FALSE,
                   use_log_scale = FALSE, facet_wrap_scales = "fixed",
                   color_by = "cluster", y_lab = "Expression")
#plot_cell_trajectory(cds_subset)

cds_subset1 <- choose_cells(cds_subset)
cds_subset1<- preprocess_cds(cds_subset1)

tail(cds_subset1@colData$Size_Factor)

pr_graph_test_res <- graph_test(cds_subset1, neighbor_graph="knn", cores=8)
pr_deg_ids <- row.names(subset(pr_graph_test_res, morans_I > 0.01 & q_value < 0.05))

gene_module_df <- find_gene_modules(cds_subset1[pr_deg_ids,], resolution=1e-3)
rownames(cds_subset1@preprocess_aux$gene_loadings)

cell_group_df <- tibble::tibble(cell=row.names(colData(cds_subset1)), 
                                cell_group=colData(cds_subset1)$cell.type)
agg_mat <- aggregate_gene_expression(cds_subset1, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")







#plot_genes_branched_heatmap(cds_subset)
#plot_pseudotime_heatmap(cds_subset)

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)

modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(pr_graph_test_res, q_value == 0 & morans_I > 0.25))
genes

pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes;
#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  km = 2,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

print(htkm)
print(hthc)




##### this is how subsetting cds in monocle3. 
cds_subset <- cds[,colData(cds)$orig.ident %in% c("correction.w21", "apex.w21")] 

plot_cells(cds_subset,
           color_cells_by = "pseudotime",
           label_cell_groups=F,
           label_leaves=FALSE,
           label_branch_points=F,
           graph_label_size=4)
library(monocle3)
lung <- load_lung()
plot_cell_trajectory(lung, color_by = "Time")



