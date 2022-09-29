library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(fgsea)
library(msigdbr)
library(enrichplot)
# choose your msigdb collection of interest
hs_gsea_c5 <- msigdbr(species = "Homo sapiens", 
                      category = "C5",
                      subcategory = "BP"
) %>% 
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 
hs_gsea_c5<- hs_gsea_c5 %>% mutate(gs_name = str_replace_all(gs_name, 'GOBP_', ''))

organoid.combined.cca.sct_1_rename<- readRDS("w10_brain_2017_sct_cca_200_30_rename_8.6.22.rds")
DimPlot(organoid.combined.cca.sct_1_rename, reduction = "umap", label = T, repel = T)
Idents(object = organoid.combined.cca.sct_1_rename) <- "CellType"
### Example for generating rnk files for oRG cell type, Chap WT/I135L vs Chap WT/WT at week 10. 
ABC<- subset(x = organoid.combined.cca.sct_1_rename, idents = c( "oRG"))

Idents(object = ABC) <- "orig.ident"
ABC<- PrepSCTFindMarkers(ABC, assay = "SCT", verbose = T)

Chap_p_vs_chap <- FindMarkers(ABC, 
                              logfc.threshold = -Inf,  # equivalent to using 0
                              assay = "SCT",
                              min.pct = 0.1,
                              ident.1 = "Chap WT/I135L.w10", ident.2 = "Chap WT/WT.w10",
                              test.use = "wilcox"
                              )

Chap_p_vs_chap.df <- Chap_p_vs_chap %>%
  as_tibble(rownames = "geneID")

#create ranks file
ranks_chap_p.chap = Chap_p_vs_chap.df$avg_log2FC
ranks_chap_p.chap <- cbind(Chap_p_vs_chap.df$geneID, ranks_chap_p.chap) 
colnames(ranks_chap_p.chap) <- c("GeneName","rank")
head(ranks_chap_p.chap)
#sort ranks in decreasing order
ranks_chap_p.chap <- ranks_chap_p.chap[order(as.numeric(ranks_chap_p.chap[,2]),decreasing = TRUE),]
write.table(ranks_chap_p.chap,file="ranks_oRG.induction_vs_chap_cellranger_w10_sct_cca_wilcox_lfc.rnk",col.name = TRUE, sep="\t", row.names = FALSE, quote = FALSE)

