library(clusterProfiler)
library(enrichplot)
library(tidyverse)
library(readr)
# choose your msigdb collection of interest
hs_gsea_c5 <- msigdbr(species = "Homo sapiens", 
                      category = "C5",
                      subcategory = "BP"
) %>% 
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 
hs_gsea_c5<- hs_gsea_c5 %>% mutate(gs_name = str_replace_all(gs_name, 'GOBP_', ''))
#for  S_w10
S_w10 = read_tsv("ranks_S.apex_vs_correction_cellranger_w10_sct_cca_wilcox_lfc.rnk")
S_w10.gsea <- S_w10$rank
names(S_w10.gsea) <- as.character(S_w10$GeneName)
S_w10.gsea <- sort(S_w10.gsea, decreasing = T)
#for  S_w21
S_w21 = read_tsv("ranks_S_phase.apex_vs_correction_cellranger_w21_sct_cca_wilcox_lfc.rnk")
S_w21.gsea <- S_w21$rank
names(S_w21.gsea) <- as.character(S_w21$GeneName)
S_w21.gsea <- sort(S_w21.gsea, decreasing = T)
#for  M_w10
M_w10 = read_tsv("ranks_M.apex_vs_correction_cellranger_w10_sct_cca_wilcox_lfc.rnk")
M_w10.gsea <- M_w10$rank
names(M_w10.gsea) <- as.character(M_w10$GeneName)
M_w10.gsea <- sort(M_w10.gsea, decreasing = T)
#for  M_w21
M_w21 = read_tsv("ranks_M_phase.apex_vs_correction_cellranger_w21_sct_cca_wilcox_lfc.rnk")
M_w21.gsea <- M_w21$rank
names(M_w21.gsea) <- as.character(M_w21$GeneName)
M_w21.gsea <- sort(M_w21.gsea, decreasing = T)
#for  vRG_w10
vRG_w10 = read_tsv("ranks_vRG.apex_vs_correction_cellranger_w10_sct_cca_wilcox_lfc.rnk")
vRG_w10.gsea <- vRG_w10$rank
names(vRG_w10.gsea) <- as.character(vRG_w10$GeneName)
vRG_w10.gsea <- sort(vRG_w10.gsea, decreasing = T)
#for  tRG_w10
tRG_w10 = read_tsv("ranks_tRG.apex_vs_correction_cellranger_w10_sct_cca_wilcox_lfc.rnk")
tRG_w10.gsea <- tRG_w10$rank
names(tRG_w10.gsea) <- as.character(tRG_w10$GeneName)
tRG_w10.gsea <- sort(tRG_w10.gsea, decreasing = T)
#for  tRG_w21
tRG_w21 = read_tsv("ranks_tRG.apex_vs_correction_cellranger_w21_sct_cca_wilcox_lfc.rnk")
tRG_w21.gsea <- tRG_w21$rank
names(tRG_w21.gsea) <- as.character(tRG_w21$GeneName)
tRG_w21.gsea <- sort(tRG_w21.gsea, decreasing = T)
#for  IPC_w10
IPC_w10 = read_tsv("ranks_IPC.apex_vs_correction_cellranger_w10_sct_cca_wilcox_lfc.rnk")
IPC_w10.gsea <- IPC_w10$rank
names(IPC_w10.gsea) <- as.character(IPC_w10$GeneName)
IPC_w10.gsea <- sort(IPC_w10.gsea, decreasing = T)
#for  IPC_w21
IPC_w21 = read_tsv("ranks_IPC.apex_vs_correction_cellranger_w21_sct_cca_wilcox_lfc.rnk")
IPC_w21.gsea <- IPC_w21$rank
names(IPC_w21.gsea) <- as.character(IPC_w21$GeneName)
IPC_w21.gsea <- sort(IPC_w21.gsea, decreasing = T)
#for  oRG_w10
oRG_w10 = read_tsv("ranks_oRG.apex_vs_correction_cellranger_w10_sct_cca_wilcox_lfc.rnk")
oRG_w10.gsea <- oRG_w10$rank
names(oRG_w10.gsea) <- as.character(oRG_w10$GeneName)
oRG_w10.gsea <- sort(oRG_w10.gsea, decreasing = T)
#for  oRG_w21
oRG_w21 = read_tsv("ranks_oRG.apex_vs_correction_cellranger_w21_sct_cca_wilcox_lfc.rnk")
oRG_w21.gsea <- oRG_w21$rank
names(oRG_w21.gsea) <- as.character(oRG_w21$GeneName)
oRG_w21.gsea <- sort(oRG_w21.gsea, decreasing = T)


inputList <- list("S_w10" = S_w10.gsea,
                  "S_w21" = S_w21.gsea,
                  "M_w10" = M_w10.gsea,
                  "M_w21" = M_w21.gsea,
                  "vRG_w10" = vRG_w10.gsea,
                  "tRG_w10" = tRG_w10.gsea,
                  "tRG_w21" = tRG_w21.gsea,
                  "IPC_w10" = IPC_w10.gsea,
                  "IPC_w21" = IPC_w21.gsea,
                  "oRG_w10" = oRG_w10.gsea,
                  "oRG_w21" = oRG_w21.gsea
)


set.seed(42)
test.out<- compareCluster(geneCluster=inputList, fun="GSEA", pvalueCutoff=0.05,
                          pAdjustMethod="BH", TERM2GENE = hs_gsea_c5, 
                          minGSSize = 10, maxGSSize = 500,
                          seed = T)
saveRDS(test.out, "w10_w21_NPC_cluster_apex_vs_correction_comparecluster_10_500_8.11.22.rds")
test.out <- readRDS("w10_w21_NPC_cluster_induction_vs_chap_comparecluster_10_500.rds")
selected_GO <- c("REGULATION_OF_NEURON_PROJECTION_DEVELOPMENT", 
               "AXON_DEVELOPMENT",
               "FOREBRAIN_DEVELOPMENT", 
               "REGULATION_OF_AXONOGENESIS",
               "TELENCEPHALON_DEVELOPMENT",
               "SENSORY_SYSTEM_DEVELOPMENT",
               "CAMERA_TYPE_EYE_DEVELOPMENT",
               "REGULATION_OF_NEUROGENESIS",
               "MITOCHONDRIAL_TRANSLATION",
               "MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",
               "ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
               "STEROL_BIOSYNTHETIC_PROCESS",
               "PROTEIN_FOLDING",
               "ATP_METABOLIC_PROCESS",
               "RIBOSOME_BIOGENESIS",
               "SMALL_MOLECULE_BIOSYNTHETIC_PROCESS",
               "OXIDATIVE_PHOSPHORYLATION",
               "CYTOPLASMIC_TRANSLATION",
               "EMBRYONIC_ORGAN_MORPHOGENESIS",
               "CEREBRAL_CORTEX_DEVELOPMENT",
               "MICROTUBULE_BASED_MOVEMENT",
               "MICROTUBULE_POLYMERIZATION",
               "REGULATION_OF_GLIAL_CELL_DIFFERENTIATION",
               "GLIAL_CELL_DIFFERENTIATION",
               "REGULATION_OF_GLIOGENESIS"
)



write.csv(test.out, "test_apex_vs_corr_NPC_subtype_comparecluster.csv")
dotplot(test.out, showCategory=selected_GO, split=".sign", label_format = 60, includeAll = T) + facet_grid(.~.sign) + theme(axis.text.x=element_text(angle=90, hjust=1)) 
