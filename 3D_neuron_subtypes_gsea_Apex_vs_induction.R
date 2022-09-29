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
#for  EN_deep_w10
EN_deep_w10 = read_tsv("ranks_EN_deep.apex_vs_induction_cellranger_w10_sct_cca_wilcox_lfc.rnk")
EN_deep_w10.gsea <- EN_deep_w10$rank
names(EN_deep_w10.gsea) <- as.character(EN_deep_w10$GeneName)
EN_deep_w10.gsea <- sort(EN_deep_w10.gsea, decreasing = T)
#for  EN_deep_w21
EN_deep_w21 = read_tsv("ranks_EN_deep.apex_vs_induction_cellranger_w21_sct_cca_wilcox_lfc.rnk")
EN_deep_w21.gsea <- EN_deep_w21$rank
names(EN_deep_w21.gsea) <- as.character(EN_deep_w21$GeneName)
EN_deep_w21.gsea <- sort(EN_deep_w21.gsea, decreasing = T)
#for  EN_upper_w10
EN_upper_w10 = read_tsv("ranks_EN_upper.apex_vs_induction_cellranger_w10_sct_cca_wilcox_lfc.rnk")
EN_upper_w10.gsea <- EN_upper_w10$rank
names(EN_upper_w10.gsea) <- as.character(EN_upper_w10$GeneName)
EN_upper_w10.gsea <- sort(EN_upper_w10.gsea, decreasing = T)
#for  EN_upper_w21
EN_upper_w21 = read_tsv("ranks_EN_upper.apex_vs_induction_cellranger_w21_sct_cca_wilcox_lfc.rnk")
EN_upper_w21.gsea <- EN_upper_w21$rank
names(EN_upper_w21.gsea) <- as.character(EN_upper_w21$GeneName)
EN_upper_w21.gsea <- sort(EN_upper_w21.gsea, decreasing = T)
#for  IN_w10
IN_w10 = read_tsv("ranks_IN.apex_vs_induction_cellranger_w10_sct_cca_wilcox_lfc.rnk")
IN_w10.gsea <- IN_w10$rank
names(IN_w10.gsea) <- as.character(IN_w10$GeneName)
IN_w10.gsea <- sort(IN_w10.gsea, decreasing = T)
#for  IN_w21
IN_w21 = read_tsv("ranks_IN.apex_vs_induction_cellranger_w21_sct_cca_wilcox_lfc.rnk")
IN_w21.gsea <- IN_w21$rank
names(IN_w21.gsea) <- as.character(IN_w21$GeneName)
IN_w21.gsea <- sort(IN_w21.gsea, decreasing = T)

inputList <- list(
  "EN_deep_w10" = EN_deep_w10.gsea,
  "EN_deep_w21" = EN_deep_w21.gsea,
  "EN_upper_w10" = EN_upper_w10.gsea,
  "EN_upper_w21" = EN_upper_w21.gsea,
  "IN_w10" = IN_w10.gsea,
  "IN_w21" = IN_w21.gsea
  
)


set.seed(42)
test.out<- compareCluster(geneCluster=inputList, fun="GSEA", pvalueCutoff=0.05,
                          pAdjustMethod="BH", TERM2GENE = hs_gsea_c5, 
                          minGSSize = 10, maxGSSize = 500,
                          seed = T)
saveRDS(test.out, "w10_w21_EN_d_u_IN_neuron_cluster_apex_vs_induction_comparecluster_10_500_8.11.22.rds")

selected_GO <- c("SYNAPTIC_SIGNALING", 
               "RIBOSOME_BIOGENESIS",
               "ACTION_POTENTIAL",
               "SYNAPTIC_TRANSMISSION_GLUTAMATERGIC",
               "REGULATION_OF_TRANS_SYNAPTIC_SIGNALING",
               "CYTOPLASMIC_TRANSLATION",
               "NEURONAL_ACTION_POTENTIAL",
               "SYNAPTIC_TRANSMISSION_GLUTAMATERGIC",
               "ATP_METABOLIC_PROCESS",
               "MITOCHONDRIAL_ELECTRON_TRANSPORT_CYTOCHROME_C_TO_OXYGEN",
               "ELECTRON_TRANSPORT_CHAIN",
               "MITOCHONDRIAL_TRANSLATION",
               "MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_ASSEMBLY",
               "ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
               "REGULATION_OF_SYNAPTIC_PLASTICITY"
               
               
)

dotplot(test.out, showCategory=selected_GO, split=".sign", label_format = 60) + facet_grid(.~.sign) + theme(axis.text.x=element_text(angle=90, hjust=1)) 



