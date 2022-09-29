library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(readr)
library(tidyverse)
library(msigdbr)
# choose your msigdb collection of interest
hs_gsea_c5 <- msigdbr(species = "Homo sapiens", 
                      category = "C5",
                      subcategory = "BP"
) %>% 
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 
hs_gsea_c5<- hs_gsea_c5 %>% mutate(gs_name = str_replace_all(gs_name, 'GOBP_', ''))

#loading rnk files for  P3
p3 = read_tsv("ranks_p3_induction_vs_chap_t.rnk")
p3.gsea <- p3$rank
names(p3.gsea) <- as.character(p3$GeneName)
p3.gsea <- sort(p3.gsea, decreasing = T)
#loading rnk files for  P4
p4 = read_tsv("ranks_p4_induction_vs_chap_t.rnk")
p4.gsea <- p4$rank
names(p4.gsea) <- as.character(p4$GeneName)
p4.gsea <- sort(p4.gsea, decreasing = T)
# loading rnk files for  P5
p5 = read_tsv("ranks_p5_induction_vs_chap_t.rnk")
p5.gsea <- p5$rank
names(p5.gsea) <- as.character(p5$GeneName)
p5.gsea <- sort(p5.gsea, decreasing = T)


inputList <- list("p3" = p3.gsea,
                  "p4" = p4.gsea,
                  "p5" = p5.gsea
)

# set.seed to make code reproducible between each run. 
set.seed(42)
test.out<- compareCluster(geneCluster=inputList, fun="GSEA", pvalueCutoff=0.05,
                          pAdjustMethod="BH", TERM2GENE = hs_gsea_c5, 
                          minGSSize = 10, maxGSSize = 500,
                          seed = T)
saveRDS(test.out, "2D_NPC_induction_vs_chap_comparecluster_Tx_10_500.rds")
test.out <- readRDS("2D_NPC_induction_vs_chap_comparecluster_Tx_10_500.rds")

dotplot(test.out, showCategory=5, split=".sign", label_format = 60) + facet_grid(.~.sign) + theme(axis.text.x=element_text(angle=45, hjust=1)) 






