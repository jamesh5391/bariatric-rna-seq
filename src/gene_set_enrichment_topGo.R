library(topGO)
library(org.Hs.eg.db)
library(biomaRt)

# Load genes and convert IDs
genes <- read.csv("results/significant_genes_T3_vs_T0.csv")$gene
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
conversion <- getBM(c("ensembl_gene_id", "entrezgene_id"), "ensembl_gene_id", genes, ensembl)
entrez_ids <- conversion$entrezgene_id[!is.na(conversion$entrezgene_id)]

# Create gene list
all_genes <- keys(org.Hs.eg.db, "ENTREZID")  
gene_list <- factor(as.integer(all_genes %in% entrez_ids))
names(gene_list) <- all_genes

# Run GO enrichment
GOdata <- new("topGOdata", ontology = "BP", allGenes = gene_list, 
              geneSel = function(x) x == 1, annot = annFUN.org, 
              mapping = "org.Hs.eg.db", ID = "entrez")
result <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
results_all <- GenTable(GOdata, fisher = result, topNodes = length(usedGO(GOdata)))
results_top10 <- head(results_all, 10)

print(results_top10)

# Create results directory and save all results
if (!dir.exists("results")) dir.create("results")
write.csv(results_all, "results/topGO_enrichment_results.csv", row.names = FALSE)
cat(paste("All", nrow(results_all), "results saved to: results/topGO_enrichment_results.csv\n"))
