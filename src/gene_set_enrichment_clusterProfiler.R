if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

required_packages <- c("clusterProfiler", "org.Hs.eg.db", "biomaRt")
for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        BiocManager::install(pkg, ask = FALSE)
        library(pkg, character.only = TRUE)
    }
}

# Load data
sig_genes <- read.csv("results/significant_genes_T3_vs_T0.csv", stringsAsFactors = FALSE)
significant_gene_ids <- sig_genes$gene

# Convert Ensembl IDs to Entrez IDs
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_conversion <- getBM(
    attributes = c("ensembl_gene_id", "entrezgene_id"),
    filters = "ensembl_gene_id",
    values = significant_gene_ids,
    mart = ensembl
)

gene_conversion <- gene_conversion[!is.na(gene_conversion$entrezgene_id), ]
entrez_ids <- as.character(gene_conversion$entrezgene_id)

# GO Biological Process enrichment
go_bp <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 1.0,  # Set to 1.0 to get ALL results
                  readable = TRUE)

# Print results
cat("GO Biological Process enrichment results:\n")
print(head(go_bp@result, 10))

# Create results directory and save all results
if (!dir.exists("results")) dir.create("results")
write.csv(go_bp@result, "results/clusterProfiler_enrichment_results.csv", row.names = FALSE)
cat(paste("All", nrow(go_bp@result), "results saved to: results/clusterProfiler_enrichment_results.csv\n"))
