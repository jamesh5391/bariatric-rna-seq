if (!require(gprofiler2, quietly = TRUE)) {
  install.packages("gprofiler2")
  library(gprofiler2)
}

# Load genes 
genes <- read.csv("results/significant_genes_T3_vs_T0.csv")$gene

# Run GO enrichment using gprofiler2
gostres <- gost(query = genes, 
                organism = "hsapiens",
                ordered_query = FALSE,
                multi_query = FALSE,
                significant = FALSE,  
                exclude_iea = FALSE,
                measure_underrepresentation = FALSE,
                evcodes = FALSE,
                user_threshold = 1.0,  
                correction_method = "g_SCS",
                domain_scope = "annotated",
                custom_bg = NULL,
                numeric_ns = "",
                sources = c("GO:BP", "GO:MF", "GO:CC"))

# Extract results
if(!is.null(gostres) && !is.null(gostres$result)) {
  results_all <- gostres$result
  
  # Convert list columns to character strings for CSV export
  if("intersection" %in% colnames(results_all)) {
    results_all$intersection <- sapply(results_all$intersection, function(x) paste(x, collapse = ";"))
  }
  if("parents" %in% colnames(results_all)) {
    results_all$parents <- sapply(results_all$parents, function(x) paste(x, collapse = ";"))
  }
  
  # Print top 10 results
  if(nrow(results_all) > 0) {
    results_top10 <- head(results_all, 10)
    print(results_top10[, c("term_id", "term_name", "p_value", "term_size", "query_size", "intersection_size")])
    
    # Create results directory and save ALL results
    if (!dir.exists("results")) dir.create("results")
    write.csv(results_all, "results/gprofiler2_enrichment_results.csv", row.names = FALSE)
    cat(paste("All", nrow(results_all), "results saved to: results/gprofiler2_enrichment_results.csv\n"))
  } else {
    cat("No enrichment results found.\n")
  }
} else {
  cat("No enrichment results returned from gprofiler2.\n")
}
