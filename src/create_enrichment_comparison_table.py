"""
Create a comprehensive comparison table of gene set enrichment results
from multiple methods (clusterProfiler, gprofiler2, topGO).

This script combines results from all three methods into a single table
where each row represents a unique gene set/term and columns show results
from each method along with summary statistics.
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path

def parse_p_value(p_val_str):
    """
    Parse p-value strings that might contain '<' or other symbols.
    Returns float value or NaN if cannot parse.
    """
    if pd.isna(p_val_str):
        return np.nan
    
    p_str = str(p_val_str).strip()
    
    # Handle cases like "< 1e-30"
    if p_str.startswith('<'):
        p_str = p_str[1:].strip()
    
    try:
        return float(p_str)
    except ValueError:
        return np.nan

def load_clusterprofiler_results(file_path):
    """Load and process clusterProfiler results."""
    df = pd.read_csv(file_path)
    
    # Standardize column names and select relevant columns
    result = pd.DataFrame({
        'term_id': df['ID'],
        'term_name': df['Description'],
        'clusterprofiler_pvalue': df['pvalue'],
        'clusterprofiler_qvalue': df['qvalue'],
        'clusterprofiler_p_adjust': df['p.adjust'],
        'clusterprofiler_fold_enrichment': df['FoldEnrichment'],
        'clusterprofiler_gene_count': df['Count'],
        'clusterprofiler_gene_ratio': df['GeneRatio'],
        'clusterprofiler_genes': df['geneID']
    })
    
    # Add significance flag (typically q-value < 0.05)
    result['clusterprofiler_significant'] = result['clusterprofiler_qvalue'] < 0.05
    
    return result

def load_gprofiler2_results(file_path):
    """Load and process gprofiler2 results."""
    df = pd.read_csv(file_path)
    
    # Filter for GO terms only (you can modify this if needed)
    df = df[df['source'] == 'GO:BP'].copy()
    
    result = pd.DataFrame({
        'term_id': df['term_id'],
        'term_name': df['term_name'],
        'gprofiler2_pvalue': df['p_value'],
        'gprofiler2_significant': df['significant'],
        'gprofiler2_term_size': df['term_size'],
        'gprofiler2_intersection_size': df['intersection_size'],
        'gprofiler2_precision': df['precision'],
        'gprofiler2_recall': df['recall']
    })
    
    return result

def load_topgo_results(file_path):
    """Load and process topGO results."""
    df = pd.read_csv(file_path)
    
    result = pd.DataFrame({
        'term_id': df['GO.ID'],
        'term_name': df['Term'],
        'topgo_pvalue_raw': df['fisher'],
        'topgo_annotated': df['Annotated'],
        'topgo_significant': df['Significant'],
        'topgo_expected': df['Expected']
    })
    
    # Parse p-values (handle "< 1e-30" format)
    result['topgo_pvalue'] = result['topgo_pvalue_raw'].apply(parse_p_value)
    
    # Add significance flag (typically p-value < 0.05)
    result['topgo_significant_flag'] = result['topgo_pvalue'] < 0.05
    
    return result

def create_comparison_table(clusterprofiler_file, gprofiler2_file, topgo_file, significance_threshold=0.05):
    """
    Create comprehensive comparison table combining all three methods.
    """
    print("Loading results from each method...")
    
    # Load results from each method
    cp_results = load_clusterprofiler_results(clusterprofiler_file)
    gp_results = load_gprofiler2_results(gprofiler2_file)
    tg_results = load_topgo_results(topgo_file)
    
    print(f"Loaded {len(cp_results)} clusterProfiler results")
    print(f"Loaded {len(gp_results)} gprofiler2 results")
    print(f"Loaded {len(tg_results)} topGO results")
    
    # Get all unique terms across all methods
    all_terms = set()
    all_terms.update(cp_results['term_id'].dropna())
    all_terms.update(gp_results['term_id'].dropna())
    all_terms.update(tg_results['term_id'].dropna())
    
    print(f"Found {len(all_terms)} unique terms across all methods")
    
    # Create master dataframe
    master_df = pd.DataFrame({'term_id': list(all_terms)})
    
    # Merge results from each method
    master_df = master_df.merge(cp_results, on='term_id', how='left', suffixes=('', '_cp'))
    master_df = master_df.merge(gp_results, on='term_id', how='left', suffixes=('', '_gp'))
    master_df = master_df.merge(tg_results, on='term_id', how='left', suffixes=('', '_tg'))
    
    # Resolve term names (use the first non-null name found)
    def get_best_term_name(row):
        for col in ['term_name', 'term_name_gp', 'term_name_tg']:
            if col in row and pd.notna(row[col]):
                return row[col]
        return 'Unknown'
    
    master_df['term_name_final'] = master_df.apply(get_best_term_name, axis=1)
    
    # Calculate summary statistics
    # Methods that included this term in analysis
    master_df['methods_total'] = (
        (~master_df['clusterprofiler_pvalue'].isna()).astype(int) +
        (~master_df['gprofiler2_pvalue'].isna()).astype(int) +
        (~master_df['topgo_pvalue'].isna()).astype(int)
    )
    
    # Methods that found this term significantly enriched
    cp_sig = master_df['clusterprofiler_significant'].fillna(False)
    gp_sig = master_df['gprofiler2_significant'].fillna(False)
    tg_sig = master_df['topgo_significant_flag'].fillna(False)
    
    master_df['methods_significant'] = (
        cp_sig.astype(int) +
        gp_sig.astype(int) +
        tg_sig.astype(int)
    )
    
    # Reorganize columns for better readability
    final_columns = [
        'term_id', 'term_name_final',
        'methods_significant', 'methods_total',
        'clusterprofiler_pvalue', 'clusterprofiler_qvalue', 'clusterprofiler_p_adjust',
        'clusterprofiler_fold_enrichment', 'clusterprofiler_gene_count', 'clusterprofiler_significant',
        'gprofiler2_pvalue', 'gprofiler2_significant', 'gprofiler2_term_size',
        'gprofiler2_intersection_size', 'gprofiler2_precision', 'gprofiler2_recall',
        'topgo_pvalue', 'topgo_pvalue_raw', 'topgo_significant_flag',
        'topgo_annotated', 'topgo_significant', 'topgo_expected',
        'clusterprofiler_genes'
    ]
    
    # Keep only columns that exist
    available_columns = [col for col in final_columns if col in master_df.columns]
    result_df = master_df[available_columns].copy()
    
    # Rename columns for clarity
    result_df = result_df.rename(columns={
        'term_name_final': 'term_name',
        'clusterprofiler_genes': 'clusterprofiler_gene_ids'
    })
    
    # Sort by number of methods finding significance (descending), then by best p-value
    result_df['min_pvalue'] = result_df[['clusterprofiler_pvalue', 'gprofiler2_pvalue', 'topgo_pvalue']].min(axis=1)
    result_df = result_df.sort_values(['methods_significant', 'methods_total', 'min_pvalue'], 
                                     ascending=[False, False, True])
    
    # Drop the helper column
    result_df = result_df.drop('min_pvalue', axis=1)
    
    return result_df

def main():
    """Main function to create the enrichment comparison table."""
    
    # File paths
    base_dir = Path(__file__).parent.parent
    results_dir = base_dir / "results"
    
    clusterprofiler_file = results_dir / "clusterProfiler_enrichment_results.csv"
    gprofiler2_file = results_dir / "gprofiler2_enrichment_results.csv"
    topgo_file = results_dir / "topGO_enrichment_results.csv"
    
    # Check if all files exist
    for file_path in [clusterprofiler_file, gprofiler2_file, topgo_file]:
        if not file_path.exists():
            print(f"Error: File not found - {file_path}")
            return
    
    print("Creating comprehensive enrichment comparison table...")
    print("=" * 60)
    
    # Create comparison table
    comparison_df = create_comparison_table(
        clusterprofiler_file, 
        gprofiler2_file, 
        topgo_file
    )
    
    # Save results
    output_file = results_dir / "enrichment_methods_comparison_table.csv"
    comparison_df.to_csv(output_file, index=False)
    
    print(f"\nComparison table saved to: {output_file}")
    print(f"Total terms in comparison table: {len(comparison_df)}")
    
    # Print summary statistics
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    
    print(f"Terms found by all 3 methods: {len(comparison_df[comparison_df['methods_total'] == 3])}")
    print(f"Terms found by 2 methods: {len(comparison_df[comparison_df['methods_total'] == 2])}")
    print(f"Terms found by 1 method: {len(comparison_df[comparison_df['methods_total'] == 1])}")
    
    print(f"\nTerms significant in all 3 methods: {len(comparison_df[comparison_df['methods_significant'] == 3])}")
    print(f"Terms significant in 2 methods: {len(comparison_df[comparison_df['methods_significant'] == 2])}")
    print(f"Terms significant in 1 method: {len(comparison_df[comparison_df['methods_significant'] == 1])}")
    print(f"Terms not significant in any method: {len(comparison_df[comparison_df['methods_significant'] == 0])}")
    
    # Show top consensus results
    print(f"\n" + "=" * 60)
    print("TOP 10 CONSENSUS RESULTS (significant in multiple methods)")
    print("=" * 60)
    
    consensus_results = comparison_df[comparison_df['methods_significant'] >= 2].head(10)
    if len(consensus_results) > 0:
        display_cols = ['term_id', 'term_name', 'methods_significant', 'methods_total',
                       'clusterprofiler_pvalue', 'gprofiler2_pvalue', 'topgo_pvalue']
        available_display_cols = [col for col in display_cols if col in consensus_results.columns]
        print(consensus_results[available_display_cols].to_string(index=False))
    else:
        print("No terms found significant in multiple methods")
    
    # Show method-specific results
    print(f"\n" + "=" * 60)
    print("METHOD-SPECIFIC SIGNIFICANT RESULTS (top 5 each)")
    print("=" * 60)
    
    cp_only = comparison_df[(comparison_df['clusterprofiler_significant'] == True) & 
                           (comparison_df['methods_significant'] == 1)]
    if len(cp_only) > 0:
        print(f"\nclusterProfiler-specific (top 5):")
        print(cp_only[['term_id', 'term_name', 'clusterprofiler_pvalue']].head().to_string(index=False))
    
    gp_only = comparison_df[(comparison_df['gprofiler2_significant'] == True) & 
                           (comparison_df['methods_significant'] == 1)]
    if len(gp_only) > 0:
        print(f"\ngprofiler2-specific (top 5):")
        print(gp_only[['term_id', 'term_name', 'gprofiler2_pvalue']].head().to_string(index=False))
    
    tg_only = comparison_df[(comparison_df['topgo_significant_flag'] == True) & 
                           (comparison_df['methods_significant'] == 1)]
    if len(tg_only) > 0:
        print(f"\ntopGO-specific (top 5):")
        print(tg_only[['term_id', 'term_name', 'topgo_pvalue']].head().to_string(index=False))

    print(f"\n" + "=" * 60)
    print("Analysis complete! Check the output file for the full comparison table.")
    print("=" * 60)

if __name__ == "__main__":
    main()
