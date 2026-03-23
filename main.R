library('tidyverse')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('testthat')
library('fgsea')

#' Function to generate a SummarizedExperiment object with counts and coldata
#' to use in DESeq2
#'
#' @param csv_path (str): path to the file verse_counts.tsv
#' @param metafile (str): path to the metadata sample_metadata.csv
#' @param selected_times (list): list of sample timepoints to use
#' 
#'   
#' @return SummarizedExperiment object with subsetted counts matrix
#'   and sample data. Ensure that the timepoints column used as input 
#'   to the model design has 'vP0' set as the reference factor level. Your 
#'   colData dataframe should have columns named samplename and timepoint.
#' @export
#'
#' @examples se <- make_se('verse_counts.tsv', 'sample_metadata.csv', c('vP0', 'vAd'))
make_se <- function(counts_csv, metafile_csv, selected_times) {
  # Read in counts matrix
  counts <- read.table(counts_csv, header = TRUE, sep = "\t", row.names = 1)
  
  # Read in metadata
  metadata <- read.csv(metafile_csv, header = TRUE)
  
  # Subset metadata to selected timepoints
  metadata <- metadata %>%
    filter(timepoint %in% selected_times) %>%
    dplyr::select(samplename, timepoint)
  
  # Subset counts matrix to only include the selected samples
  counts <- counts[, metadata$samplename]
  
  # Set vP0 as the reference factor level
  metadata$timepoint <- factor(metadata$timepoint, levels = c("vP0", "vAd"))
  
  # Convert counts to matrix
  counts_matrix <- as.matrix(counts)
  
  # Build SummarizedExperiment object
  se <- SummarizedExperiment(
    assays = list(counts = counts_matrix),
    colData = metadata
  )
  
  return(se)
}

#' Function that runs DESeq2 and returns a named list containing the DESeq2
#' results as a dataframe and the dds object returned by DESeq2
#'
#' @param se (obj): SummarizedExperiment object containing counts matrix and
#' coldata
#' @param design: the design formula to be used in DESeq2
#'
#' @return list with DESeqDataSet object after running DESeq2 and results from
#'   DESeq2 as a dataframe
#' @export
#'
#' @examples results <- return_deseq_res(se, ~ timepoint)
return_deseq_res <- function(se, design) {
  # Create DESeqDataSet from SummarizedExperiment
  dds <- DESeqDataSet(se, design = design)
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Extract results as dataframe
  res <- results(dds)
  res_df <- as.data.frame(res)
  
  # Return both as a named list
  return(list(dds = dds, results = res_df))
}

#' Function that takes the DESeq2 results dataframe, converts it to a tibble and
#' adds a column to denote plotting status in volcano plot. Column should denote
#' whether gene is either 1. Significant at padj < .10 and has a positive log
#' fold change, 2. Significant at padj < .10 and has a negative log fold change,
#' 3. Not significant at padj < .10. Have the values for these labels be UP,
#' DOWN, NS, respectively. The column should be named `volc_plot_status`. Ensure
#' that the column name for your rownames is called "genes". 
#'
#' @param deseq2_res (df): results from DESeq2 
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return Tibble with all columns from DESeq2 results and one additional column
#'   labeling genes by significant and up-regulated, significant and
#'   downregulated, and not significant at padj < .10.
#'   
#' @export
#'
#' @examples labeled_results <- label_res(res, .10)
label_res <- function(deseq2_res, padj_threshold) {
  labeled <- deseq2_res %>%
    rownames_to_column("genes") %>%
    as_tibble() %>%
    mutate(volc_plot_status = case_when(
      padj < padj_threshold & log2FoldChange > 0 ~ "UP",
      padj < padj_threshold & log2FoldChange < 0 ~ "DOWN",
      TRUE ~ "NS"
    ))
  
  return(labeled)
}

#' Function to plot the unadjusted p-values as a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#'
#' @return ggplot: a histogram of the raw p-values from the DESeq2 results
#' @export
#'
#' @examples pval_plot <- plot_pvals(labeled_results)
plot_pvals <- function(labeled_results) {
  ggplot(labeled_results, aes(x = pvalue)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "white") +
    theme_bw() +
    labs(
      title = "Distribution of Raw P-values",
      x = "P-value",
      y = "Count"
    )
}

#' Function to plot the log2foldchange from DESeq2 results in a histogram
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one additional
#' column denoting status in volcano plot
#' @param padj_threshold (float): threshold for considering significance (padj)
#'
#' @return ggplot: a histogram of log2FC values from genes significant at padj 
#' threshold of 0.1
#' @export
#'
#' @examples log2fc_plot <- plot_log2fc(labeled_results, .10)
plot_log2fc <- function(labeled_results, padj_threshold) {
  sig_genes <- labeled_results %>%
    filter(padj < padj_threshold)
  
  ggplot(sig_genes, aes(x = log2FoldChange)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "white") +
    theme_bw() +
    labs(
      title = "Distribution of Log2 Fold Changes for DE Genes",
      x = "Log2 Fold Change",
      y = "Count"
    )
}

#' Function to make scatter plot of normalized counts for top ten genes ranked
#' by ascending padj
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param dds_obj (obj): The object returned by running DESeq (dds) containing
#' the updated DESeqDataSet object with test results
#' @param num_genes (int): Number of genes to plot
#'
#' @return ggplot: a scatter plot with the normalized counts for each sample for
#' each of the top ten genes ranked by ascending padj
#' @export
#'
#' @examples norm_counts_plot <- scatter_norm_counts(labeled_results, dds, 10)
scatter_norm_counts <- function(labeled_results, dds_obj, num_genes){
  # Get top genes ranked by ascending padj
  top_genes <- labeled_results %>%
    filter(!is.na(padj)) %>%
    arrange(padj) %>%
    head(num_genes) %>%
    pull(genes)
  
  # Extract normalized counts
  norm_counts <- counts(dds_obj, normalized = TRUE)
  
  # Subset to top genes and convert to long format
  norm_counts_df <- norm_counts[top_genes, ] %>%
    as.data.frame() %>%
    rownames_to_column("genes") %>%
    pivot_longer(-genes, names_to = "samplename", values_to = "counts")
  
  ggplot(norm_counts_df, aes(x = genes, y = log10(counts + 1), color = samplename)) +
    geom_jitter(width = 0.2, size = 3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Normalized Counts for Top 10 DE Genes",
      x = "Gene",
      y = "Log10 Normalized Counts",
      color = "Sample"
    )
}

#' Function to generate volcano plot from DESeq2 results
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#'
#' @return ggplot: a scatterplot (volcano plot) that displays log2foldchange vs
#'   -log10(padj) and labeled by status
#' @export
#'
#' @examples volcano_plot <- plot_volcano(labeled_results)
#' 
plot_volcano <- function(labeled_results) {
  ggplot(labeled_results, aes(x = log2FoldChange, y = -log10(padj), color = volc_plot_status)) +
    geom_point(alpha = 0.5, size = 1) +
    scale_color_manual(values = c("UP" = "red", "DOWN" = "blue", "NS" = "grey")) +
    theme_bw() +
    labs(
      title = "Volcano Plot of DESeq2 Results",
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-value",
      color = "Status"
    )
}

#' Function to generate a named vector ranked by log2FC descending
#'
#' @param labeled_results (tibble): Tibble with DESeq2 results and one
#'   additional column denoting status in volcano plot
#' @param id2gene_path (str): Path to the file containing the mapping of
#' ensembl IDs to MGI symbols
#'
#' @return Named vector with gene symbols as names, and log2FoldChange as values
#' ranked in descending order
#' @export
#'
#' @examples rnk_list <- make_ranked_log2fc(labeled_results, 'data/id2gene.txt')

make_ranked_log2fc <- function(labeled_results, id2gene_path) {
  # Read in the id2gene mapping file without header
  id2gene <- read.table(id2gene_path, header = FALSE, sep = "\t",
                        col.names = c("ensembl_gene_id", "mgi_symbol"))
  
  # Join the labeled results with the id2gene mapping
  labeled_results <- labeled_results %>%
    left_join(id2gene, by = c("genes" = "ensembl_gene_id"))
  
  # Remove duplicates and NAs
  ranked <- labeled_results %>%
    filter(!is.na(log2FoldChange)) %>%
    filter(!is.na(mgi_symbol)) %>%
    filter(mgi_symbol != "") %>%
    distinct(mgi_symbol, .keep_all = TRUE) %>%
    arrange(desc(log2FoldChange))
  
  # Create named vector
  rnk_list <- ranked$log2FoldChange
  names(rnk_list) <- ranked$mgi_symbol
  
  return(rnk_list)
}

#' Function to run fgsea with arguments for min and max gene set size
#'
#' @param gmt_file_path (str): Path to the gene sets of interest in GMT format
#' @param rnk_list (named vector): Named vector generated previously with gene 
#' symbols and log2Fold Change values in descending order
#' @param min_size (int): Minimum number of genes in gene sets to be allowed
#' @param max_size (int): Maximum number of genes in gene sets to be allowed
#'
#' @return Tibble of results from running fgsea
#' @export
#'
#' @examples fgsea_results <- run_fgsea('data/m2.cp.v2023.1.Mm.symbols.gmt', rnk_list, 15, 500)
run_fgsea <- function(gmt_file_path, rnk_list, min_size, max_size) {
  # Read in the GMT file
  pathways <- gmtPathways(gmt_file_path)
  
  # Run fgsea
  fgsea_res <- fgsea(pathways = pathways,
                     stats = rnk_list,
                     minSize = min_size,
                     maxSize = max_size)
  
  # Convert to tibble
  fgsea_res <- as_tibble(fgsea_res)
  
  return(fgsea_res)
}

#' Function to plot top ten positive NES and top ten negative NES pathways
#' in a barchart
#'
#' @param fgsea_results (tibble): the fgsea results in tibble format returned by
#'   the previous function
#' @param num_paths (int): the number of pathways for each direction (top or
#'   down) to include in the plot. Set this at 10.
#'
#' @return ggplot with a barchart showing the top twenty pathways ranked by positive
#' and negative NES
#' @export
#'
#' @examples fgsea_plot <- top_pathways(fgsea_results, 10)
top_pathways <- function(fgsea_results, num_paths){
  # Get top positive and negative NES pathways
  top_pos <- fgsea_results %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    head(num_paths)
  
  top_neg <- fgsea_results %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    head(num_paths)
  
  # Combine and add direction column
  top_combined <- bind_rows(top_pos, top_neg) %>%
    mutate(direction = ifelse(NES > 0, "Positive", "Negative")) %>%
    arrange(NES)
  
  # Plot
  ggplot(top_combined, aes(x = reorder(pathway, NES), y = NES, fill = direction)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Positive" = "red", "Negative" = "blue")) +
    coord_flip() +
    theme_bw() +
    theme(axis.text.y = element_text(size = 6)) +
    labs(
      title = "Top 10 Positive and Negative NES Pathways",
      x = "Pathway",
      y = "Normalized Enrichment Score (NES)",
      fill = "Direction"
    )
}

