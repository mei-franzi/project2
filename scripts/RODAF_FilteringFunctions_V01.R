# MinCount Filter
min_count_filter <- function(dds, min_count) {
  counts_matrix <- counts(dds, normalized = FALSE)  # Use raw counts
  keep <- rowSums(counts_matrix >= min_count) > 0
  dds_filtered <- dds[keep, ]
  removed <- nrow(dds) - nrow(dds_filtered)
  print(paste("MinCount filter removed:", removed, "genes. Remaining:", nrow(dds_filtered)))
  return(dds_filtered)
}


relevance_filter <- function(dds, design, min_count) {
  # Ensure size factors are calculated
  if (is.null(sizeFactors(dds))) {
    dds <- estimateSizeFactors(dds)
  }
  
  # Get normalized counts (matches original script)
  counts_matrix <- counts(dds, normalized = TRUE)
  cpm_matrix <- cpm(counts_matrix)
  
  # Get sample group sizes
  group_counts <- table(design$treatment)
  
  # Apply relevance filter
  pass_filter <- sapply(rownames(counts_matrix), function(gene) {
    CountsPass <- sapply(names(group_counts), function(group) {
      sample_cols <- which(design$treatment == group)
      sum(cpm_matrix[gene, sample_cols] >= min_count) >= 0.75 * group_counts[group]
    })
    
    return(any(CountsPass))  # Keep gene if at least one group passes
  })
  
  dds_filtered <- dds[pass_filter, ]
  removed <- nrow(dds) - nrow(dds_filtered)
  print(paste("Relevance filter removed:", removed, "genes. Remaining:", nrow(dds_filtered)))
  
  return(dds_filtered)
}

spurious_spikes_filter <- function(DECounts, design) {
  # Extract normalized count matrix
  counts_matrix <- counts(DECounts, normalized = TRUE)
  
  # Get group names
  treatment_groups <- dimnames(table(design$treatment))[[1]]
  
  # Initialize filter matrix (only one column now: spike)
  Filter <- matrix(data = NA, ncol = 1, nrow = nrow(counts_matrix))
  rownames(Filter) <- rownames(counts_matrix)
  colnames(Filter) <- c("spike")
  
  # Apply spike check
  spike_check <- apply(counts_matrix, 1, function(gene_counts) {
    spikePass <- sapply(treatment_groups, function(group) {
      sampleCols <- grep(group, design$treatment)
      if (max(gene_counts[sampleCols]) == 0) {
        return(FALSE)
      } else {
        return((max(gene_counts[sampleCols]) / sum(gene_counts[sampleCols])) >= 1.4 * (length(sampleCols))^(-0.66))
      }
    })
    return(sum(spikePass) <= 1)  # Spike filter fails if >1 group has a spike
  })
  Filter[, 1] <- as.numeric(spike_check)
  
  # Keep genes that pass the spike check
  keep_genes <- Filter[, "spike"] == 1
  DECounts_filtered <- DECounts[keep_genes, ]
  removed <- nrow(DECounts) - nrow(DECounts_filtered)
  
  print(paste("Spurious spikes filter removed:", removed, "genes. Remaining:", nrow(DECounts_filtered)))
  return(DECounts_filtered)
}

