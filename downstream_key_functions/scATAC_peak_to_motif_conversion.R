library("parallel")
library("motifmatchr")
library("chromVAR")
library("Matrix")
library("SummarizedExperiment")
library("GenomicRanges")
library("BSgenome.Hsapiens.UCSC.hg38")

###function converting peak characters to a grange obj
from_char_to_GRanges <- function(input_peak_char){
    ###split
    peak_chrom <- sapply(input_peak_char, function(x) strsplit(x, split = ":", fixed = T)[[1]][1])
    peak_start <- sapply(input_peak_char, function(x) as.numeric(strsplit(strsplit(x, split = ":", fixed = T)[[1]][2], split = "-", fixed = T)[[1]][1]))
    peak_end <- sapply(input_peak_char, function(x) as.numeric(strsplit(strsplit(x, split = ":", fixed = T)[[1]][2], split = "-", fixed = T)[[1]][2]))
    
    ###reformat into granges                   
    out_granges <- GRanges(seqnames = peak_chrom, ranges = IRanges(start = peak_start, end = peak_end), names = input_peak_char)
    return(out_granges)                   
}

###write the function to obtain the intersected ATAC cells
trajectory_peak_motif_input_prep <- function(ctrl_RNA_meta_with_pseudotime, ATAC_meta, ATAC_count_mat, n_cells=50){
    
    ###intersect the ATAC meta
    ATAC_meta_with_pseudotime <- ATAC_meta %>% filter(concensus_cell_names %in% ctrl_RNA_meta_with_pseudotime$concensus_cell_names)
    ATAC_meta_with_pseudotime$pseudotime <- ctrl_RNA_meta_with_pseudotime[ATAC_meta_with_pseudotime$concensus_cell_names, "corrected_pseudotime"]    
    ###get peak counts matrix
    ATAC_with_pseudotime_count_mat <- ATAC_count_mat[,ATAC_meta_with_pseudotime$ATAC_names_with_prefix]
    colnames(ATAC_with_pseudotime_count_mat) <- ATAC_meta_with_pseudotime$concensus_cell_names
    rownames(ATAC_meta_with_pseudotime) <- ATAC_meta_with_pseudotime$concensus_cell_names
    
    ###only keep peaks detected in at least N cells
    binarized_counts_mat <- ATAC_with_pseudotime_count_mat
    binarized_counts_mat@x[binarized_counts_mat@x > 0] <- 1

    n_cell_expr <- Matrix::rowSums(binarized_counts_mat)
    ATAC_with_pseudotime_count_mat <- ATAC_with_pseudotime_count_mat[n_cell_expr >= n_cells,]
    
    ###get rowdata and rowranges
    rowdata_df <- data.frame(peak_names = rownames(ATAC_with_pseudotime_count_mat), peak_ranges=rownames(ATAC_with_pseudotime_count_mat))
    rowranges <- from_char_to_GRanges(rownames(ATAC_with_pseudotime_count_mat))
    
    ###create a SE obj
    ATAC_SE_obj <- SummarizedExperiment(assays = list(counts = ATAC_with_pseudotime_count_mat), 
                                        rowRanges = rowranges, colData = ATAC_meta_with_pseudotime)
    
    ###convert peaks to motif counts
    motifs <- getJasparMotifs()
    motif_ix <- matchMotifs(motifs, ATAC_SE_obj, 
                            genome = BSgenome.Hsapiens.UCSC.hg38, out="scores")
    
    ###convert motif counts
    motif_counts_dense_mat <- as.matrix(assays(motif_ix)[["motifCounts"]])
    #return(motif_counts_dense_mat)
    ###get single cell-motif matrix
    uniq_col_index <- unique(assays(ATAC_SE_obj)[["counts"]]@j)
    col_index <- (assays(ATAC_SE_obj)[["counts"]]@j)
    row_index <- (assays(ATAC_SE_obj)[["counts"]]@i)
    
    cell_spec_peaks_indices <- parallel::mclapply(uniq_col_index, mc.cores=16, function(x){
        cell_spec_index <- which(col_index == x)
        peak_index <- row_index[cell_spec_index]
        return(peak_index)
    })
    names(cell_spec_peaks_indices) <- as.character(uniq_col_index)
    #return(cell_spec_peaks_indices)
    cell_spec_motif_counts <- parallel::mclapply(as.character(uniq_col_index), mc.cores=16, function(x){
        peak_row_index <- cell_spec_peaks_indices[[x]]
        motif_counts <- colSums(motif_counts_dense_mat[peak_row_index,])
        return(motif_counts)
    })
    motif_cell_mat <- do.call(cbind, cell_spec_motif_counts)
    rownames(motif_cell_mat) <- colnames(motif_counts_dense_mat)
    colnames(motif_cell_mat) <- colnames(assays(ATAC_SE_obj)[["counts"]])
    
    ###get out df
    out_list <- list(motif_cell_mat=motif_cell_mat,
                     ATAC_meta_with_pseudotime=ATAC_meta_with_pseudotime)
    
    return(out_list)
}

#### write a function to modify the row names of the matrix
simplify_motif_mat <- function(input_motif_matrix) {
    
    # Extract the original row names
    input_motif_vector <- rownames(input_motif_matrix)
    
    # Simplify the gene names
    trimmed_gene_names <- sapply(input_motif_vector, function(x) {
        first_round_splitting <- strsplit(x, split = "_", fixed = TRUE)[[1]][2]
        sec_round_splitting <- strsplit(first_round_splitting, split = "(", fixed = TRUE)[[1]][1]
        return(sec_round_splitting)
    })
    
    # Identify rows containing "::" and split them into two elements
    expanded_rows <- lapply(seq_along(trimmed_gene_names), function(i) {
        name <- trimmed_gene_names[i]
        if (grepl("::", name)) {
            # Split by "::"
            split_names <- strsplit(name, "::", fixed = TRUE)[[1]]
            # Duplicate the row in the matrix for each element in the split
            duplicated_rows <- lapply(split_names, function(new_name) {
                row_data <- as.numeric(input_motif_matrix[i, , drop = FALSE])
                return(row_data)
            })
            return(do.call(rbind, duplicated_rows))
        } else {
            row_data <- as.numeric(input_motif_matrix[i, , drop = FALSE])
            return(row_data)
        }
    })
    expanded_matrix <- do.call(rbind, expanded_rows)
    
    split_rownames <- lapply(seq_along(trimmed_gene_names), function(i) {
        name <- trimmed_gene_names[i]
        if (grepl("::", name)) {
            # Split by "::"
            split_names <- strsplit(name, "::", fixed = TRUE)[[1]]
            return(split_names)
            } else {
            # If no "::", return the row as is
            split_names <- name
            return(split_names)
        }
    })
    split_rownames <- do.call(c, split_rownames)
    
    # Average rows with the same label
    unique_names <- unique(split_rownames)
    averaged_matrix_list <- lapply(unique_names, function(name) {
        rows_to_average <- expanded_matrix[split_rownames == name, , drop = FALSE]
        message(name)
        if (nrow(rows_to_average) > 1) {
            return(colMeans(rows_to_average))
        } else {
            return(as.vector(rows_to_average))
        }
    })
    
    # Convert the result back to a matrix and set row names
    averaged_matrix <- do.call(rbind, averaged_matrix_list)
    rownames(averaged_matrix) <- unique_names
    
    return(averaged_matrix)
}

