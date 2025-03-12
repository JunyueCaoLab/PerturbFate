library("Matrix")
library("dplyr")
library("parallel")

GRN_validation.step1 <- function(perturbations_to_try, dual_omics_meta_df, tri_omics_meta_df, steady_norm_whole_tx_expr_mat, norm_nascent_mat, combined_ATAC_mat, motif_gene_df, peak_to_motif_mat, peak_to_gene_df){
    
    ###extract cells with spec perturbations and merge into pseudobulk
    ###steady state expr
    binarized_steady_expr <- steady_norm_whole_tx_expr_mat
    binarized_steady_expr@x[binarized_steady_expr@x > 0] <- 1
    steady_state_genes_passed <- rownames(binarized_steady_expr)[Matrix::rowSums(binarized_steady_expr)>100]
    
    pseudobulk_steady_state_expr <- parallel::mclapply(perturbations_to_try, mc.cores=16, function(x){
        cells_to_merge <- (dual_omics_meta_df %>% filter(sgRNA_identity == x))$concensus_cell_names
        norm_perturbation_wise_gex <- as.numeric(Matrix::rowMeans(steady_norm_whole_tx_expr_mat[,cells_to_merge]))
        return(norm_perturbation_wise_gex)
    })
    pseudobulk_steady_state_expr_mat <- do.call(cbind, pseudobulk_steady_state_expr)
    colnames(pseudobulk_steady_state_expr_mat) <- perturbations_to_try
    rownames(pseudobulk_steady_state_expr_mat) <- rownames(steady_norm_whole_tx_expr_mat)
    pseudobulk_steady_state_expr_mat <- pseudobulk_steady_state_expr_mat[steady_state_genes_passed,]
    message("Pseudobulk aggregation on steady state RNA expression matrix is done.")
    
    ###extract cells with spec perturbations and merge the ATAC data
    binarized_ATAC_mat <- combined_ATAC_mat
    binarized_ATAC_mat@x[binarized_ATAC_mat@x > 0] <- 1
    ATAC_peaks_passed <- rownames(binarized_ATAC_mat)[Matrix::rowSums(binarized_ATAC_mat)>100]
    
    pseudobulk_ATAC_expr_mat <- parallel::mclapply(perturbations_to_try, mc.cores=16, function(x){
        SciCAR_cells_to_merge <- paste0("SciCAR.", (dual_omics_meta_df %>% filter(sgRNA_identity == x))$ATAC_cell_names)
        TriSci_cells_to_merge <- paste0("TriSci.", (tri_omics_meta_df %>% filter(sgRNA_identity == x))$ATAC_cell_names)
        combined_cells <- intersect(c(SciCAR_cells_to_merge, TriSci_cells_to_merge), colnames(combined_ATAC_mat))
        
        ###for atac, merge cells and then perform normalization
        raw_count_col <- Matrix::rowSums(combined_ATAC_mat[,combined_cells])
        norm_count_col <- log((raw_count_col*1e6/sum(raw_count_col))+1)
        return(norm_count_col)
    })
    pseudobulk_ATAC_expr_mat <- do.call(cbind, pseudobulk_ATAC_expr_mat)
    colnames(pseudobulk_ATAC_expr_mat) <- perturbations_to_try
    rownames(pseudobulk_ATAC_expr_mat) <- rownames(combined_ATAC_mat)
    pseudobulk_ATAC_expr_mat <- pseudobulk_ATAC_expr_mat[ATAC_peaks_passed,]
    message("Pseudobulk aggregation on ATAC matrix is done.")
    
    ###extract cells with spec perturbations and merge nascent expr data
    binarized_nascent_expr <- norm_nascent_mat
    binarized_nascent_expr@x[binarized_nascent_expr@x > 0] <- 1
    nascent_genes_passed <- rownames(binarized_nascent_expr)[Matrix::rowSums(binarized_nascent_expr)>100]
    
    pseudobulk_nascent_expr <- parallel::mclapply(perturbations_to_try, mc.cores=16, function(x){
        cells_to_merge <- (tri_omics_meta_df %>% filter(sgRNA_identity == x))$concensus_cell_names
        norm_perturbation_wise_nascent_gex <- as.numeric(Matrix::rowMeans(norm_nascent_mat[,cells_to_merge]))
        return(norm_perturbation_wise_nascent_gex)
    })
    pseudobulk_nascent_expr_mat <- do.call(cbind, pseudobulk_nascent_expr)
    colnames(pseudobulk_nascent_expr_mat) <- perturbations_to_try
    rownames(pseudobulk_nascent_expr_mat) <- rownames(norm_nascent_mat)
    pseudobulk_nascent_expr_mat <- pseudobulk_nascent_expr_mat[nascent_genes_passed,]
    message("Pseudobulk aggregation on nascent RNA expression matrix is done.")
    
    return(list(pseudobulk_steady_state_expr_mat=pseudobulk_steady_state_expr_mat,
                pseudobulk_ATAC_expr_mat=pseudobulk_ATAC_expr_mat,
                pseudobulk_nascent_expr_mat=pseudobulk_nascent_expr_mat))
}


GRN_validation.step2 <- function(step1_output_list, perturbations_to_try, dual_omics_meta_df, tri_omics_meta_df, steady_norm_whole_tx_expr_mat, norm_nascent_mat, combined_ATAC_mat, motif_gene_df, peak_to_motif_mat, peak_to_gene_df){
    
    pseudobulk_steady_state_expr_mat <- step1_output_list$pseudobulk_steady_state_expr_mat
    pseudobulk_ATAC_expr_mat <- step1_output_list$pseudobulk_ATAC_expr_mat
    pseudobulk_nascent_expr_mat <- step1_output_list$pseudobulk_nascent_expr_mat
    
    ###get unique motifs
    unique_motifs <- intersect(rownames(pseudobulk_steady_state_expr_mat), motif_gene_df$motif)
    
    ###for each motif, find target genes
    TF_regulon_score_list <- parallel::mclapply(unique_motifs, mc.cores=16, function(x){
        
        message(paste0("Start processing regulon: ", x))
        
        ###get the steady state expr of TF
        TF_steady_state_expr <- as.numeric(pseudobulk_steady_state_expr_mat[x,perturbations_to_try])
        message("TF expression retrieved.")
        
        ###get the gene targets
        target_genes <- intersect((motif_gene_df %>% filter(motif == x))$gene, rownames(pseudobulk_nascent_expr_mat))
        target_genes_nascent_expr_mat <- (pseudobulk_nascent_expr_mat[target_genes, perturbations_to_try])
        message("Target gene nascent expression retrieved.")
        
        if(nrow(target_genes_nascent_expr_mat) > 0){
        
            ###get the peak containing motifs
            target_gene_wise_ATAC_expr <- lapply(target_genes, function(y){
                ###get peaks that contain the motif
                peaks_associated <- intersect(peak_to_gene_df[peak_to_gene_df$genes == y, "peaks"], rownames(pseudobulk_ATAC_expr_mat))

                if(length(peaks_associated) > 0){
                    peaks_with_motifs <- peaks_associated[as.logical(peak_to_motif_mat[peaks_associated, x])]

                    if(length(peaks_with_motifs) > 0){
                        if(length(peaks_with_motifs) > 1){
                            avg_ATAC_expr <- as.numeric(colMeans(pseudobulk_ATAC_expr_mat[peaks_with_motifs, perturbations_to_try]))
                            return(avg_ATAC_expr)
                        }else{
                            avg_ATAC_expr <- as.numeric(pseudobulk_ATAC_expr_mat[peaks_with_motifs, perturbations_to_try])
                            return(avg_ATAC_expr)
                        }
                    }else{
                        return(NA)
                    }
                }else{
                    return(NA)
                }
            })

            if(length(target_gene_wise_ATAC_expr) > 1){
                target_gene_wise_ATAC_expr_mat <- do.call(rbind, target_gene_wise_ATAC_expr)
                rownames(target_gene_wise_ATAC_expr_mat) <- target_genes
                colnames(target_gene_wise_ATAC_expr_mat) <- perturbations_to_try
            }else{
                target_gene_wise_ATAC_expr_mat <- matrix(target_gene_wise_ATAC_expr[[1]], nrow=1)
                rownames(target_gene_wise_ATAC_expr_mat) <- target_genes
                colnames(target_gene_wise_ATAC_expr_mat) <- perturbations_to_try
            }

            target_gene_wise_ATAC_expr_mat <- na.omit(target_gene_wise_ATAC_expr_mat)
            message("Associated peak expression retrieved.")

            ###get the final score of the regulon
            if(nrow(target_gene_wise_ATAC_expr_mat) > 0){
                regulon_score <- colMeans(t(t(target_genes_nascent_expr_mat[rownames(target_gene_wise_ATAC_expr_mat),] * target_gene_wise_ATAC_expr_mat) * TF_steady_state_expr))
                message("Regulon score calculated.")
                return(regulon_score)
            }else{
                message("Regulon score has to be NA.")
                return(NA)
            }
        }else{
            return(NA)
        }        
    })
    names(TF_regulon_score_list) <- unique_motifs
    
    TF_regulon_score_mat <- do.call(rbind, na.omit(TF_regulon_score_list))
    rownames(TF_regulon_score_mat) <- unique_motifs
    colnames(TF_regulon_score_mat) <- perturbations_to_try
    
    return(TF_regulon_score_mat)
}