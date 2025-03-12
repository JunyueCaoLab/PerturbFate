library("Matrix")
library("dplyr")

### function for normalizing single-cell nascent RNA
normalize_nascent_counts_all_batches <- function(input_tri_omics_obj, perturbations_to_keep, NTC_name="NO-TARGET", baseline_condition="DMSO_treated_cells"){
    
    ###keep perturbations passing the sciCAR filter
    meta_trimmed <- input_tri_omics_obj$meta.data %>% filter(sgRNA_identity %in% perturbations_to_keep)
    newRNA_mat <- input_tri_omics_obj$newRNA_count_mat[,meta_trimmed$concensus_cell_names]
    
    ###get the scaling factor    
    baseline_median_nascent_UMI <- sapply(perturbations_to_keep, function(y){
        medium_UMI <- median((meta_trimmed %>% filter((sgRNA_identity == y) & (conditions == baseline_condition)))$newRNA_UMI_counts)
        return(medium_UMI)
    })
    names(baseline_median_nascent_UMI) <- perturbations_to_keep
    baseline <- baseline_median_nascent_UMI[NTC_name]
    baseline_ratio <- baseline_median_nascent_UMI/baseline

    treatment_median_nascent_UMI <- sapply(perturbations_to_keep, function(y){
        medium_UMI <- median((meta_trimmed %>% filter((sgRNA_identity == y) & (conditions != baseline_condition)))$newRNA_UMI_counts)
        return(medium_UMI)
    })
    names(treatment_median_nascent_UMI) <- perturbations_to_keep
    treatment_ratio <- treatment_median_nascent_UMI/baseline
    
    scaling_factor_list <- list(baseline_ratio, treatment_ratio)
    
    ###make a ref table
    treatment_condition <- setdiff(unique(meta_trimmed$conditions), baseline_condition)
    condition_vector <- c(baseline_condition, treatment_condition)
    
    ref_df_list <- lapply(1:2, function(m){
        cond <- condition_vector[m]
        scaling_factor <- scaling_factor_list[[m]]
        cod_perturbation <- paste0(cond, "_", names(scaling_factor))
        ref_df <- data.frame(condition_perturbation=cod_perturbation,
                             scaling_factor=scaling_factor)
    })
    ref_df <- as.data.frame(do.call(rbind, ref_df_list))

    ###add scaling factor into meta data
    meta_trimmed$condition_perturbation <- paste0(meta_trimmed$conditions, "_", meta_trimmed$sgRNA_identity)
    meta_trimmed <- dplyr::left_join(x = meta_trimmed, y = ref_df, by = "condition_perturbation")

    ###normalize the newCount matrix
    rep_time_vector <- diff(newRNA_mat@p)
    
    norm_factor_for_x <- lapply((1:nrow(meta_trimmed)), function(a){
        scaling_factor_a <- meta_trimmed$scaling_factor[a]
        rep_time_a <- rep(scaling_factor_a, rep_time_vector[a])
        return(rep_time_a)
    })
    norm_factor_for_x <- do.call(c, norm_factor_for_x)
    
    total_UMI_for_x <- sapply((1:nrow(meta_trimmed)), function(a){
        newRNA_UMI_a <- meta_trimmed$newRNA_UMI_counts[a]
        rep_time_a <- rep(newRNA_UMI_a, rep_time_vector[a])
        return(rep_time_a)
    })
    total_UMI_for_x <- do.call(c, total_UMI_for_x)
    
    log1p_norm_newRNA_mat_x <- log1p((newRNA_mat@x*1e5*norm_factor_for_x)/total_UMI_for_x)
    norm_newRNA_mat_x <- (newRNA_mat@x*1e5*norm_factor_for_x)/total_UMI_for_x

    log_norm_newRNA_mat <- Matrix::sparseMatrix(i = newRNA_mat@i, p = newRNA_mat@p, x = log1p_norm_newRNA_mat_x, 
                                                dims = c(nrow(newRNA_mat), ncol(newRNA_mat)))
    norm_newRNA_mat <- Matrix::sparseMatrix(i = newRNA_mat@i, p = newRNA_mat@p, x = norm_newRNA_mat_x, 
                                                dims = c(nrow(newRNA_mat), ncol(newRNA_mat)))
    
    colnames(log_norm_newRNA_mat) <- colnames(newRNA_mat)
    rownames(log_norm_newRNA_mat) <- rownames(newRNA_mat)
    
    colnames(norm_newRNA_mat) <- colnames(newRNA_mat)
    rownames(norm_newRNA_mat) <- rownames(newRNA_mat)

    return(list(log_norm_newRNA_mat=log_norm_newRNA_mat,
                norm_newRNA_mat=norm_newRNA_mat,
                meta_table=meta_trimmed))
}

### function for identifying global nascent transcription changes upon perturbation
synthesis_top_genes_perturbation_wise_check.v2 <- function(Tri_omics_meta, log_transformed_nascent_expr_mat, min_cells, condition_to_choose, NTC_name="NO-TARGET", top_n_genes = 1000){
    
    ###get qualified populations
    conditional_meta <- Tri_omics_meta %>% filter((conditions == condition_to_choose) & (concensus_cell_names %in% colnames(log_transformed_nascent_expr_mat)))
    
    qualified_sgRNA <- names(table(conditional_meta$sgRNA_identity))[table(conditional_meta$sgRNA_identity) > min_cells]
    qualified_sgRNA_wo_NTC <- setdiff(qualified_sgRNA, NTC_name)
    
    ###get topN genes from NTC
    NTC_cell_names <- (conditional_meta %>% filter(sgRNA_identity == NTC_name))$concensus_cell_names
    NTC_all_expr <- Matrix::rowMeans(log_transformed_nascent_expr_mat[, NTC_cell_names])
    topN_genes <- names(NTC_all_expr)[order(NTC_all_expr, decreasing = T)][1:top_n_genes]
    NTC_dist <- NTC_all_expr[topN_genes]
    
    ###compare
    newRNA_UMI_comparison_df_list <- lapply(qualified_sgRNA_wo_NTC, function(x){
        KD_cell_names <- (conditional_meta %>% filter(sgRNA_identity == x))$concensus_cell_names
        KD_all_expr <- Matrix::rowMeans(log_transformed_nascent_expr_mat[, KD_cell_names])
        KD_dist <- KD_all_expr[topN_genes]
        ks_test_obj <- wilcox.test(x = NTC_dist, y = KD_dist)
        mean_NTC <- mean(NTC_dist)
        mean_KD <- mean(KD_dist)
        eff_size <- mean_KD - mean_NTC
        if(eff_size < 0){
            return(data.frame(sgRNA = x,
                       eff_size = eff_size,
                       direction = "Down",
                       pval = ks_test_obj$p.value))
        }else{
            return(data.frame(sgRNA = x,
                       eff_size = eff_size,
                       direction = "Up",
                       pval = ks_test_obj$p.value))
        }
    })
    
    newRNA_UMI_comparison_df <- as.data.frame(do.call(rbind, newRNA_UMI_comparison_df_list))
    newRNA_UMI_comparison_df$FDR <- p.adjust(newRNA_UMI_comparison_df$pval, method = "BH")
    
    return(newRNA_UMI_comparison_df)
}