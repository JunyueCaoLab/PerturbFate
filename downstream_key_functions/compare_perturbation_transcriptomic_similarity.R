library("dplyr")

### function for comparing the transcriptomic similarity between different 
choose_perturbation_and_compare_similarity <- function(DEG_list, all_gene_list, perturbation_choose, top_ngenes=200, NTC_names="NO-TARGET"){
    
    ###get all shared genes included in perturbations mentioned
    shared_genes_list <- lapply(perturbation_choose, function(x){
        genes <- rownames(all_gene_list[[x]])
        return(genes)
    })
    shared_genes_vector <- do.call(c, shared_genes_list)
    gene_occurence <- table(shared_genes_vector)
    gene_keep <- names(gene_occurence)[gene_occurence == length(perturbation_choose)]
    
    ###get top DE genes for each perturbation
    topDE_genes <- lapply(perturbation_choose, function(x){
        
        DEG_df <- DEG_list[[x]]
        DEG_KD_up <- DEG_df %>% filter((max.tissue == x))
        DEG_NTC_up <- DEG_df %>% filter((max.tissue == NTC_names))
        n_KD_up_genes <- nrow(DEG_KD_up)
        n_NTC_up_genes <- nrow(DEG_NTC_up)
        
        if(n_KD_up_genes <= top_ngenes){
            KD_up_genes_selected <- rownames(DEG_KD_up)
        }else{
            DEG_KD_up_ordered <- DEG_KD_up[order(DEG_KD_up$fold.change, decreasing = T),]
            KD_up_genes_selected <- rownames(DEG_KD_up)[1:top_ngenes]
        }
        if(n_NTC_up_genes <= top_ngenes){
            NTC_up_genes_selected <- rownames(DEG_NTC_up)
        }else{
            DEG_NTC_up_ordered <- DEG_NTC_up[order(DEG_NTC_up$fold.change, decreasing = T),]
            NTC_up_genes_selected <- rownames(DEG_NTC_up)[1:top_ngenes]
        }
        
        return(c(KD_up_genes_selected, NTC_up_genes_selected))
    })
    all_top_genes <- unique(do.call(c, topDE_genes))
    all_top_genes_included <- intersect(all_top_genes, gene_keep)
    
    ###get FC of these top DE genes 
    topDE_FC_list <- lapply(perturbation_choose, function(x){
        topDE_genes_perturbation <- all_gene_list[[x]] %>% filter(gene_short_name %in% all_top_genes_included)
        topDE_genes_perturbation[topDE_genes_perturbation$max.tissue == NTC_names, "fold.change"] <- 1/topDE_genes_perturbation[topDE_genes_perturbation$max.tissue == NTC_names, "fold.change"]
        return(topDE_genes_perturbation$fold.change)
    })
    topDE_FC_mat <- do.call(cbind, topDE_FC_list)
    colnames(topDE_FC_mat) <- perturbation_choose
    rownames(topDE_FC_mat) <- all_top_genes_included
    
    return(topDE_FC_mat)
}

### function for calculating the pseudobulk FC between perturbations and NTC 
calculate_pseudobulk_FC_gene_features <- function(merged_perturbfate_obj, gene_signature, condition_to_choose="DMSO_treated_cells"){
    
    ###get all sgRNA
    sgRNAs_included <- unique(merged_perturbfate_obj$meta.data$sgRNA_identity)
    pseudobulk_expr_list <- lapply(sgRNAs_included, function(x){
        sgRNA_cells <- (merged_perturbfate_obj$meta.data %>% filter(sgRNA_identity==x & conditions == condition_to_choose))$concensus_cell_names
        pseudobulk_expr <- rowSums(merged_perturbfate_obj$RNA_count_mat[,sgRNA_cells])
        norm_pseudobulk <- (pseudobulk_expr * 1e6)/sum(pseudobulk_expr)
        return(norm_pseudobulk)
    })
    pseudobulk_expr_mat <- do.call(cbind, pseudobulk_expr_list)
    rownames(pseudobulk_expr_mat) <- rownames(merged_perturbfate_obj$RNA_count_mat)
    colnames(pseudobulk_expr_mat) <- sgRNAs_included
    
    ###only keep genes in the signature
    pseudobulk_expr_mat <- pseudobulk_expr_mat[intersect(rownames(pseudobulk_expr_mat), gene_signature),]
    
    binarized_pseudobulk_expr_mat <- apply(pseudobulk_expr_mat, 1, function(mm){mm != 0})
    
    trimmed_pseudobulk_expr_mat <- pseudobulk_expr_mat[rowSums(binarized_pseudobulk_expr_mat) > ncol(binarized_pseudobulk_expr_mat)*2/3,]
    trimmed_NTC_expr <- as.numeric(trimmed_pseudobulk_expr_mat[,"NO-TARGET"])
    trimmed_pseudobulk_expr_mat_no_NTC <- trimmed_pseudobulk_expr_mat[,-match("NO-TARGET", colnames(trimmed_pseudobulk_expr_mat))]
    trimmed_pseudobulk_expr_mat_no_NTC_FC <- log2((trimmed_pseudobulk_expr_mat_no_NTC+0.1)/(trimmed_NTC_expr+0.1))
    
    return(trimmed_pseudobulk_expr_mat_no_NTC_FC)
}

### function for checking the overlap of DEGs with a gene set across perturbations
check_DEGs_overlap_and_overall_correlations <- function(sgRNAs_to_compare, input_DEG_list, ref_geneset_log2FC_df, merged_perturbfate_obj, gene_signature, condition_to_choose="DMSO_treated_cells"){
    
    ###check n DEGs overlap
    n_DEGs_overlap <- sapply(sgRNAs_to_compare, function(x){
        sgRNA_spec_DEGs <- rownames(input_DEG_list$DEG_list[[x]])
        return(length(intersect(sgRNA_spec_DEGs, gene_signature_loose)))
     })
    names(n_DEGs_overlap) <- sgRNAs_to_compare
    
    ###check correlation with the signature
    all_perturbations_to_NTC_gene_log2FC <- calculate_pseudobulk_FC_gene_features(merged_perturbfate_obj, gene_signature, condition_to_choose=condition_to_choose)
    
    ###get corr coef between perturbations and ref KD log2FC profile
    intersected_gene_signature <- intersect(rownames(all_perturbations_to_NTC_gene_log2FC), gene_signature)
    sgRNA_KD_log2FC_subdf <- all_perturbations_to_NTC_gene_log2FC[intersected_gene_signature, sgRNAs_to_compare]
    ref_log2FC_avg <- rowMeans(ref_geneset_log2FC_df[intersected_gene_signature,])
    
    corr_coef <- sapply(sgRNAs_to_compare, function(z){
        sgRNA_spec_DEGs <- cor(as.numeric(sgRNA_KD_log2FC_subdf[,z]), ref_log2FC_avg)
        return(sgRNA_spec_DEGs)
     })
    names(corr_coef) <- sgRNAs_to_compare
    
    ###make the output df
    output_df <- data.frame(sgRNA=sgRNAs_to_compare, n_DEGs_overlap=n_DEGs_overlap, corr_coef=corr_coef)
    
    return(output_df)
}