library("chromVAR")
library("dplyr")

### function for identifying relative enrichment/depletion of TF after perturbations
check_TF_scores_in_perturbations <- function(ATAC_meta, condition_to_use, perturbations_to_try, motif_score_mat){
    
    ###get unique programs
    uniq_motifs <- rownames(motif_score_mat)
    
    ###get cells to analyze in each bin
    cells_in_perturbations <- lapply(perturbations_to_try, function(x){
        
        cell_names <- (ATAC_meta %>% filter((sgRNA_identity == x) & (conditions == condition_to_use)))$ATAC_names_with_prefix
        
        return(cell_names)
    })
    names(cells_in_perturbations) <- perturbations_to_try
    
    ###get nascent UMI counts of cells in each bin
    motif_score_test <- lapply(uniq_motifs, function(z){
        
        message(z)
        program_NTC_score <- motif_score_mat[z, cells_in_perturbations[["NO-TARGET"]]]
        
        sgRNA_wise_distribution <- lapply(setdiff(perturbations_to_try, "NO-TARGET"), function(aa){
            
            message(aa)
            
            program_perturbation_score <- motif_score_mat[z, cells_in_perturbations[[aa]]]
            #message(sum(is.na(program_perturbation_score)))
            
            test_pval <- wilcox.test(x = program_NTC_score, y = program_perturbation_score)$p.value
            eff_size <- median(program_perturbation_score) - median(program_NTC_score)
            
            return(data.frame(sgRNA=aa, program=z, z_score_diff=eff_size, pval=test_pval))
        })
        sgRNA_wise_distribution_df <- do.call(rbind, sgRNA_wise_distribution)
        return(sgRNA_wise_distribution_df)
    })
    agg_motif_score_df <- do.call(rbind, motif_score_test)
    agg_motif_score_df$FDR <- p.adjust(agg_motif_score_df$p.value, method="BH")
    
    return(agg_motif_score_df)
}

### function for identifying enriched/depleted TFs with consistent gex trends 
load_TF_with_right_direction <- function(RNA_meta, condition, whole_tx_DE_path, motif_chromVAR_test_out_df){
    
    ###get qualified perturbations
    all_perturbations <- unique((RNA_meta %>% filter(conditions == condition))$sgRNA_identity)
    all_perturbations <- setdiff(all_perturbations, "NO-TARGET")
    
    ###for each perturbation
    perturbation_wise_active_TF <- lapply(all_perturbations, function(x){
        
        DE_motif_out <- motif_chromVAR_test_out_df %>% filter(sgRNA == x & z_score_diff > 0)
        
        DE_RNA_path <- paste0(whole_tx_DE_path, "/", x, "_whole.RDS")
        DE_RNA <- readRDS(DE_RNA_path) %>% filter((qval < 0.01) & (max.tissue == x) & (fold.change > 1.2) & (max.expr > 5))
        
        intersected_TF <- intersect(DE_motif_out$program, DE_RNA$gene_id)
        
        motif_out_for_FDR <- DE_motif_out %>% filter(program %in% intersected_TF)
        
        return(motif_out_for_FDR)
    })
    
    perturbation_wise_inactive_TF <- lapply(all_perturbations, function(x){
        
        DE_motif_out <- motif_chromVAR_test_out_df %>% filter(sgRNA == x & z_score_diff < 0)
        
        DE_RNA_path <- paste0(whole_tx_DE_path, "/", x, "_whole.RDS")
        DE_RNA <- readRDS(DE_RNA_path) %>% filter((qval < 0.01) & (max.tissue != x) & (fold.change > 1.2) & (max.expr > 5))
        
        intersected_TF <- intersect(DE_motif_out$program, DE_RNA$gene_id)
        
        motif_out_for_FDR <- DE_motif_out %>% filter(program %in% intersected_TF)
        
        return(motif_out_for_FDR)
    })
    
    perturbation_wise_active_TF_df <- do.call(rbind, perturbation_wise_active_TF)
    perturbation_wise_inactive_TF_df <- do.call(rbind, perturbation_wise_inactive_TF)
    total_df <- rbind(perturbation_wise_active_TF_df, perturbation_wise_inactive_TF_df)
    
    total_df$FDR <- p.adjust(total_df$pval, method = "BH")
    
    return(total_df)
}