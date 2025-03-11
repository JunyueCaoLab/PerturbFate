###identify phenotypic changes caused by perturbation
pseudotime_perturbation_wise_check <- function(meta_with_pseudotime, min_cells=20, NTC_name="NO-TARGET"){
    
    ###get qualified populations
    qualified_sgRNA <- names(table(meta_with_pseudotime$sgRNA_identity))[table(meta_with_pseudotime$sgRNA_identity) > min_cells]
    qualified_sgRNA_wo_NTC <- setdiff(qualified_sgRNA, NTC_name)
    
    ###get the distribution
    NTC_dist <- (meta_with_pseudotime %>% filter(sgRNA_identity == NTC_name))$corrected_pseudotime
    
    ###compare
    pseudotime_comparison_df_list <- lapply(qualified_sgRNA_wo_NTC, function(x){
        KD_dist <- (meta_with_pseudotime %>% filter(sgRNA_identity == x))$corrected_pseudotime
        ks_test_obj <- ks.test(x = NTC_dist, y = KD_dist)
        median_NTC <- median(NTC_dist)
        median_KD <- median(KD_dist)
        eff_size <- median_KD - median_NTC
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
    
    pseudotime_comparison_df <- as.data.frame(do.call(rbind, pseudotime_comparison_df_list))
    pseudotime_comparison_df$FDR <- p.adjust(pseudotime_comparison_df$pval, method = "BH")
    
    return(pseudotime_comparison_df)
}

###identify relative enrichments of perturbations compared to NTC
state_specific_cell_enrichment <- function(cell_meta_with_cell_state){
    
    ###get all cell states
    all_cell_states <- unique(cell_meta_with_cell_state$cell_state)
    
    ###for each cell state, perform proportion test
    all_cells <- nrow(cell_meta_with_cell_state)
    all_NTC_cells <- sum(cell_meta_with_cell_state$sgRNA_identity == "NO-TARGET")
    
    prop_test <- lapply(all_cell_states, function(x){
        state_sub_df <- cell_meta_with_cell_state %>% filter(cell_state == x)
        state_cells <- nrow(state_sub_df)
        state_NTC_cells <- sum(state_sub_df$sgRNA_identity == "NO-TARGET")
        
        state_sgRNA <- unique(state_sub_df$sgRNA_identity)
        sgRNA_prop_test <- lapply(state_sgRNA, function(y){
            
            all_sgRNA_cells <- sum(cell_meta_with_cell_state$sgRNA_identity == y)
            state_sgRNA_cells <- sum(state_sub_df$sgRNA_identity == y)
            
            ###get test significance
            pval <- (prop.test(c(state_sgRNA_cells, state_NTC_cells), 
                              c(all_sgRNA_cells, all_NTC_cells), alternative = "greater", correct=T))$p.value
            enrich_FC <- (state_sgRNA_cells/all_sgRNA_cells)/(state_NTC_cells/all_NTC_cells)
            return(data.frame(pval=pval,
                              enrich_FC=enrich_FC,
                              sgRNA=y))
        })
        sgRNA_prop_test_df <- do.call(rbind, sgRNA_prop_test)
        sgRNA_prop_test_df$cell_state <- x
        return(sgRNA_prop_test_df)
    })
    all_test_df <- as.data.frame(do.call(rbind, prop_test))
    
    return(all_test_df)
}

### identify relative state shift of perturbed cells to that of NTC
compare_pseudotime_distribution_between_perturbations_and_conditions <- function(meta_df_with_both_conditions_and_phenptype_score, baseline_condition="DMSO_treated_cells", trt_condition="PLX_treated_cells ", NTC_sgRNA="NO-TARGET"){
    
    ###get sgRNA identities
    baseline_sgRNAs <- unique(meta_df_with_both_conditions_and_phenptype_score %>% filter((conditions == baseline_condition)))$sgRNA_identity
    trt_sgRNAs <- unique(meta_df_with_both_conditions_and_phenptype_score %>% filter((conditions == trt_condition)))$sgRNA_identity
    sgRNA_keep <- intersect(baseline_sgRNAs, trt_sgRNAs)
    
    ###go through all sgRNAs
    rel_pseudotime <- lapply(sgRNA_keep, function(x){
        
        pseudotime_baseline <- (meta_df_with_both_conditions_and_phenptype_score %>% filter((conditions == baseline_condition) & (sgRNA_identity == x)))$corrected_pseudotime
        mean_pseudotime_baseline <- mean(pseudotime_baseline)
        pseudotime_trt <- (meta_df_with_both_conditions_and_phenptype_score %>% filter((conditions == trt_condition) & (sgRNA_identity == x)))$corrected_pseudotime
        
        rel_pseudotime_shift <- pseudotime_trt - mean_pseudotime_baseline
        
        return(rel_pseudotime_shift)
    })
    names(rel_pseudotime) <- sgRNA_keep
    
    ###compare all other sgRNAs with NTC
    sgRNA_wo_NTC <- setdiff(sgRNA_keep, NTC_sgRNA)
    
    compare_pseudotime <- lapply(sgRNA_wo_NTC, function(x){
        NTC_rel_pseudotime <- rel_pseudotime[[NTC_sgRNA]]
        target_rel_pseudotime <- rel_pseudotime[[x]]
        
        pval <- (wilcox.test(NTC_rel_pseudotime, target_rel_pseudotime))$p.value
        eff_size <- mean(target_rel_pseudotime) - mean(NTC_rel_pseudotime)
        
        return(data.frame(sgRNA = x, eff_size = eff_size, pval = pval))
    })
    
    pseudotime_compare_df <- as.data.frame(do.call(rbind, compare_pseudotime))
    pseudotime_compare_df$FDR <- p.adjust(pseudotime_compare_df$pval, method = "BH")
    
    return(pseudotime_compare_df)
}