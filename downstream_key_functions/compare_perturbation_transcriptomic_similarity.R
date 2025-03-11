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