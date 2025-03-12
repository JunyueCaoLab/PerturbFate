library("dplyr")

### function for identifying steady-state and nascent gex
identify_closest_programs <- function(scaled_whole_tx_traj_mat, scaled_nascent_tx_traj_mat, whole_tx_clustering_df, nascent_tx_clustering_df){
    
    ###Get unique clusters
    whole_clusters <- unique(whole_tx_clustering_df[,1])
    nascent_clusters <- unique(nascent_tx_clustering_df[,1])
    
    ###For each cluster, get gene members and averaged values
    whole_averaged_program_val <- lapply(whole_clusters, function(x){
        gene_members <- rownames(whole_tx_clustering_df)[whole_tx_clustering_df[,1] == x]
        member_mat <- scaled_whole_tx_traj_mat[gene_members,]
        avg_val <- colMeans(member_mat)
        return(avg_val)
    })
    names(whole_averaged_program_val) <- whole_clusters
    
    nascent_averaged_program_val <- lapply(nascent_clusters, function(x){
        gene_members <- rownames(nascent_tx_clustering_df)[nascent_tx_clustering_df[,1] == x]
        member_mat <- scaled_nascent_tx_traj_mat[gene_members,]
        avg_val <- colMeans(member_mat)
        return(avg_val)
    })
    names(nascent_averaged_program_val) <- nascent_clusters
    
    ###check correlation between average values
    cor_mat <- matrix(0, nrow=length(whole_clusters), ncol=length(nascent_clusters))
    rownames(cor_mat) <- whole_clusters
    colnames(cor_mat) <- nascent_clusters
    
    for(whole_i in whole_clusters){
        for(nascent_i in nascent_clusters){
            cor_out <- cor(whole_averaged_program_val[[whole_i]], nascent_averaged_program_val[[nascent_i]])
            cor_mat[whole_i, nascent_i] <- cor_out
        }
    }
    
    return(cor_mat)
}

### function for identifying genes following the transcriptional dynamics
check_nascent_whole_temporal_gex <- function(whole_to_nascent_link_list, whole_search_area_list, whole_gene_program_df, nascent_gene_program_df, scaled_whole_tx_traj_mat, scaled_nascent_tx_traj_mat){
    
    ###remove colnames of matrices for better matching
    colnames(scaled_whole_tx_traj_mat) <- NULL
    colnames(scaled_nascent_tx_traj_mat) <- NULL
    
    ###get whole programs included
    unique_whole_programs <- unique(names(whole_to_nascent_link_list))
    
    ###for each slot, only keep genes
    whole_program_wise_genes_passed <- lapply(unique_whole_programs, function(x){
        nascent_partner_program <- whole_to_nascent_link_list[[x]]
        whole_gene_candidates <- rownames(whole_gene_program_df)[whole_gene_program_df[,1] %in% x]
        nascent_gene_candidates <- rownames(nascent_gene_program_df)[nascent_gene_program_df[,1] %in% nascent_partner_program]
        intersected_genes <- intersect(whole_gene_candidates, nascent_gene_candidates)
        
        ###get the bins to search
        search_area <- whole_search_area_list[[x]]
        refined_scaled_whole_tx_traj_mat <- scaled_whole_tx_traj_mat[intersected_genes, search_area]
        refined_scaled_nascent_tx_traj_mat <- scaled_nascent_tx_traj_mat[intersected_genes, search_area]
        
        ###get genes passed
        whole_max_bins <- lapply(1:nrow(refined_scaled_whole_tx_traj_mat), function(x){
            max_val <- max(refined_scaled_whole_tx_traj_mat[x,])
            max_bin <- which(refined_scaled_whole_tx_traj_mat[x,] == max_val)
            
            if(length(max_bin)>1){ 
                max_bin <- min(max_bin)
            }
            return(max_bin)
        })
        names(whole_max_bins) <- intersected_genes
        
        nascent_max_bins <- lapply(1:nrow(refined_scaled_nascent_tx_traj_mat), function(x){
            max_val <- max(refined_scaled_nascent_tx_traj_mat[x,])
            max_bin <- which(refined_scaled_nascent_tx_traj_mat[x,] == max_val)
            
            if(length(max_bin)>1){ 
                max_bin <- sample(max_bin, size = 1)
            }
            return(max_bin)
        })
        names(nascent_max_bins) <- intersected_genes
        
        ###compare bin val and only keep genes passed
        whole_nascent_bin_filtering <- lapply(intersected_genes, function(x){
            
            if(nascent_max_bins[[x]] <= whole_max_bins[[x]]){
                return(data.frame(gene=x, 
                                  nascent_max_bin=search_area[nascent_max_bins[[x]]], 
                                  whole_max_bin=search_area[whole_max_bins[[x]]]))
            }else{
                return(NA)
            }
        })
        ###reformat into a dataframe
        clean_whole_nascent_bin_filtering <- whole_nascent_bin_filtering[!sapply(whole_nascent_bin_filtering, function(x) all(is.na(x)))]
        passed_whole_nascent_bin_df <- as.data.frame(do.call(rbind, clean_whole_nascent_bin_filtering))
        passed_whole_nascent_bin_df$whole_gene_cluster <- x
        passed_whole_nascent_bin_df$nascent_gene_cluster <- nascent_gene_program_df[match(passed_whole_nascent_bin_df$gene, rownames(nascent_gene_program_df)), 1]
                                                                                 
        return(passed_whole_nascent_bin_df)                                                                         
    })
    final_out_df <- as.data.frame(do.call(rbind, whole_program_wise_genes_passed))
    return(final_out_df)                                                                     
}
