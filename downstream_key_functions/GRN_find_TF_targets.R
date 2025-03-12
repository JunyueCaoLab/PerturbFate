
library("dplyr")

### function for identifying close motif and nascent gex programs
identify_closest_programs.nascent_motif <- function(scaled_nascent_tx_traj_mat, scaled_motif_traj_mat, nascent_tx_clustering_df, motif_clustering_df){
    
    ###Get unique clusters
    nascent_clusters <- unique(nascent_tx_clustering_df[,1])
    motif_clusters <- unique(motif_clustering_df[,1])
    
    ###For each cluster, get gene members and averaged values
    nascent_averaged_program_val <- lapply(nascent_clusters, function(x){
        gene_members <- rownames(nascent_tx_clustering_df)[nascent_tx_clustering_df[,1] == x]
        message(length(gene_members))
        member_mat <- scaled_nascent_tx_traj_mat[gene_members,]
        avg_val <- colMeans(member_mat)
        return(avg_val)
    })
    names(nascent_averaged_program_val) <- nascent_clusters
    
    motif_averaged_program_val <- lapply(motif_clusters, function(x){
        gene_members <- rownames(motif_clustering_df)[motif_clustering_df[,1] == x]
        member_mat <- scaled_motif_traj_mat[gene_members,]
        if(length(gene_members) > 1){
            avg_val <- colMeans(member_mat)
        }else{
            avg_val <- as.numeric(member_mat)
        }
        return(avg_val)
    })
    names(motif_averaged_program_val) <- motif_clusters

    ###check correlation between average values
    cor_mat <- matrix(0, nrow=length(nascent_clusters), ncol=length(motif_clusters))
    rownames(cor_mat) <- nascent_clusters
    colnames(cor_mat) <- motif_clusters
    
    for(nascent_i in nascent_clusters){
        for(motif_i in motif_clusters){
            cor_out <- cor(nascent_averaged_program_val[[nascent_i]], motif_averaged_program_val[[motif_i]])
            cor_mat[nascent_i, motif_i] <- cor_out
        }
    }
    
    return(cor_mat)
}

### function for identifying motif-nascent gex pairs 
check_motif_nascent_temporal_gex <- function(gene_to_motif_mat, passed_gene_df, TF_to_nascent_link_list, search_area_list, motif_program_df, nascent_gene_program_df, scaled_motif_traj_mat, scaled_nascent_tx_traj_mat){
    
    colnames(scaled_motif_traj_mat) <- NULL
    colnames(scaled_nascent_tx_traj_mat) <- NULL
    
    ###get whole programs included
    unique_TF_motif_programs <- unique(names(TF_to_nascent_link_list))
    
    ###for each slot, only keep genes
    TF_program_wise_genes_passed <- lapply(unique_TF_motif_programs, function(x){
        nascent_partner_program <- TF_to_nascent_link_list[[x]]
        TF_motif_candidates <- rownames(motif_program_df)[motif_program_df[,1] %in% x]
        ###only keep genes that passed the gene level filtering
        nascent_gene_candidates <- intersect(passed_gene_df$gene, rownames(nascent_gene_program_df)[nascent_gene_program_df[,1] %in% nascent_partner_program])
        message(TF_motif_candidates)
        message(length(nascent_gene_candidates))
        ###get the bins to search
        search_area <- search_area_list[[x]]
        refined_scaled_motif_traj_mat <- scaled_motif_traj_mat[TF_motif_candidates, search_area]
        refined_scaled_nascent_tx_traj_mat <- scaled_nascent_tx_traj_mat[nascent_gene_candidates, search_area]
        
        ###search through each TF-gene pair: if the gene contains motif from the TF, then consider it
        TF_wise_gene_scan <- lapply(TF_motif_candidates, function(y){
            
            max_val <- max(refined_scaled_motif_traj_mat[y,])
            max_bin <- which(refined_scaled_motif_traj_mat[y,] == max_val)
            message(max_bin)
            
            if(length(max_bin)>1){ 
                max_bin <- min(max_bin)
            }
            return(max_bin)
        })
        names(TF_wise_gene_scan) <- TF_motif_candidates
        
        nascent_max_bins <- lapply(nascent_gene_candidates, function(z){
            max_val <- max(refined_scaled_nascent_tx_traj_mat[z,])
            max_bin <- which(refined_scaled_nascent_tx_traj_mat[z,] == max_val)
            
            if(length(max_bin)>1){ 
                max_bin <- sample(max_bin, size = 1)
            }
            return(max_bin)
        })
        names(nascent_max_bins) <- nascent_gene_candidates
        
        out_list=list()
        
        for(each_TF in TF_motif_candidates){
            message(each_TF)
            motif_spec_max_bin <- TF_wise_gene_scan[[each_TF]]
            message(motif_spec_max_bin)
            genes_with_motif <- rownames(gene_to_motif_mat)[as.logical(gene_to_motif_mat[,each_TF])]
            nascent_gene_candidate_with_motif <- intersect(nascent_gene_candidates, genes_with_motif)
            message(length(nascent_gene_candidate_with_motif))
            
            if(length(nascent_gene_candidate_with_motif)>0){
                
                for(each_nascent_gene in nascent_gene_candidate_with_motif){
                    nascent_spec_max_bin <- nascent_max_bins[[each_nascent_gene]]
                    
                    if(motif_spec_max_bin <= nascent_spec_max_bin){
                        out_list[[paste0(each_TF, "_", each_nascent_gene)]] <- data.frame(motif=each_TF,
                                                                                          gene=each_nascent_gene,
                                                                                          motif_max_bin=search_area[motif_spec_max_bin],
                                                                                          nascent_gex_max_bin=search_area[nascent_spec_max_bin])
                    }
                }
            }
        }
        
        ###reformat into a dataframe
        out_df <- as.data.frame(do.call(rbind, out_list))
        out_df$TF_motif_cluster <- x
        out_df$nascent_gene_cluster <- nascent_gene_program_df[match(out_df$gene, rownames(nascent_gene_program_df)), 1]
                                                                                 
        return(out_df)                                                                         
    })
    
    final_out_df <- as.data.frame(do.call(rbind, TF_program_wise_genes_passed))
    return(final_out_df)                                                                     
}