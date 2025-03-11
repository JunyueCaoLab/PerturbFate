library(dplyr)
library(tidyr)
library(ggplot2)

#### load single-cell nascent/pre-existing RNA and filter
filter_Trimodal_RNA_cells_2lig <- function(file_directory, condition_table_directory, UMI_threshold, gene_threshold, unmatched_rate_threshold){
    
    load(file_directory)
    df_cell$sample <- as.character(df_cell$sample)
    df_gene$gene_id <- as.character(df_gene$gene_id)
    rownames(gene_count_all) <- df_gene$gene_id
    colnames(gene_count_all) <- df_cell$sample

    #load condition table
    condition_table <- read.csv(condition_table_directory, header = T)

    #filter cells based on UMI-cutoff, gene-cutoff, and unmatched rate
    #this method is much much faster than sapply!!!
    binarized_gene_count_all <- gene_count_all
    binarized_gene_count_all@x[binarized_gene_count_all@x > 0] <- 1
    genes_per_cell <- Matrix::colSums(binarized_gene_count_all)
    #genes_per_cell <- sapply(1:ncol(gene_count_all), function(x) sum(gene_count_all[,x] > 0) )
    gene_based_whitelist <- names(genes_per_cell)[genes_per_cell > gene_threshold]
    UMI_based_whitelist <- df_cell[df_cell$UMI_count > UMI_threshold, "sample"]
    unmatched_rate_based_whitelist <- df_cell[df_cell$unmatched_rate < unmatched_rate_threshold, "sample"]

    whitelist <- intersect(gene_based_whitelist, intersect(UMI_based_whitelist, unmatched_rate_based_whitelist))

    whitelisted_df_cell <- df_cell[df_cell$sample %in% whitelist, ]
    whitelisted_gene_count_all <- gene_count_all[, whitelist]

    #split the cell names for recognition of different conditions
    separated_CB_df <- tidyr::separate(data = data.frame(full_name = whitelist), col = "full_name", into = c('PCR_group', 'lig_CBs'), sep = '\\.', remove = F)
    separated_CB_df$lig_cb1 <- substr(separated_CB_df$lig_CBs, start = 1, stop = 10)
    separated_CB_df$lig_cb2 <- substr(separated_CB_df$lig_CBs, start = 11, stop = 20)

    #separate cells in different groups
    UMI_counts <- c()
    gene_number <- c()
    conditions <- c()
    cell_names <- c()
    
    for (each_condition in unique(condition_table$Conditions)){
    available_lig_bc1 <- condition_table[condition_table$Conditions == each_condition, "barcode"]
    cells_in_this_condition <- separated_CB_df[separated_CB_df$lig_cb1 %in% available_lig_bc1, "full_name"]
    PCR_group <- separated_CB_df[separated_CB_df$lig_cb1 %in% available_lig_bc1, "PCR_group"]                             
    UMI_counts <- append(UMI_counts, whitelisted_df_cell[whitelisted_df_cell$sample %in% cells_in_this_condition, "UMI_count"])
    gene_number <- append(gene_number, genes_per_cell[cells_in_this_condition])                               
    conditions <- append(conditions, rep(each_condition, length(cells_in_this_condition)))
    cell_names <- append(cell_names, cells_in_this_condition)}
    
    UMI_per_group_df <- data.frame(cell_names = cell_names, 
                                   UMI_counts = UMI_counts, 
                                   gene_number = gene_number,
                                   conditions = conditions)
    
    p <- ggplot(UMI_per_group_df, 
                   aes(x=conditions, y=UMI_counts, fill=conditions))+geom_violin()+geom_boxplot(width=0.2)                              

    list(filtered_merged_matrix = whitelisted_gene_count_all,
       filtered_merged_cells_bc_tb = UMI_per_group_df,
       UMI_per_condition_plot = p)
}

#### change gene ids to gene symbols
gene_id2gene_names <- function(input_gene_id_vector, gtf_dir){
    #load the gtf file
    gtf <- rtracklayer::import(gtf_dir, format = 'gtf')
    gtf <- gtf[gtf$type == "gene"]
    id_conversion_table <- data.frame(gene_id = gtf$gene_id, symbol = gtf$gene_name)
    remove(gtf)

    id_conversion_table[match(input_gene_id_vector, id_conversion_table$gene_id), "symbol"]
}

### merge nascent and preexisting tx of single cells
merge_expr_matrix <- function(matrix1, matrix2, ordered_cell_names1, ordered_cell_names2, concensus_cell_name_list){
    
    n_cell = ncol(matrix1)
    
    ###get the concensus gene list
    concensus_genes <- unique(c(rownames(matrix1), rownames(matrix2)))
    matrix1_gene_append <- setdiff(concensus_genes, rownames(matrix1))
    matrix2_gene_append <- setdiff(concensus_genes, rownames(matrix2))
    
    if(length(matrix1_gene_append) > 0){
        append_mat_for_1 <- Matrix::sparseMatrix(i=c(1), j = c(1), x=c(0), dims = c(length(matrix1_gene_append),n_cell))
        rownames(append_mat_for_1) <- matrix1_gene_append
        new_mat1 <- rbind(matrix1, append_mat_for_1)
        new_mat1 <- new_mat1[concensus_genes, ordered_cell_names1]
    }else{
        new_mat1 <- matrix1
        new_mat1 <- new_mat1[concensus_genes, ordered_cell_names1]
    }
    
    if(length(matrix2_gene_append) > 0){
        append_mat_for_2 <- Matrix::sparseMatrix(i=c(1), j = c(1), x=c(0), dims = c(length(matrix2_gene_append),n_cell))
        rownames(append_mat_for_2) <- matrix2_gene_append
        new_mat2 <- rbind(matrix2, append_mat_for_2)
        new_mat2 <- new_mat2[concensus_genes, ordered_cell_names2]
    }else{
        new_mat2 <- matrix2
        new_mat2 <- new_mat2[concensus_genes, ordered_cell_names2]
    }
    
    merged_matrix <- new_mat1 + new_mat2
    colnames(merged_matrix) <- concensus_cell_name_list
    return(list(new_matrix_out1 = new_mat1, 
                new_matrix_out2 = new_mat2, 
                merged_matrix = merged_matrix))
}


merge_nascent_preexisting_tx_Trimodal_sc <- function(old_tx_obj, nascent_tx_obj, additional_old_condition_table_col_2include, additional_old_cols_rename, additional_new_condition_table_col_2include, additional_new_cols_rename, link_table){
    
    new_nascent_tx_obj <- nascent_tx_obj
    new_nascent_tx_obj$filtered_merged_matrix <- new_nascent_tx_obj$filtered_merged_matrix[,new_nascent_tx_obj$filtered_merged_cells_bc_tb$cell_names]
    
    ###consolidate cell names between old and new RNA data
    for(row_i in 1:nrow(link_table)){
        nascent_suffix <- link_table[row_i, "nascent_tx"]
        old_suffix <- link_table[row_i, "old_tx"]
        cells_indices <- grep(x = new_nascent_tx_obj$filtered_merged_cells_bc_tb$cell_names, pattern = nascent_suffix)
        cell_names_to_replace <- new_nascent_tx_obj$filtered_merged_cells_bc_tb[cells_indices, "cell_names"]
        new_nascent_tx_obj$filtered_merged_cells_bc_tb[cells_indices, "cell_names"] <- gsub(x = cell_names_to_replace, pattern = nascent_suffix, replacement = old_suffix)
    }
    
    ###change names in the nascent count matrix
    colnames(new_nascent_tx_obj$filtered_merged_matrix) <- new_nascent_tx_obj$filtered_merged_cells_bc_tb$cell_names
    
    ###merge both data
    intersected_cell_names <- intersect(new_nascent_tx_obj$filtered_merged_cells_bc_tb$cell_names,
                                        old_tx_obj$filtered_merged_cells_bc_tb$cell_names)
    
    ###get an integrated condition table
    integrated_condition_table <- old_tx_obj$filtered_merged_cells_bc_tb[match(intersected_cell_names, old_tx_obj$filtered_merged_cells_bc_tb$cell_names), 
                                                                         c(c("cell_names", "UMI_counts", "gene_number", "conditions"), additional_old_condition_table_col_2include)]
    colnames(integrated_condition_table)[match(additional_old_condition_table_col_2include, colnames(integrated_condition_table))] <- additional_old_cols_rename
    colnames(integrated_condition_table)[colnames(integrated_condition_table) == "UMI_counts"] <- "old_UMI_counts"
    integrated_condition_table$new_UMI_counts <- new_nascent_tx_obj$filtered_merged_cells_bc_tb[match(intersected_cell_names, new_nascent_tx_obj$filtered_merged_cells_bc_tb$cell_names), "UMI_counts"]
    integrated_condition_table[, additional_new_cols_rename] <- new_nascent_tx_obj$filtered_merged_cells_bc_tb[match(intersected_cell_names, new_nascent_tx_obj$filtered_merged_cells_bc_tb$cell_names), 
                                                                                                                          additional_new_condition_table_col_2include]
    
    rownames(integrated_condition_table) <- integrated_condition_table$cell_names
    
    ###also get a raw nascent ratio
    integrated_condition_table$raw_nascent_ratio <- integrated_condition_table$new_UMI_counts/(integrated_condition_table$new_UMI_counts+integrated_condition_table$old_UMI_counts)
    
    ###integrate nascent and old matrices
    old_mat <- old_tx_obj$filtered_merged_matrix[,integrated_condition_table$cell_names]
    new_mat <- new_nascent_tx_obj$filtered_merged_matrix[,integrated_condition_table$cell_names]

    ####considering the gene list might have difference, create a concensus matrix
    concensus_mat_list <- merge_expr_matrix(old_mat, new_mat, 
                                          integrated_condition_table$cell_names, 
                                          integrated_condition_table$cell_names, 
                                          integrated_condition_table$cell_names)
    
    output_obj <- list(condition_table = integrated_condition_table,
                       new_old_merged_count_matrix = concensus_mat_list$merged_matrix,
                       old_only_count_matrix = concensus_mat_list$new_matrix_out1,
                       new_only_count_matrix = concensus_mat_list$new_matrix_out2) 
    
    return(output_obj)                                         
}

### read in ATAC cells identified by snapATAC2
readin_2lig_ATAC <- function(ATAC_meta_from_snapATAC2_path, ATAC_ncount_cutoff, ATAC_TSSE_cutoff){
    
    ###reconstruct a meta table
    ATAC_meta <- (read.csv(ATAC_meta_from_snapATAC2_path))[,2:8]

    ###rename data frame columns                
    colnames(ATAC_meta) <- c("sample_batch", "ATAC_counts", "ATAC_Mito_counts", "ATAC_TSSE", "ATAC_cell_names") 

    ###set up cutoff
    ATAC_meta <- ATAC_meta %>% filter(ATAC_counts >= ATAC_ncount_cutoff & ATAC_TSSE >= ATAC_TSSE_cutoff)
            
    return(ATAC_meta)
}

### load sgRNA of cells and reformat cell names
gRNA_cell_reformatting <- function(input_gRNA_summary_rdata_dir, reformatting_df, gRNA_PCR_condition_name_2bused = 'No'){
  
  load(input_gRNA_summary_rdata_dir)
  input_gRNA_count_mat <- gRNA_count
  
  #the first step: unify cell name format between gRNA lib and whole txme lib
  gRNA_cell_names_df <- data.frame(gRNA_cell_names=colnames(input_gRNA_count_mat)) %>% tidyr::separate(col="gRNA_cell_names", into=c("gRNA_PCR_group", "bc"), sep="\\.", remove=FALSE)
  
  if (gRNA_PCR_condition_name_2bused != "No"){
    gRNA_cell_names_df <- dplyr::inner_join(x = gRNA_cell_names_df, y = reformatting_df, by = "gRNA_PCR_group")
    new_gRNA_count_mat <- input_gRNA_count_mat[,gRNA_cell_names_df$gRNA_cell_names]
    
  } else {
    gRNA_cell_names_df <- dplyr::left_join(x = gRNA_cell_names_df, y = reformatting_df, by = "gRNA_PCR_group")
    gRNA_cell_names_df <- na.omit(gRNA_cell_names_df)  
    new_gRNA_count_mat <- input_gRNA_count_mat[,gRNA_cell_names_df$gRNA_cell_names]
  }
  new_cell_names <- paste(gRNA_cell_names_df$whole_txme_PCR_group, gRNA_cell_names_df$bc, sep = ".")
  colnames(new_gRNA_count_mat) <- new_cell_names
  new_gRNA_count_mat
}

### the function for identifying sgRNA identity
sgRNA_singlets_identification_stringency <- function(sgRNA_RNA_ATAC_cocapture_df, sgRNA_expr_matrix, lower_sgRNA_UMI_cutoff=5, max_to_all_prop=0.5, sec_to_max_prop=0.3){
    
    ###remove cells with low sgRNA UMI
    clean_sgRNA_RNA_ATAC_df <- sgRNA_RNA_ATAC_cocapture_df %>% filter(total_sgRNA_UMI >= lower_sgRNA_UMI_cutoff)
    
    ###get sgRNA expr info
    max_sgRNA_UMI <- sapply(clean_sgRNA_RNA_ATAC_df$concensus_cell_names, function(x){
        return(max(sgRNA_expr_matrix[,x]))
    })
    sec_sgRNA_UMI <- sapply(clean_sgRNA_RNA_ATAC_df$concensus_cell_names, function(x){
        sgRNA_UMI_vector <- sgRNA_expr_matrix[,x]
        return(sgRNA_UMI_vector[order(sgRNA_UMI_vector, decreasing = T)][2])
        })
    clean_sgRNA_RNA_ATAC_df$max_sgRNA_UMI <- max_sgRNA_UMI
    clean_sgRNA_RNA_ATAC_df$sec_sgRNA_UMI <- sec_sgRNA_UMI
    clean_sgRNA_RNA_ATAC_df$max_to_total_sgRNA_prop <- clean_sgRNA_RNA_ATAC_df$max_sgRNA_UMI/clean_sgRNA_RNA_ATAC_df$total_sgRNA_UMI
    clean_sgRNA_RNA_ATAC_df$sec_to_max_sgRNA_prop <- clean_sgRNA_RNA_ATAC_df$sec_sgRNA_UMI/clean_sgRNA_RNA_ATAC_df$max_sgRNA_UMI
    
    ###check the primary sgRNA cells get
    perturbation_ident <- parallel::mclapply(clean_sgRNA_RNA_ATAC_df$concensus_cell_names, mc.cores=8, function(x){
        sgRNA_expr_vector <- sgRNA_expr_matrix[,x]
        max_sgRNA_UMI <- max(sgRNA_expr_vector)
        max_sgRNA_name <- names(sgRNA_expr_vector)[sgRNA_expr_vector == max_sgRNA_UMI]
        
        if(length(max_sgRNA_name) > 1){
            return("Mixed")
        }else{
            sgRNA_target_gene <- strsplit(max_sgRNA_name, split = "_", fixed = TRUE)[[1]][1]
            return(sgRNA_target_gene)
        }
    })
    perturbation_ident_vector <- do.call(c, perturbation_ident)
    
    ###add identities
    clean_sgRNA_RNA_ATAC_df$sgRNA_identity <- perturbation_ident_vector
    clean_sgRNA_RNA_ATAC_df[clean_sgRNA_RNA_ATAC_df$max_to_total_sgRNA_prop < max_to_all_prop | 
                            clean_sgRNA_RNA_ATAC_df$sec_to_max_sgRNA_prop > sec_to_max_prop, "sgRNA_identity"] <- "Mixed"
    
    return(clean_sgRNA_RNA_ATAC_df) 
}

### merge multi-modal cells
merge_Trimodal_data <- function(newRNA_obj, oldRNA_obj, new_old_link_table, ATAC_meta, sgRNA_mat, lower_sgRNA_UMI_cutoff=5, max_to_all_prop=0.4, sec_to_max_prop=0.3, ATAC_ncount_cutoff=200, ATAC_TSSE_cutoff=1.5){
    
    ###merge new and old RNA first
    new_old_RNA_merge <- merge_nascent_preexisting_tx_Trimodal_sc(old_tx_obj = oldRNA_obj, 
                                                             nascent_tx_obj = newRNA_obj, 
                                                             link_table = new_old_link_table,
                                                             additional_new_cols_rename = c(),
                                                             additional_new_condition_table_col_2include = c(),
                                                             additional_old_cols_rename = c(),
                                                             additional_old_condition_table_col_2include = c())
    message(paste0("new+preexisting RNA cells: ", nrow(new_old_RNA_merge$condition_table)))
    
    ###merge ATAC data
    RNA_meta <- new_old_RNA_merge$condition_table
    ATAC_meta <- ATAC_meta %>% filter(ATAC_counts >= ATAC_ncount_cutoff & ATAC_TSSE >= ATAC_TSSE_cutoff)
    RNA_meta$concensus_cell_names <- gsub(x = RNA_meta$cell_names, pattern = "RNA_", replacement = "")
    ATAC_meta$concensus_cell_names <- gsub(x = ATAC_meta$cell_names, pattern = "ATAC_", replacement = "")
    combined_meta <- dplyr::inner_join(x = RNA_meta, y = ATAC_meta, by = "concensus_cell_names")
    message(paste0("new+old RNA+ATAC cells: ", nrow(combined_meta)))
    
    ###merge sgRNA data
    intersected_cell_names <- intersect(combined_meta$concensus_cell_names, colnames(sgRNA_mat))
    combined_meta <- combined_meta %>% filter(concensus_cell_names %in% intersected_cell_names)
    intersected_sgRNA_mat <- sgRNA_mat[,intersected_cell_names]
    combined_meta$total_sgRNA_UMI <- Matrix::colSums(intersected_sgRNA_mat[,combined_meta$concensus_cell_names])
    message(paste0("new+preexisting RNA+ATAC+sgRNA cells: ", nrow(combined_meta)))
    
    ###identify sgRNA singlets
    singlets_meta <- sgRNA_singlets_identification_stringency(sgRNA_RNA_ATAC_cocapture_df = combined_meta, sgRNA_expr_matrix = intersected_sgRNA_mat, 
                                                             lower_sgRNA_UMI_cutoff = lower_sgRNA_UMI_cutoff, 
                                                             max_to_all_prop = max_to_all_prop, sec_to_max_prop = sec_to_max_prop)
    n_singlets <- sum(singlets_meta$sgRNA_identity != "Mixed")
    message(paste0("sgRNA singlets: ", n_singlets))
    
    ###return singlet data
    singlets_meta <- singlets_meta %>% filter(sgRNA_identity != "Mixed")
    
    newRNA_counts_mat <- new_old_RNA_merge$new_only_count_matrix[,singlets_meta$cell_names.x]
    preexistingRNA_counts_mat <- new_old_RNA_merge$old_only_count_matrix[,singlets_meta$cell_names.x]
    colnames(newRNA_counts_mat) <- singlets_meta$concensus_cell_names
    colnames(preexistingRNA_counts_mat) <- singlets_meta$concensus_cell_names
    
    final_sgRNA_mat <- intersected_sgRNA_mat[,singlets_meta$concensus_cell_names]
    
    
    ###reorganize the meta data
    reformatted_singlets_meta_df <- data.frame(concensus_cell_names=singlets_meta$concensus_cell_names,
                                               conditions=singlets_meta$Conditions,
                                               ATAC_cell_names=singlets_meta$cell_names.y,
                                               RNA_cell_names=singlets_meta$cell_names.x,
                                               ATAC_read_counts=singlets_meta$ATAC_counts,
                                               ATAC_TSSE=singlets_meta$ATAC_TSSE,
                                               newRNA_UMI_counts=singlets_meta$new_UMI_counts,
                                               preexistingRNA_UMI_counts=singlets_meta$old_UMI_counts,
                                               RNA_gene_num=singlets_meta$gene_number,
                                               sgRNA_total_UMI=singlets_meta$total_sgRNA_UMI,
                                               primary_sgRNA_UMI=singlets_meta$max_sgRNA_UMI,
                                               second_sgRNA_UMI=singlets_meta$sec_sgRNA_UMI,
                                               top_to_total_sgRNA_prop=singlets_meta$max_to_total_sgRNA_prop,
                                               sec_to_top_sgRNA_prop=singlets_meta$sec_to_max_sgRNA_prop,
                                               sgRNA_identity=singlets_meta$sgRNA_identity,
                                               first_lig_plate=singlets_meta$Plate_ID,
                                               first_lig_row=singlets_meta$row,
                                               first_lig_col=singlets_meta$col,
                                               first_lig_well=singlets_meta$order_bc)
    rownames(reformatted_singlets_meta_df) <- reformatted_singlets_meta_df$concensus_cell_names                                                 
    
    ###get the output list
    out_list <- list(newRNA_count_mat = newRNA_counts_mat,
                     preexistingRNA_count_mat = preexistingRNA_counts_mat,
                     sgRNA_count_mat = final_sgRNA_mat,
                     meta.data = reformatted_singlets_meta_df)
    
    return(out_list)
}
