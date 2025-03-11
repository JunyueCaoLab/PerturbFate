library(dplyr)
library(tidyr)
library(ggplot2)

### load and filter RNA cells
filter_perturbfate_cells_2lig_coassay <- function(file_directory, condition_table_directory, UMI_threshold, gene_threshold, unmatched_rate_threshold){
    
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

### change gene ids to gene symbols
gene_id2gene_names <- function(input_gene_id_vector, gtf_dir){
    #load the gtf file
    gtf <- rtracklayer::import(gtf_dir, format = 'gtf')
    gtf <- gtf[gtf$type == "gene"]
    id_conversion_table <- data.frame(gene_id = gtf$gene_id, symbol = gtf$gene_name)
    remove(gtf)

    id_conversion_table[match(input_gene_id_vector, id_conversion_table$gene_id), "symbol"]
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

### A wrapper function to load RNA, ATAC, sgRNA, and merge cells
perturbfate_load_dual_modal_data <- function(RNA_obj_directory, condition_table_directory, RNA_UMI_threshold, RNA_ngene_threshold,  RNA_unmatched_rate_threshold, ATAC_meta_path, sgRNA_count_matrix, ATAC_ncount_cutoff=500, ATAC_TSSE_cutoff=2, sgRNA_UMI_cutoff=5, sec_to_max_prop=0.3, max_to_all_prop=0.5){
    
    ###load RNA layer
    perturbfate_RNA <- filter_perturbfate_cells_2lig_coassay(file_directory = RNA_obj_directory, 
                                                   	     condition_table_directory = condition_table_directory, 
                                                   	     UMI_threshold = RNA_UMI_threshold, 
                                                   	     gene_threshold = RNA_ngene_threshold, 
                                                   	     unmatched_rate_threshold = RNA_unmatched_rate_threshold)
    
    message(paste0("RNA cells detected: ", nrow(perturbfate_RNA$filtered_merged_cells_bc_tb)))
    
    ###load ATAC layer
    perturbfate_ATAC <- readin_2lig_ATAC(ATAC_meta_path, ATAC_ncount_cutoff=ATAC_ncount_cutoff, ATAC_TSSE_cutoff=ATAC_TSSE_cutoff)
    
    message(paste0("ATAC cells detected: ", nrow(perturbfate_ATAC)))
    
    ###load sgRNA layer
    perturbfate_sgRNA <- sgRNA_count_matrix
    
    ###merge cells
    ###merge ATAC and RNA
    perturbfate_ATAC$RNA_ATAC_concensus_names <- gsub(x = perturbfate_ATAC$cell_names, pattern = "ATAC_", replacement = "")
    perturbfate_RNA$filtered_merged_cells_bc_tb$RNA_ATAC_concensus_names <- gsub(x = perturbfate_RNA$filtered_merged_cells_bc_tb$cell_names, pattern = "RNA_", replacement = "")
    RNA_ATAC_cocapture_df <- dplyr::inner_join(x = perturbfate_ATAC, y = perturbfate_RNA$filtered_merged_cells_bc_tb, by = "RNA_ATAC_concensus_names")
    
    message(paste0("Shared RNA and ATAC cells: ", nrow(RNA_ATAC_cocapture_df)))
    
    ###merge ATAC, RNA, and sgRNA
    intersected_cell_names <- intersect(colnames(sperturbfate_sgRNA), RNA_ATAC_cocapture_df$RNA_ATAC_concensus_names)
    sgRNA_RNA_ATAC_cocapture_df <- RNA_ATAC_cocapture_df[RNA_ATAC_cocapture_df$RNA_ATAC_concensus_names %in% intersected_cell_names,]

    sgRNA_RNA_ATAC_cocapture_df$total_sgRNA_UMI <- Matrix::colSums(perturbfate_sgRNA[,sgRNA_RNA_ATAC_cocapture_df$RNA_ATAC_concensus_names])
    sgRNA_RNA_ATAC_cocapture_df$max_sgRNA_UMI <- apply(X = perturbfate_sgRNA[,sgRNA_RNA_ATAC_cocapture_df$RNA_ATAC_concensus_names], MARGIN = 2, FUN = max)
    sgRNA_RNA_ATAC_cocapture_df$sec_sgRNA_UMI <- apply(perturbfate_sgRNA[, sgRNA_RNA_ATAC_cocapture_df$RNA_ATAC_concensus_names],2,function(x){
        sorted_counts <- sort(x, decreasing = T)
        return(sorted_counts[2])
    })
    sgRNA_RNA_ATAC_cocapture_df$sec_to_max_sgRNA_prop <- sgRNA_RNA_ATAC_cocapture_df$sec_sgRNA_UMI/sgRNA_RNA_ATAC_cocapture_df$max_sgRNA_UMI
    sgRNA_RNA_ATAC_cocapture_df$max_to_total_sgRNA_prop <- sgRNA_RNA_ATAC_cocapture_df$max_sgRNA_UMI/sgRNA_RNA_ATAC_cocapture_df$total_sgRNA_UMI
    
    message(paste0("Shared RNA, ATAC, sgRNA cells: ", nrow(sgRNA_RNA_ATAC_cocapture_df)))
    
    ###Identify perturbation identity
    ###check the primary sgRNA cells get
    clean_sgRNA_RNA_ATAC_df <- sgRNA_RNA_ATAC_cocapture_df %>% filter(total_sgRNA_UMI >= sgRNA_UMI_cutoff)
    perturbation_ident <- lapply(clean_sgRNA_RNA_ATAC_df$RNA_ATAC_concensus_names, function(x){
        sgRNA_expr_vector <- perturbfate_sgRNA[,x]
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
    
    ###add identities and remove doublets
    clean_sgRNA_RNA_ATAC_df$perturbation_identities <- perturbation_ident_vector
    clean_sgRNA_RNA_ATAC_df[clean_sgRNA_RNA_ATAC_df$sec_to_max_sgRNA_prop > sec_to_max_prop | 
                            clean_sgRNA_RNA_ATAC_df$max_to_total_sgRNA_prop < max_to_all_prop, "perturbation_identities"] <- "Mixed"
    
    n_singlets <- nrow(clean_sgRNA_RNA_ATAC_df %>% filter(perturbation_identities != "Mixed"))
    
    message(paste0("Final multi-modal singlets number: ", n_singlets))
    
    ###return the final project
    final_RNA_matrix <- perturbfate_RNA$filtered_merged_matrix[,clean_sgRNA_RNA_ATAC_df$cell_names.y]
    colnames(final_RNA_matrix) <- clean_sgRNA_RNA_ATAC_df$RNA_ATAC_concensus_names
    final_sgRNA_matrix <- perturbfate_sgRNA[,clean_sgRNA_RNA_ATAC_df$RNA_ATAC_concensus_names]
    
    out_list <- list(RNA_counts_mat = final_RNA_matrix,
                     sgRNA_counts_mat = final_sgRNA_matrix,
                     meta_data = clean_sgRNA_RNA_ATAC_df)
    
    return(out_list)
}