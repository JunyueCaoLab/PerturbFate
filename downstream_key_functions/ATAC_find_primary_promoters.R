library("Matrix")
library("rtracklayer")
library("parallel")
library("GenomicRanges")

### the function for changing gene ids to gene symbols
gene_id2gene_names_cross_versions <- function(input_gene_id_vector, gtf_dir="/file/to/path/human_v38/gencode.v38.primary_assembly.annotation.gtf"){
    #load the gtf file
    gtf <- rtracklayer::import(gtf_dir, format = 'gtf')
    gtf <- gtf[gtf$type == "gene"]
    gene_id_no_version <- sapply(gtf$gene_id, function(x) strsplit(x, split = ".", fixed = T)[[1]][1])
    id_conversion_table <- data.frame(gene_id = gene_id_no_version, symbol = gtf$gene_name)
    remove(gtf)
                                 
    input_gene_id_no_version <- sapply(input_gene_id_vector, function(x) strsplit(x, split = ".", fixed = T)[[1]][1])                             

    id_conversion_table[match(input_gene_id_no_version, id_conversion_table$gene_id), "symbol"]
}

### the function for finding the primary promoter regions of genes
ATAC_promoter_identification <- function(perturbfate_ATAC_promoter_counts_obj, sgRNA_alignment_bed, sgRNA_identity_colname="sgRNA_identity", cell_cutoff=30){
    
    ###get all genes and split promoters into gene categories
    all_genes_in_mat_list <- parallel::mclapply(perturbfate_ATAC_promoter_counts_obj$promoter_annotation_df$gene, mc.cores = 4, function(x) 
                                                strsplit(x = x, split = ",")[[1]])
    all_genes_in_mat <- do.call(c, all_genes_in_mat_list)                                            
    all_genes_in_mat <- unique(all_genes_in_mat)
    message("Gene id are transformed.")
    
    ###identify promoters for each gene and return a list
    gene_wise_promoter <- parallel::mclapply(all_genes_in_mat, mc.cores=16, function(a){
        indices <- grep(pattern = a, x = perturbfate_ATAC_promoter_counts_obj$promoter_annotation_df$gene)
        return(perturbfate_ATAC_promoter_counts_obj$promoter_annotation_df$peak[indices]) 
    })
    names(gene_wise_promoter) <- all_genes_in_mat
    message("Gene specific promoter list is identified.")
    
    ###make it into a GRanges object
    gene_wise_promoter_df_list <- parallel::mclapply(names(gene_wise_promoter), mc.cores=16, function(b){
        promoter_peaks <- gene_wise_promoter[[b]]
        promoter_peak_subdf <- perturbfate_ATAC_promoter_counts_obj$promoter_annotation_df %>% filter(peak %in% promoter_peaks)
        promoter_peak_subdf$gene_id <- b
        return(promoter_peak_subdf)
    })
    gene_wise_promoter_df <- as.data.frame(do.call(rbind, gene_wise_promoter_df_list))
    gene_wise_promoter_gr <- GenomicRanges::GRanges(seqnames = gene_wise_promoter_df$chr, 
                                                    ranges = IRanges(start = gene_wise_promoter_df$start, end = gene_wise_promoter_df$end))
    gene_wise_promoter_gr$gene_id <- gene_wise_promoter_df$gene_id
    gene_wise_promoter_gr$peak <- gene_wise_promoter_df$peak
    message("Gene specific promoter GRange object is made.")
    
    ###convert gene id to gene symbols
    gene_wise_promoter_gr$gene_symbol <- gene_id2gene_names_cross_versions(input_gene_id_vector = gene_wise_promoter_gr$gene_id)
    gene_wise_promoter_gr <- gene_wise_promoter_gr[!is.na(gene_wise_promoter_gr$gene_symbol)]
    message("Gene id to gene symbol conversion is made.")
    
    ###overlap sgRNA alignment with promoter regions
    sgRNA_alignment <- rtracklayer::import(sgRNA_alignment_bed, format = "bed")
    sgRNA_alignment$gene_symbol <- sapply(sgRNA_alignment$name, function(x) strsplit(x, split = "_", fixed = T)[[1]][1])
    sgRNA_promoter_overlap <- findOverlaps(query = sgRNA_alignment, subject = gene_wise_promoter_gr)                                      
    
    ###only keep promoters with at least sgRNA targeting it 
    sgRNA_promoter_overlap_df <- as.data.frame(sgRNA_promoter_overlap)
                                          
    gene_matching <- sgRNA_alignment$gene_symbol[sgRNA_promoter_overlap_df$queryHits] == gene_wise_promoter_gr$gene_symbol[sgRNA_promoter_overlap_df$subjectHits]                                      
                                          
    sgRNA_matched_promoters <- gene_wise_promoter_gr[unique(sgRNA_promoter_overlap_df$subjectHits[gene_matching])]                                                                         
    message("Only promoters with sgRNA targeted are kept.")                                   
                                          
    ###get normalized counts of promoters of each population
    cell_pop <- table(perturbfate_ATAC_promoter_counts_obj$meta.data[,sgRNA_identity_colname])
    qualified_cell_pop <- names(cell_pop)[cell_pop > cell_cutoff]
    ###only include perturbation which have corresponding promoters                                       
    qualified_cell_pop <- intersect(qualified_cell_pop, unique(sgRNA_matched_promoters$gene_symbol))                                      
    qualified_cell_pop_w_ntc <- c(qualified_cell_pop, "NO-TARGET")
                                          
    ###get normalized counts for qualified cell pops
    norm_pop_level_in_prom_counts <- parallel::mclapply(qualified_cell_pop_w_ntc, mc.cores=8, function(x){
        all_cells_in_pop <- perturbfate_ATAC_promoter_counts_obj$meta.data[perturbfate_ATAC_promoter_counts_obj$meta.data[,sgRNA_identity_colname]==x, "ATAC_cell_names"]
        accum_in_prom_counts <- Matrix::rowSums(perturbfate_ATAC_promoter_counts_obj$ATAC_in_promoter_Mat[,all_cells_in_pop])
        accum_total_counts <- sum(perturbfate_ATAC_promoter_counts_obj$meta.data[match(all_cells_in_pop, perturbfate_ATAC_promoter_counts_obj$meta.data$ATAC_cell_names), "ATAC_read_counts"])
        norm_accum_in_prom_counts <- accum_in_prom_counts*1e6/accum_total_counts
        return(norm_accum_in_prom_counts)
    })
    names(norm_pop_level_in_prom_counts) <- qualified_cell_pop_w_ntc
    message("Normalized ATAC counts of promoters are retrieved.")
    
    #return(norm_pop_level_in_prom_counts)                                      
                                          
    ###choose promoters with highest ATAC signal from each gene
    gene_npromoter <- table(sgRNA_matched_promoters$gene_symbol)                                     
    genes_with_multiple_promoters <- names(gene_npromoter)[gene_npromoter > 1]
    genes_with_one_promoters <- names(gene_npromoter)[gene_npromoter == 1]
                                          
    highest_promoter_annotation_df1 <- data.frame(region=sgRNA_matched_promoters$peak[sgRNA_matched_promoters$gene_symbol %in% genes_with_one_promoters],
                                                  gene_symbol=sgRNA_matched_promoters$gene_symbol[sgRNA_matched_promoters$gene_symbol %in% genes_with_one_promoters])                                      
                                          
    highest_atac_promoter <- parallel::mclapply(genes_with_multiple_promoters, mc.cores=8, function(x){
        
        peak_regions <- sgRNA_matched_promoters$peak[sgRNA_matched_promoters$gene_symbol == x]

        signals <- norm_pop_level_in_prom_counts[["NO-TARGET"]][peak_regions]
        max_signal <- max(signals)
        if(max(signals)==0){
            return(NA)
        }
        highest_region <- names(signals)[signals == max_signal]
        if(length(highest_region) > 1){
            return(sample(highest_region, size = 1))
        }else{
            return(highest_region)
        }
    })                                   
                                          
    highest_atac_promoter_vector <- do.call(c, highest_atac_promoter)
                                          
    highest_promoter_annotation_df2 <- data.frame(region=highest_atac_promoter_vector, 
                                                  gene_symbol=genes_with_multiple_promoters)                                         
    highest_promoter_annotation_df <- as.data.frame(rbind(highest_promoter_annotation_df1, 
                                                    highest_promoter_annotation_df2))                                      
    
    rownames(highest_promoter_annotation_df) <- highest_promoter_annotation_df$gene_symbol
    message("The main promoter of each gene is chosen.")
                                          
    ###Calculate the log2FC
    calculable_genes <- intersect(highest_promoter_annotation_df$gene_symbol, qualified_cell_pop)
    regions <- highest_promoter_annotation_df[calculable_genes, "region"]
    
    atac_logfc_kd_to_ctrl <- lapply(calculable_genes, function(x){
        fc <- norm_pop_level_in_prom_counts[[x]][regions]/norm_pop_level_in_prom_counts[["NO-TARGET"]][regions]
        return(fc)
    })
    atac_fc_kd_to_ctrl_mat <- do.call(cbind, atac_logfc_kd_to_ctrl)
                                   
    colnames(atac_fc_kd_to_ctrl_mat) <- calculable_genes
    rownames(atac_fc_kd_to_ctrl_mat) <- calculable_genes
    message("FC calculation is done.")
                                          
    return(list(ATAC_fc_mat=atac_fc_kd_to_ctrl_mat,
                promoter_info=highest_promoter_annotation_df))
}