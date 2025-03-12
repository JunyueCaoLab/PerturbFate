library("rtracklayer")
library("GenomicRanges")
library("motifmatchr")
library("chromVAR")

### function for identifying candidate genes-associated ATAC peaks
find_genes_related_peaks <- function(figR_filtered_cis_df, input_gene_level_passed_gene_df, ATAC_peak_names, 
                                     hg38_gtf_input="/rugpfs/fs0/cao_lab/scratch/zxu/gencode_ref_anno/human_v38_050521/gencode.v38.primary_assembly.annotation.gtf"){
    
    ###get genes
    all_genes <- unique(input_gene_level_passed_gene_df$gene)
    
    ###get gtf and human promoters
    hg38_gtf <- rtracklayer::import(hg38_gtf_input, format = "gtf")
    hg38_protein_coding_genes <- hg38_gtf[hg38_gtf$type == "gene" & hg38_gtf$gene_type == "protein_coding"]
    
    ###only keep protein coding genes in the list too
    all_protein_coding_genes <- intersect(all_genes, hg38_protein_coding_genes$gene_name)
    
    ###get human promoter regions
    hg38_protein_coding_genes_promoters <- promoters(hg38_protein_coding_genes, upstream = 2000, downstream = 100)
    
    ###get peaks overlapping promoters of goi
    goi_promoters <- hg38_protein_coding_genes_promoters[hg38_protein_coding_genes_promoters$gene_name %in% all_protein_coding_genes]
    
    ###convert peakrange characters to granges
    peak_ranges <- from_char_to_GRanges(input_peak_char = ATAC_peak_names)
    
    ###intersect with promoter region of genes to get promoter regions
    peak_goi_overlap <- findOverlaps(query = peak_ranges, subject = goi_promoters)
    promoter_peak_gene_link_df <- data.frame(Overlapped_promoter_peaks=ATAC_peak_names[queryHits(peak_goi_overlap)],
                                    Gene=goi_promoters$gene_name[subjectHits(peak_goi_overlap)])
    
    ###get the rest of CRE peaks
    CRE_df <- figR_filtered_cis_df %>% filter((Gene %in% all_protein_coding_genes)&(!(PeakRanges %in% promoter_peak_gene_link_df$Overlapped_promoter_peaks)))
    
    ###combine into one df
    related_peaks <- data.frame(type=c(rep("Promoter", nrow(promoter_peak_gene_link_df)), rep("CRE", nrow(CRE_df))),
                                peaks=c(promoter_peak_gene_link_df$Overlapped_promoter_peaks, CRE_df$PeakRanges),
                                genes=c(promoter_peak_gene_link_df$Gene, CRE_df$Gene))
    return(related_peaks)
}

### function for scanning TF motifs within peaks associated to candidate genes
convert_peaks_to_motifs <- function(GRN_gene_df){
    
    ###get unique genes in the list
    unique_genes <- unique(GRN_gene_df$genes)
    
    ###convert peaks into granges
    GRN_peaks_granges <- from_char_to_GRanges(unique(GRN_gene_df$peaks))
    
    ###use chromVAR to convert
    ###convert peaks to motif counts
    motifs <- getJasparMotifs()
    motif_ix <- matchMotifs(motifs, GRN_peaks_granges, 
                            genome = BSgenome.Hsapiens.UCSC.hg38, out="scores")
    motif_counts_mat <- assays(motif_ix)[["motifCounts"]]
    binarized_motif_counts_mat <- motif_counts_mat
    binarized_motif_counts_mat@x[binarized_motif_counts_mat@x > 0] <- 1
    dense_binarized_motif_counts_mat <- as.matrix(binarized_motif_counts_mat)
    rownames(dense_binarized_motif_counts_mat) <- GRN_peaks_granges$names
    
    #return(dense_binarized_motif_counts_mat)
    
    ###get gene level motif count matrix
    gene_wise_motif_occur <- lapply(unique_genes, function(x){
        peaks_involved <- GRN_gene_df[GRN_gene_df$genes == x, "peaks"]
        if(length(peaks_involved) > 1){
            gene_level_motif_counts <- as.numeric(colSums(dense_binarized_motif_counts_mat[peaks_involved,]))
            gene_level_motif_counts[gene_level_motif_counts > 0] <- 1
            return(gene_level_motif_counts)
        }else{
            return(as.numeric(dense_binarized_motif_counts_mat[peaks_involved,]))
        }
    })
    gene_wise_motif_occur_mat <- do.call(rbind, gene_wise_motif_occur)
    colnames(gene_wise_motif_occur_mat) <- colnames(dense_binarized_motif_counts_mat)
    rownames(gene_wise_motif_occur_mat) <- unique_genes
    
    #return(gene_wise_motif_occur_mat)
    
    ###split complex
    TFs <- as.character(sapply(colnames(gene_wise_motif_occur_mat), function(x) strsplit(strsplit(x, split = "_", fixed = T)[[1]][2], split = "(", fixed = T)[[1]][1]))
    complexes_indices <- grep(pattern = "::", x = TFs)
    add_on_colnames_list <- lapply(TFs[complexes_indices], function(x) strsplit(x, split = "::", fixed = T)[[1]])
    names(add_on_colnames_list) <- TFs[complexes_indices]
    
    original_col_TFs <- TFs[-complexes_indices]                              
                                   
    add_on_TFs_list <- lapply(complexes_indices, function(x){
        return(add_on_colnames_list[[TFs[x]]])
    })                           
    add_on_TFs <- do.call(c, add_on_TFs_list)
                                   
    add_on_mat_list <- lapply(complexes_indices, function(x){
        n_elm <- length(add_on_colnames_list[[TFs[x]]])
        sub_mat_to_add <- gene_wise_motif_occur_mat[,rep(x, n_elm)]
        colnames(sub_mat_to_add) <- NULL
        return(sub_mat_to_add)
    })
    add_on_mat <- do.call(cbind, add_on_mat_list)                               
    new_peak_to_motif_mat <- cbind(gene_wise_motif_occur_mat[,-complexes_indices], add_on_mat)                               
    new_peak_to_motif_mat_colnames <- c(original_col_TFs, add_on_TFs)                           
                                   
    ###finally remove duplicated columns (maybe more than once)
    dup_TFs_boolean <- duplicated(new_peak_to_motif_mat_colnames)                               
    dup_TFs <- unique(new_peak_to_motif_mat_colnames[dup_TFs_boolean])
                                   
    combined_TF_cols_list <- lapply(dup_TFs, function(x){
        merged_TF_col <- rowSums(new_peak_to_motif_mat[,which(new_peak_to_motif_mat_colnames == x)])
        merged_TF_col[merged_TF_col > 0] <- 1
        return(merged_TF_col)
    })
    combined_TF_cols <- do.call(cbind, combined_TF_cols_list)                               
    final_TF_mat <- cbind(new_peak_to_motif_mat[,!(new_peak_to_motif_mat_colnames %in% dup_TFs)], combined_TF_cols)
    colnames(final_TF_mat) <- c(new_peak_to_motif_mat_colnames[!(new_peak_to_motif_mat_colnames %in% dup_TFs)], dup_TFs)                            
    rownames(final_TF_mat) <- unique_genes
    
    return(final_TF_mat)                                                             
}

