library(ggplot2)
library(dplyr)

### function for processing bulk CRISPR screen raw counts exported from the upstream script
bulk_CRISPR_screen_libraries_processing <- function(raw_counts_url_vector, condition_names){
    
    ###for each library
    every_lib <- lapply(1:length(raw_counts_url_vector), function(xx){
        x <- raw_counts_url_vector[xx]
        raw_lib <- read.csv(x, header=T)
        
        depth_norm_counts <- (raw_lib$Lib_counts)*1e6/sum(raw_lib$Lib_counts)
        
        ###get NTC count sum
        NTC_count_sum <- sum(depth_norm_counts[raw_lib$gene_symbol == "NO-TARGET"])
        
        norm_NTC_rel_counts <- depth_norm_counts/NTC_count_sum
        out_subdf <- data.frame(gene_symbol=raw_lib$gene_symbol,
                                sgRNA_id=raw_lib$id,
                                norm_NTC_rel_counts=norm_NTC_rel_counts)
        colnames(out_subdf)[3] <- paste0(condition_names[xx], "_norm_NTC_rel_counts")
        return(out_subdf)
    })
    
    base_df <- every_lib[[1]]
    
    for(condition_i in 2:length(raw_counts_url_vector)){
        base_df <- dplyr::inner_join(x=base_df, y=every_lib[[condition_i]][,2:ncol(every_lib[[condition_i]])], by="sgRNA_id")
    }
    
    return(base_df)
}