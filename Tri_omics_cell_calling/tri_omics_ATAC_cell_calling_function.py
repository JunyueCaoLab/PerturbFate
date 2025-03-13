import snapatac2 as snap
import subprocess
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse

def tri_omics_ATAC_tsv_to_peak_count_mat(tsv_file_paths, batch_h5ad_out_path, full_anndataset_h5ad_out_path, full_anndata_h5ad_out_path, ATAC_meta_out_path, batch_names, peaks_bed, chrom_sizes=snap.genome.hg38, cores=16):
    
    ###step1. load all batches of tsv
    tri_omics_ATAC = snap.pp.import_data(fragment_file=tsv_file_paths,
                                          file=batch_h5ad_out_path,
                                          chrom_sizes=chrom_sizes, 
                                          sorted_by_barcode=False, n_jobs=cores)
    
    ###step2. load into an anndataset
    tri_omics_ATAC_anndataset = snap.AnnDataSet(adatas=list(zip(batch_names, batch_h5ad_out_path)),
                                                 filename=full_anndataset_h5ad_out_path,backend=None)
    tri_omics_ATAC_anndataset.close()
    
    ###step3. make a anndata for downstream processing
    tri_omics_ATAC_anndataset = snap.read_dataset(full_anndataset_h5ad_out_path)
    tri_omics_ATAC_anndataset.obsm['fragment_single'] = tri_omics_ATAC_anndataset.adatas.obsm['fragment_single']
    tri_omics_ATAC_anndataset.obs['n_counts'] = tri_omics_ATAC_anndataset.adatas.obs['n_fragment']
    tri_omics_ATAC_anndataset.obs['n_mito_counts'] = tri_omics_ATAC_anndataset.adatas.obs['frac_mito']
    tri_omics_ATAC_anndata = tri_omics_ATAC_anndataset.to_adata()
    
    
    ###step4. make a peak count matrix
    snap.pp.make_peak_matrix(adata=tri_omics_ATAC_anndata, inplace=True, peak_file=peaks_bed)
    
    ###step5. calculate single-cell tsse
    snap.metrics.tsse(tri_omics_ATAC_anndata, snap.genome.hg38)
    
    ###step6. export files
    tri_omics_ATAC_anndata.obs.to_csv(ATAC_meta_out_path)
    tri_omics_ATAC_anndata.write_h5ad(full_anndata_h5ad_out_path)

 
