# PerturbFate
A combinatorial-indexing single-cell platform for concurrent profiling of chromatin accessibility, nascent and pre-existing RNA, or steady-state transcriptomic phenotypes, along with sgRNA identities

# PerturbFate Dual-omics, Tri-omics, and Bulk Screen Pipelines

---

## Dual-omics

The `Dual-omics` folder contains the main preprocessing scripts and daughter scripts that process FASTQ files into crude count matrices.

### RNA Pipeline
**Main script:**  
`Dual_omics_processing_scripts/RNA_modality/perturbFate_dual_omics_RNA.main.sh`

**General workflow**
1. Fastq renaming  
2. RNA molecular identifier recognition and UMI/cell barcode extraction  
3. Adapter trimming  
4. STAR reference genome mapping  
5. Filter mapped reads  
6. Deduplication based on cell barcode, UMI, and genomic coordinates  
7. Single-cell SAM file splitting  
8. Feature counting and matrix reformatting  

**Required inputs**
1. Fastq files named by PCR well identifiers  
2. Fastq sample ID file — single-column `.txt` without header  
3. Ligation barcode error-correction dictionaries (≤1 mismatch)  
4. Reference genome index for STAR  
5. Genomic annotation GTF for feature counting  

---

### ATAC Pipeline
**Main script:**  
`Dual_omics_processing_scripts/ATAC_modality/perturbFate_dual_omics_ATAC.main.sh`

**General workflow**
1. Fastq renaming  
2. ATAC mosaic end sequence recognition and cell barcode extraction  
3. Adapter trimming  
4. STAR reference genome mapping  
5. Filter mapped reads  
6. PCR-well-level deduplication  
7. BAM merging and conversion  
8. TSV creation and indexing for snapATAC2  

**Required inputs**
1. Fastq files named by PCR well identifiers  
2. Fastq sample ID file — single-column `.txt` without header  
3. Ligation barcode error-correction dictionaries (≤1 mismatch)  
4. Reference genome index for STAR  
5. Chromosome length file — `.txt` with chromosome names and lengths (must match reference genome)  

---

### sgRNA Pipeline
**Main script:**  
`Dual_omics_processing_scripts/sgRNA_modality/perturbFate_dual_omics_sgRNA.main.sh`

**General workflow**
1. Fastq renaming  
2. Amplicon constant region matching  
3. sgRNA UMI parsing  
4. Count matrix reformatting  

**Required inputs**
1. Fastq files named by PCR well identifiers  
2. Fastq sample ID file — single-column `.txt` without header  
3. Ligation barcode error-correction dictionaries (≤1 mismatch)  
4. Inner i7 barcode error-correction dictionary (≤1 mismatch)  
5. sgRNA annotation file — `.csv` with columns: `gRNA_seq`, `names` (sgRNA identifiers), `gene_symbol`  
6. sgRNA correction file — dictionary of sgRNA sequences (≤1 mismatch)  

---

### Dual-omics Cell Calling

The `Dual_omics_cell_calling` folder contains key functions for parsing count matrices across modalities.

#### ATAC Cell Calling
**Script:**  
`Dual_omics_cell_calling/dual_omics_ATAC_cell_calling_function.py`

**General workflow**
1. Read SnapATAC2 fragment files  
2. Combine batches into a unified AnnData object  
3. Perform peak counting  
4. Calculate TSS enrichment scores  
5. Export outputs  

**Required inputs**
1. `tsv.gz` — outputs from ATAC preprocessing  
2. `peaks_bed` — ATAC peak BED file  

#### Cross-modality Matching
**Script:**  
`Dual_omics_cell_calling/dual_omics_RNA_ATAC_sgRNA_matching_function.py`

**General workflow**
1. Read single-cell RNA gene count matrix  
2. Match treatment conditions split across ligation wells via barcodes  
3. Read ATAC metadata  
4. Read sgRNA count matrix  
5. Convert sgRNA cell names to consensus using barcode matching (`inner i7 + sgRNA i7 ↔ RNA/ATAC i7`)  
6. Return matched cell names and outputs  

**Required inputs**
1. RNA count matrix (`.RDS`) from RNA preprocessing  
2. ATAC metadata (`.csv`) from ATAC cell calling  
3. sgRNA count matrix (`.RDS`) from sgRNA preprocessing  
4. `condition_table` (`.csv`) with columns: `Plate_ID`, `order_bc`, `row`, `col`, `barcode`, `Conditions`  

---

## Tri-omics

The `Tri-omics` folder contains preprocessing scripts for three-modality data.  
Workflows and input requirements are similar to Dual-omics.

**Notes**
- Nascent RNA is sequenced separately and processed as an independent RNA modality.  
- Barcode matching with RNA/ATAC happens downstream.  

### Tri-omics Cell Calling
The `Tri_omics_cell_calling` folder provides functions for cross-modality parsing.

**Notes**
- Requires an additional mapping file: nascent RNA i7 ↔ preexisting ATAC i7.  
- Enables integration of nascent RNA with RNA and ATAC modalities.  

---

## Bulk Screen Processing

The `bulk_screen_processing` folder contains scripts for preprocessing and normalizing PerturbFate bulk screen data.

### Preprocessing
**Script:**  
`bulk_screen_processing/perturbFate_dual_omics_sgRNA.main.sh`

**General workflow**
1. Construct dual-sgRNA mismatch dictionary  
2. Match i5, i7, and constant regions in FASTQs  
3. Assign dual sgRNAs from PE reads and count perturbation reads  

**Required inputs**
1. Fastq files containing bulk screen reads  
2. i5 barcodes (10bp) to distinguish samples  
3. `seq_gene_table` (`.csv`) containing sgRNA annotations  
   - Required columns: `sgRNA_1_seq`, `sgRNA_2_seq`  
   - Example: `bulk_screen_processing/perturbfate_dual_KD_library_with_NTC.csv`  

### Downstream Normalization
**Script:**  
`bulk_screen_processing/bulk_screen_downstream_processing.R`

**General workflow**
1. Normalize counts by sequencing depth (factor = 1e6)  
2. Normalize to NO-TARGET (divide each sgRNA by sum of NO-TARGET counts)  

**Required input**
- Raw count table from preprocessing  

---

## Downstream Key Functions

### Identify Primary Promoter Regions
**Script:**  
`downstream_key_functions/ATAC_find_primary_promoters.R`

**General workflow**
1. Identify all promoters per gene  
2. Overlap sgRNA alignments with promoters  
3. Select promoter with highest ATAC coverage in NO-TARGET cells  

**Required inputs**
- `perturbfate_ATAC_promoter_counts_obj`  
- `sgRNA_alignment_bed`  
- `sgRNA_identity_colname` (default `"sgRNA_identity"`)  

---

### Identify Perturbation-driven State Shifts
**Script:**  
`downstream_key_functions/identify_perturbation_state_changes.R`

**General workflow**
1. Filter sgRNAs with > `min_cells`  
2. Extract NTC pseudotime distribution  
3. For each sgRNA:  
   - Compare KD vs NTC pseudotime with KS test  
   - Compute effect size and direction  
4. Combine results and adjust p-values (BH FDR)  

**Required input**
- `meta_with_pseudotime` — pseudotime metadata  

---

### Normalize Nascent RNA and Assess Global Changes
**Script:**  
`downstream_key_functions/nascent_RNA_normalization_and_effect_check.R`

**General workflow**
1. Compute median nascent UMI per sgRNA (baseline and treatment)  
2. Derive scaling factors relative to NTC and join to metadata  
3. Normalize sparse matrix and create log1p-normalized version  
4. From NTC, select top expressed genes  
5. For each sgRNA:  
   - Compute mean expression of top genes  
   - Compare vs NTC with Wilcoxon test  
   - Record effect size/direction, adjust p-values  

**Required inputs**
- Raw nascent single-cell count matrix  
- Perturbation identity metadata  

---

