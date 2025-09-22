# PerturbFate
A combinatorial-indexing single-cell platform for concurrent profiling of chromatin accessibility, nascent and pre-existing RNA, or steady-state transcriptomic phenotypes, along with sgRNA identities

# Dual-omics and Tri-omics Preprocessing Pipelines

The `Dual-omics` and `Tri-omics` folders contain the main preprocessing scripts and daughter scripts that process FASTQ files into crude count matrices. Each modality (RNA, ATAC, sgRNA) has its own pipeline, and cell-calling functions are provided to unify count matrices across modalities.

---

## RNA Pipeline

**Main script:**  
`Dual_omics_processing_scripts/RNA_modality/perturbFate_dual_omics_RNA.main.sh`

### General Workflow
1. Fastq renaming  
2. RNA molecular identifier recognition and UMI/cell barcode extraction  
3. Adapter trimming  
4. STAR reference genome mapping  
5. Filter mapped reads  
6. Deduplication based on cell barcode, UMI, and genomic coordinates  
7. Single-cell SAM file splitting  
8. Feature counting and matrix reformatting  

### Required Inputs
1. **Fastq files** – named by identifiers of each PCR well  
2. **Fastq sample ID file** – single-column `.txt` file (no header) with the FASTQ prefix of each sequenced PCR well  
3. **Ligation barcode error correction dictionaries** – containing ligation barcodes with at most one mismatch  
4. **Reference genome index** for STAR  
5. **Genomic annotation GTF** for feature counting  

---

## ATAC Pipeline

**Main script:**  
`Dual_omics_processing_scripts/ATAC_modality/perturbFate_dual_omics_ATAC.main.sh`

### General Workflow
1. Fastq renaming  
2. ATAC mosaic end sequence recognition and cell barcode extraction  
3. Adapter trimming  
4. STAR reference genome mapping  
5. Filter mapped reads  
6. PCR-well-level deduplication  
7. BAM merging and conversion  
8. TSV creation and indexing for snapATAC2 input  

### Required Inputs
1. **Fastq files** – named by identifiers of each PCR well  
2. **Fastq sample ID file** – single-column `.txt` file (no header) with the FASTQ prefix of each sequenced PCR well  
3. **Ligation barcode error correction dictionaries** – containing ligation barcodes with at most one mismatch  
4. **Reference genome index** for STAR  
5. **Chromosome length file** – `.txt` file with:  
   - Column 1: chromosome names  
   - Column 2: chromosome lengths  
   (must match the reference genome version)  

---

## sgRNA Pipeline

**Main script:**  
`Dual_omics_processing_scripts/sgRNA_modality/perturbFate_dual_omics_sgRNA.main.sh`

### General Workflow
1. Fastq renaming  
2. Amplicon constant region matching  
3. sgRNA UMI parsing  
4. Count matrix reformatting  

### Required Inputs
1. **Fastq files** – named by identifiers of each PCR well  
2. **Fastq sample ID file** – single-column `.txt` file (no header) with the FASTQ prefix of each sequenced PCR well  
3. **Ligation barcode error correction dictionaries** – containing ligation barcodes with at most one mismatch  
4. **Inner i7 barcode error correction dictionary** – containing inner i7 barcodes with at most one mismatch  
5. **sgRNA annotation file** – `.csv` file with three columns:  
   - `gRNA_seq`  
   - `names` (sgRNA identifiers)  
   - `gene_symbol` (target gene)  
6. **sgRNA correction file** – dictionary containing sgRNA sequences with at most one mismatch  

---

## Dual-omics Cell Calling

The `Dual_omics_cell_calling` folder contains functions for parsing count matrices across modalities.

### ATAC Cell Calling

**Script:**  
`Dual_omics_cell_calling/dual_omics_ATAC_cell_calling_function.py`

#### General Workflow
1. Read fragment files exported from SnapATAC2 preprocessing  
2. Combine multiple batches into a unified AnnData object  
3. Perform single-cell peak counting  
4. Calculate TSS enrichment scores  
5. Export outputs  

#### Required Inputs
1. **tsv.gz** – outputs from the ATAC preprocessing step  
2. **peaks_bed** – ATAC peak BED file  

---

### Cross-modality RNA + ATAC + sgRNA Matching

**Script:**  
`Dual_omics_cell_calling/dual_omics_RNA_ATAC_sgRNA_matching_function.py`

#### General Workflow
1. Read single-cell RNA gene count sparse matrix  
2. Match treatment conditions across ligation wells using cell barcodes  
3. Read ATAC cell metadata  
4. Read single-cell sgRNA count sparse matrix  
5. Convert sgRNA cell names to consensus across modalities via barcode matching (`inner i7` + `sgRNA i7` ↔ `RNA/ATAC i7`)  
6. Return unified cell names and outputs  

#### Required Inputs
1. **RNA count matrix (.RDS)** – from RNA preprocessing  
2. **ATAC metadata (.csv)** – from ATAC cell calling  
3. **sgRNA count matrix (.RDS)** – from sgRNA preprocessing  
4. **condition_table (.csv)** – with columns:  
   - `Plate_ID`  
   - `order_bc`  
   - `row`  
   - `col`  
   - `barcode`  
   - `Conditions` (treatment condition names per ligation well)  

---

## Tri-omics Pipeline

The `Tri-omics` folder contains preprocessing scripts for three-modality data. Workflows and required inputs for each modality are similar to Dual-omics.  

**Note:**  
- Nascent RNA is sequenced separately and processed as an independent RNA modality.  
- Cell barcode matching with ATAC/RNA occurs downstream.  

---

## Tri-omics Cell Calling

The `Tri_omics_cell_calling` folder contains functions for parsing and unifying count matrices across three modalities.  

**Note:**  
- Matching requires an additional nascent RNA i7 ↔ preexisting ATAC i7 correspondence file.  
- This dataframe enables integration of nascent RNA with RNA and ATAC modalities.  
