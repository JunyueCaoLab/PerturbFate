# PerturbFate
A combinatorial-indexing single-cell platform for concurrent profiling of chromatin accessibility, nascent and pre-existing RNA, or steady-state transcriptomic phenotypes, along with sgRNA identities

# Dual-omics Preprocessing Pipelines

The `Dual-omics` folder contains the main preprocessing script and daughter scripts embedded for processing FASTQ files into crude count matrices.

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
2. **Fastq sample ID file** – a single-column `.txt` file (no header) containing the FASTQ prefix of each sequenced PCR well  
3. **Ligation barcode error correction dictionaries** – dictionaries containing ligation barcodes with at most one mismatch  
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
2. **Fastq sample ID file** – a single-column `.txt` file (no header) containing the FASTQ prefix of each sequenced PCR well  
3. **Ligation barcode error correction dictionaries** – dictionaries containing ligation barcodes with at most one mismatch  
4. **Reference genome index** for STAR  
5. **Chromosome length file** – `.txt` file with chromosome names in the first column and lengths in the second column (must match reference genome version)  

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
2. **Fastq sample ID file** – a single-column `.txt` file (no header) containing the FASTQ prefix of each sequenced PCR well  
3. **Ligation barcode error correction dictionaries** – dictionaries containing ligation barcodes with at most one mismatch  
4. **Inner i7 barcode error correction dictionary** – dictionary containing inner i7 barcodes with at most one mismatch  
5. **sgRNA annotation file** – `.csv` file with three columns:  
   - `gRNA_seq`  
   - `names` (sgRNA identifiers)  
   - `gene_symbol` (target gene)  
6. **sgRNA correction file** – dictionary containing sgRNA sequences with at most one mismatch  

---
