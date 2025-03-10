#! /bin/bash

#####set up directories and parameters
# define the fastq folder including all fastq files
fastq_folder="/file/path/to/raw_fq"

# define the PCR group sample id for each fastq file
sample_ID="/file/path/to/sampleID.txt"

#define the scripts for single cell splitting and UMI recording
main_script="/file/path/to/perturbfate_2lig_sgRNA_sc_counting.10xCS.58bp_ver.py"
reformatting_script="/file/path/to/gRNA_MM_formatting.R"

# define the output folder
all_output_folder="/file/path/to/intermediate_data"

# define the core number for parallele processing
core=4 # for cell identity and UMI identification

# define the number of gRNA UMI cutoff for keeping cells
cutoff=5

#define the folder of RT barcode dictionary
ligation_barcode_file1="/file/path/to/full_1st_lig_barcodes_96.pickle2"

#define the folder of ligation barcode dictionary
ligation_barcode_file2="/file/path/to/full_2nd_lig_barcodes_384.pickle2"

#define the folder of inner i7 barcode dictionary
inner_i7_bc_file="/file/path/to/inner_i7.pickle2"

#define the folder of gRNA barcode dictionary
gRNA_correction_file=/file/path/to/dualKD_perturbfate_sgRNA_whitelist.pickle2"

#define the folder containing the gRNA annotation file
gRNA_annotation_df="/file/path/to/dualKD_perturbfate_sgRNA_info_df.txt"

#####change the environment

#####change the file names of raw fastq.gz
echo "Changing the name of the fastq files..."
for sample in $(cat $sample_ID); do echo changing name $sample; mv $fastq_folder/*$sample*R1_001.fastq.gz $fastq_folder/$sample.R1.fastq.gz; mv $fastq_folder/*$sample*R2_001.fastq.gz $fastq_folder/$sample.R2.fastq.gz; mv $fastq_folder/*$sample*R3_001.fastq.gz $fastq_folder/$sample.R3.fastq.gz; done

#####run the main script
echo "Start identifying single cells, gRNA sequence and UMI..."
python3 $main_script $fastq_folder $sample_ID $all_output_folder $inner_i7_bc_file $ligation_barcode_file1 $ligation_barcode_file2 $gRNA_correction_file $gRNA_annotation_df $cutoff $core 

#####format the output into a sparse matrix ready for R analysis
report_output_folder=$all_output_folder/gRNA_report
echo "Start exporting the final gRNA expression matrix..."
Rscript $reformatting_script $report_output_folder

echo "All done!"
