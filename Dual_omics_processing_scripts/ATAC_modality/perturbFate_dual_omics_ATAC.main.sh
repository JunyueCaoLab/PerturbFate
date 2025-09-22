#! /bin/bash

# The script is for processing PerturbFate Dual-comics ATAC modality. The output tsv file is fully compatible with SnapATAC2

########################################################## MODIFY THIS
fastq_folder="/file/path/to/raw_fq"
sample_ID="/file/path/to/sampleID.txt"
all_output_folder="/file/path/to/intermediate_data"
core=16

###Choose the human/mouse index, GTF file, chromosome length file, snapATAC_genome_name
index="/file/path/to/human_v38/hg38_STAR_index"
#the chromosome length file should be consistent with the reference we used for mapping. That's the key thing.
chromosome_length_file='/file/path/to/human_v38/hsa_ref_chr_length.txt'

script_folder="/file/path/to/perturbfate_scripts"

R_script="/file/path/to/perturbfate_scripts/sci3_bash_input_ID_output_core.R"

ligation_1st_bc_list="/file/path/to/full_1st_lig_barcodes_96.pickle2"
ligation_2nd_bc_list="/file/path/to/full_2nd_lig_barcodes_384.pickle2"
barcodes="/file/path/to/all_lig_96x384_combinations.txt"
###########################################################


################ Reformat file names
input_folder=$fastq_folder
output_folder=$all_output_folder/UMI_attach
echo "Changing the name of the fastq files..."
for sample in $(cat $sample_ID); do echo changing name $sample; mv $input_folder/*$sample*R1*fastq.gz $input_folder/$sample.R1.fastq.gz; mv $input_folder/*$sample*R2*fastq.gz $input_folder/$sample.R2.fastq.gz; mv $input_folder/*$sample*R3*fastq.gz $input_folder/$sample.R3.fastq.gz; done

################ Barcode extraction
script="$script_folder/UMI_attach_for_perturbfate_2lig_ATAC.py"
echo "Attaching barcode...."
mkdir -p $output_folder

python $script $input_folder $sample_ID $output_folder $ligation_1st_bc_list $ligation_2nd_bc_list $core
echo "Barcode transformed."

################ Trimming the reads
echo
echo "Start trimming the read file..."
mkdir -p $all_output_folder/trimmed_fastq
trimmed_fastq=$all_output_folder/trimmed_fastq

input_folder=$all_output_folder/UMI_attach
output_folder=$trimmed_fastq

bash_script=$script_folder/perturbfate_2lig_ATAC_trim.sh
Rscript $R_script $bash_script $input_folder $sample_ID $output_folder $core

################# align the reads with STAR, filter the reads based on q > 30, and remove duplicates based on the genomic coordinates
input_folder=$trimmed_fastq
STAR_output_folder=$all_output_folder/STAR_alignment
filtered_sam_folder=$all_output_folder/filtered_sam
rmdup_sam_folder=$all_output_folder/rmdup_sam

echo "Start alignment using STAR..."
mkdir -p $STAR_output_folder
STAR --genomeDir $index --genomeLoad Remove

for sample in $(cat $sample_ID); do echo aligning $sample; STAR --runThreadN $core --outSAMstrandField intronMotif --genomeDir $index --readFilesCommand zcat --readFilesIn $input_folder/$sample*R1*.gz $input_folder/$sample*R2*.gz --outFileNamePrefix $STAR_output_folder/$sample --genomeLoad LoadAndKeep --outFilterMultimapNmax 1; done
#remove the index from the memory
STAR --genomeDir $index --genomeLoad Remove
echo "All alignment done."

################## Filter aligned reads 
echo
echo "Start filter and sort the sam files..."
mkdir -p $filtered_sam_folder

for sample in $(cat $sample_ID); do echo Filtering $sample; samtools view -bh -q 30 -F 4 -@ $core -F 0x100 $STAR_output_folder/$sample*.sam|samtools sort -@ $core -|samtools view -h - >$filtered_sam_folder/$sample.sam; done

################### remove duplicates using Picard
echo Remove duplicates
rmdup_sam_folder=$all_output_folder/rmdup_sam
mkdir -p $rmdup_sam_folder

bash_script=$script_folder/perturbfate_2lig_ATAC_rmdup.sh
Rscript $R_script $bash_script $filtered_sam_folder $sample_ID $rmdup_sam_folder $core
echo "Deduplication completed."

#################### calculate the reads number in each file
report_folder=$all_output_folder/report
mkdir -p $report_folder
fastq_folder=$fastq_folder
trimmed_folder=$trimmed_fastq
alignment=$STAR_output_folder
filtered_sam=$filtered_sam_folder
rm_dup_sam=$rmdup_sam_folder
echo sample,Total reads,Barcode matching,After trimming,Uniquely aligned reads,Filtered reads,After remove duplicates>$report_folder/read_number.csv
for sample in $(cat $sample_ID); do echo calculating $sample; echo $sample,$(expr $(zcat $fastq_folder/$sample.R[13].fastq.gz|wc -l) / 4),$(expr $(zcat $all_output_folder/UMI_attach/$sample*.gz|wc -l) / 4),$(expr $(zcat $trimmed_folder/$sample*.gz|wc -l) / 4),$(samtools view $alignment/$sample*.sam|wc -l),$(samtools view $filtered_sam/$sample.sam|wc -l),$(samtools view $rm_dup_sam/$sample.bam|wc -l)>>$report_folder/read_number.csv; done
echo "Read number calculation is done."

#################### Merge bam files
sam_folder=$all_output_folder/rmdup_sam
merged_bam_folder=$all_output_folder/merged_bam
mkdir -p $merged_bam_folder

samtools merge -f -@ $core $merged_bam_folder/merged.bam $sam_folder/*.bam
samtools sort -@ $core -n $merged_bam_folder/merged.bam > $merged_bam_folder/merged.sorted.bam

#################### Add PCR sample name to the head of read names for tsv generation
sam_folder=$all_output_folder/rmdup_sam
output_folder=$all_output_folder/bed_file
mkdir -p $output_folder

echo "Start generating the bed files for snapATAC2..."
bash_script=$script_folder/perturbfate_2lig_ATAC_rmdup_bam_to_bed.sh
Rscript $R_script $bash_script $sam_folder $sample_ID $output_folder $core

##################### make the tsv file for snapATAC2
input_folder=$all_output_folder/bed_file
output_folder=$all_output_folder/tsv_file
mkdir -p $output_folder

echo "Start generating the tsv files for snapATAC2..."

bed_to_tsv_script=$script_folder/perturbfate_bed_to_tsv_SE.py
python $bed_to_tsv_script $input_folder $sample_ID $output_folder $chromosome_length_file $core

echo "tsv file is generated!"

##merge tsv files
snapATAC2_folder=$output_folder/snapATAC2
mkdir -p $snapATAC2_folder

echo "Start merging and indexing tsv files..."
cat $output_folder/*tsv | sort -k1,1 -k2,2n > $snapATAC2_folder/merged.tsv
bgzip -@ $core -f $snapATAC2_folder/merged.tsv
tabix -0 -b 2 -e 3 -f $snapATAC2_folder/merged.tsv.gz

echo "tsv file is merged and indexed!"

