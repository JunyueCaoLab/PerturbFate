#! /bin/bash

# The script is for processing PerturbFate Dual-comics RNA modality.

# define the fastq folder including all fastq files
fastq_folder="/file/path/to/raw_fq"

# define the PCR group sample id for each fastq file, the format of cell name is usually: (250307XXXa01).(RT+barcode Ligation barcode)
sample_ID="/file/path/to/sampleID.txt"

# define the output folder
all_output_folder="/file/path/to/intermediate_data"

# define the core number for parallele processing
core=16 # for most steps
samtools_core=8 # for reads filtering and sorting - this number is normally lower than the core number used in other script

# define the number of UMI cutoff for splitting single cell, cells with UMIs less than this number will be discarded
cutoff=200

# define the location of index files for reads alignment with STAR
index="/file/path/to/human_v38/hg38_STAR_index"

# define the gtf file for gene counting
gtf_file="/file/path/to/human_v38/gencode.v38.primary_assembly.annotation.gtf.gz"

# Define the location of the sub script folder
script_folder="/file/path/to/perturbfate_scripts"

#define the bin of python
python_path="/file/path/to/bin"

# define the location of the ligation barcodes
# define the location of the RT barcodes
ligation_1st_bc_list="/file/path/to/full_1st_lig_barcodes_96.pickle2"
ligation_2nd_bc_list="/file/path/to/full_2nd_lig_barcodes_384.pickle2"
barcodes="/file/path/to/all_lig_96x384_combinations.txt"

# define the location of the R script for multi-core processing
R_script=$script_folder/sci3_bash_input_ID_output_core.R
script_folder=$script_folder

now=$(date)
echo "Current time : $now"

############ UMI attach

input_folder=$fastq_folder
output_folder=$all_output_folder/UMI_attach
script="$script_folder/UMI_attach_for_new_coassay3_RNA_shortdT_rN_combined_NexteraR2_tag.py"

echo "Changing the name of the fastq files..."
for sample in $(cat $sample_ID); do echo changing name $sample; mv $input_folder/*$sample*R1_001.fastq.gz $input_folder/$sample.R1.fastq.gz; mv $input_folder/*$sample*R2_001.fastq.gz $input_folder/$sample.R2.fastq.gz; mv $input_folder/*$sample*R3_001.fastq.gz $input_folder/$sample.R3.fastq.gz; done

echo "Attaching barcode and UMI...."
mkdir -p $output_folder

python $script $input_folder $sample_ID $output_folder $ligation_1st_bc_list $ligation_2nd_bc_list $core
echo "Barcode transformed"

################# generate a dT rN separated sample ID list
for sample in $(cat $sample_ID); do 
echo ${sample}_dT >> $all_output_folder/dT_sampleID.txt
echo ${sample}_rN >> $all_output_folder/rN_sampleID.txt
echo ${sample}_dT >> $all_output_folder/dT_rN_sampleID.txt
echo ${sample}_rN >> $all_output_folder/dT_rN_sampleID.txt
done

dT_sample_ID=$all_output_folder/dT_sampleID.txt
rN_sample_ID=$all_output_folder/rN_sampleID.txt
dT_rN_sample_ID=$all_output_folder/dT_rN_sampleID.txt

################# Trimming the read2
echo
echo "Start trimming the read2 file..."
echo $(date)

###trim dT reads
trimmed_fastq=$all_output_folder/trimmed_fastq
UMI_attached_R2=$all_output_folder/UMI_attach
bash_script=$script_folder/perturbfate_RNA_shortT_trim.sh

Rscript $R_script $bash_script $UMI_attached_R2 $dT_sample_ID $trimmed_fastq $core

echo dT trimming is done!

###trim rN reads
trimmed_fastq=$all_output_folder/trimmed_fastq
UMI_attached_R2=$all_output_folder/UMI_attach
bash_script=$script_folder/perturbfate_RNA_randomN_trim.sh

Rscript $R_script $bash_script $UMI_attached_R2 $rN_sample_ID $trimmed_fastq $core

echo rN trimming is done!

############merge dT and rN files

trimmed_fastq=$all_output_folder/trimmed_fastq
merged_trimmed_fastq=$all_output_folder/merged_trimmed_fastq

mkdir -p $merged_trimmed_fastq

for each_dT_rN in $(cat $sample_ID); do 
zcat $trimmed_fastq/$each_dT_rN*gz > $merged_trimmed_fastq/${each_dT_rN}.fq
pigz -p $core $merged_trimmed_fastq/${each_dT_rN}.fq
done

###rm -rf $trimmed_fastq

echo Merging trimmed fastq is done!

############Align the reads with STAR, filter the reads based on qscore, and remove duplicates based on UMI sequence and tagmentation site

#define the output folder for mapping
input_folder=$merged_trimmed_fastq
STAR_output_folder=$all_output_folder/STAR_alignment
filtered_sam_folder=$all_output_folder/filtered_sam
rmdup_sam_folder=$all_output_folder/rmdup_sam

#align read2 to the index file using STAR
echo "Start alignment using STAR..."
echo input folder: $input_folder
echo sample ID file: $sample_ID
echo index file: $index
echo output_folder: $STAR_output_folder
#make the output folder
mkdir -p $STAR_output_folder

#remove the index from the memory
STAR --genomeDir $index --genomeLoad Remove
#start the alignment
for sample in $(cat $sample_ID); do echo Aligning $sample;STAR --runThreadN $core --outSAMstrandField intronMotif --genomeDir $index --readFilesCommand zcat --readFilesIn $input_folder/$sample*gz --outFileNamePrefix $STAR_output_folder/$sample --genomeLoad LoadAndKeep; done
#remove the index from the memory
STAR --genomeDir $index --genomeLoad Remove
echo "All alignment done."

##############Filter and sort the sam file

echo
echo "Start filter and sort the sam files..."
echo input folder: $STAR_output_folder
echo output folder: $filtered_sam_folder
bash_script=$script_folder/sci3_filter.sh
Rscript $R_script $bash_script $STAR_output_folder $sample_ID $filtered_sam_folder $samtools_core

echo
echo "Start removing duplicates..."
echo input folder: $filtered_sam_folder
echo output folder: $rmdup_sam_folder
mkdir -p $rmdup_sam_folder
mismatch=0
bash_script=$script_folder/sci3_rmdup_nomismatch.sh
Rscript $R_script $bash_script $filtered_sam_folder $sample_ID $rmdup_sam_folder $core $mismatch


################# split the sam file based on the barcode, and mv the result to the report folder
sam_folder=$all_output_folder/rmdup_sam
output_folder=$all_output_folder/sam_splitted

echo
echo "Start splitting the dT sam file..."
echo samfile folder: $sam_folder
echo sample list: $sample_ID
echo ouput folder: $output_folder
echo barcode file: $barcodes
echo cutoff value: $cutoff

bash_script=$script_folder/sci3_split.sh
Rscript $R_script $bash_script $sam_folder $sample_ID $output_folder $core $barcodes $cutoff

cat $output_folder/*sample_list.txt>$output_folder/All_samples.txt
cp $output_folder/All_samples.txt $all_output_folder/barcode_samples.txt

mkdir -p $all_output_folder/report/barcode_read_distribution
mv $output_folder/*.txt $all_output_folder/report/barcode_read_distribution/
mv $output_folder/*.png $all_output_folder/report/barcode_read_distribution/
echo
echo "All sam file splitted."

################### calculate the reads number
fastq_folder=$fastq_folder
trimmed_folder=$merged_trimmed_fastq
UMI_attach=$UMI_attached_R2
alignment=$STAR_output_folder
filtered_sam=$filtered_sam_folder
rm_dup_sam=$rmdup_sam_folder
report_folder=$all_output_folder/report/read_num
echo
echo "Start calculating the reads number..."
#make the report folder
mkdir -p $report_folder
#calculate the read number and output the read number into the report folder
echo sample,total reads,after filtering barcode,after trimming,uniquely aligned reads,After remove duplicates > $report_folder/read_number.csv
for sample in $(cat $sample_ID); do 
echo calculating $sample
echo $sample,$(expr $(zcat $fastq_folder/$sample*R2*.gz|wc -l) / 4),$(expr $(zcat $UMI_attach/$sample*R2*.gz|wc -l) / 4),$(expr $(zcat $trimmed_folder/$sample*.gz|wc -l) / 4),$(cat $filtered_sam/$sample*.sam|grep -v "^@"|wc -l),$(cat $rm_dup_sam/$sample*.sam|grep -v "^@"|wc -l) >> $report_folder/read_number.csv
done
echo "Read number calculation is done."


################# gene count
output_folder=$all_output_folder/report/gene_count/
input_folder=$all_output_folder/sam_splitted
script=$script_folder/TriSci_count_rN_dT.merged.py
sample_ID=$all_output_folder/barcode_samples.txt

echo "Start the gene count...."
$python_path/python $script $gtf_file $input_folder $sample_ID $core

echo "Make the output folder and transfer the files..."
mkdir -p $output_folder
find $input_folder/ -name "*.count" -print0 | sort -z | xargs -r0 cat > $output_folder/count.MM
find $input_folder/ -name "*.count" -delete
find $input_folder/ -name "*.report" -print0 | sort -z | xargs -r0 cat > $output_folder/report.MM
find $input_folder/ -name "*.report" -delete
mv $input_folder/*_annotate.txt $output_folder/
echo "All output files are transferred~"

################# reformat output files for R downstream processing
R_script=$script_folder/gene_count_processing_sciRNAseq_exon_intron.R
Rscript $R_script $all_output_folder/report
conda deactivate

mkdir $all_output_folder/report/Log_files
mv Log.out $all_output_folder/report/Log_files/
mv Log.progress.out $all_output_folder/report/Log_files/
mv Aligned.out.sam $all_output_folder/report/Log_files/
mv _STARtmp $all_output_folder/report/Log_files/

now=$(date)
echo "Current time : $now"