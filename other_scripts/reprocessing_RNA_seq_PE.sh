#! /bin/sh 

############################## modify this part
input_folder="/file/path/to/intermediate_data"
output_folder="/file/path/to/intermediate_data"

fastq_folder="/file/path/to/fq_files"

cores=16
samples=`cat $input_folder/sample_id.txt`

index="/file/path/to/human_v38/hg38_STAR_index"

picard_path=/file/path/to/picard/picard.jar

annotation_gtf="/file/path/to/human_v38/gencode.v38.primary_assembly.annotation.gtf"
##############################

###1. cutadapt
trimmed_folder="$output_folder/trimmed_fq"
mkdir -p $trimmed_folder

for each_sample in $samples; do
    trim_galore -j $cores --paired --gzip $fastq_folder/${each_sample}_1.fastq $fastq_folder/${each_sample}_2.fastq -o $trimmed_folder
done

echo "Adapter trimming is done."

###2. map to the reference genome
refdir="/rugpfs/fs0/cao_lab/scratch/zxu/tools/ATAC-pipe/Data/Ref/hg38/hg38"
mapping_folder="$output_folder/mapping"
mkdir -p $mapping_folder

#remove the index from the memory
STAR --genomeDir $index --genomeLoad Remove

for each_sample in $samples; do
    STAR --runThreadN $cores --outSAMstrandField intronMotif --genomeDir $index --readFilesCommand zcat --readFilesIn $trimmed_folder/${each_sample}_1_val_1.fq.gz $trimmed_folder/${each_sample}_2_val_2.fq.gz --outFileNamePrefix $mapping_folder/${each_sample} --genomeLoad LoadAndKeep
done

#remove the index from the memory
STAR --genomeDir $index --genomeLoad Remove

echo "Reads mapping is done."

###3. samtools processing and format transformation using bedtools
filtered_folder="$output_folder/filtered_bam"
mkdir -p $filtered_folder
for each_sample in $samples; do
	samtools view -bh -q 10 -@ $cores -S $mapping_folder/${each_sample}*sam > $filtered_folder/${each_sample}.bam
	samtools sort -@ $cores $filtered_folder/${each_sample}.bam > $filtered_folder/${each_sample}.sorted.bam
    samtools index -@ $cores $filtered_folder/${each_sample}.sorted.bam
done

echo "Samtools processing is done."

###4. deduplication
rmdup_sam_folder=$output_folder/rmdup_bam
mkdir -p $rmdup_sam_folder

for each_sample in $samples; do
    java -jar $picard_path MarkDuplicates INPUT=$filtered_folder/${each_sample}.sorted.bam OUTPUT=$rmdup_sam_folder/${each_sample}.rmdup.bam REMOVE_DUPLICATES=true \
    ASSUME_SORTED=True METRICS_FILE=$rmdup_sam_folder/${each_sample}.picard_metrics_file.txt VALIDATION_STRINGENCY=LENIENT
done

echo "Deduplication processing is done."

###5. gene counting
count_dir=$output_folder/gene_counts
mkdir -p $count_dir

for each_sample in $samples; do
    featureCounts -p -a $annotation_gtf -T $cores -o $count_dir/$each_sample $rmdup_sam_folder/${each_sample}.rmdup.bam -g gene_name 
done