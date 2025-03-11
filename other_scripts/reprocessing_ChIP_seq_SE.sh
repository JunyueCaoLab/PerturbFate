#! /bin/sh 

############################### modify this part
input_folder="/file/path/to/intermediate_data"
output_folder="/file/path/to/intermediate_data"

fastq_folder="/file/path/to/fq_files"

cores=16
samples=`cat $input_folder/sample_id.txt`
bg_prefix="SRR1234567"

refdir="/file/path/to/Ref/hg38"

picard_path=/file/path/to/picard/picard.jar
################################

###1. cutadapt
mkdir -p $input_folder/trimmed_fq
trimmed_folder="$output_folder/trimmed_fq"
for each_sample in $samples; do
	trim_galore -j $cores --output_dir $trimmed_folder $fastq_folder/${each_sample}.fastq
done

echo "Adapter trimming is done."

###2. map to the reference genome
mapping_folder="$output_folder/mapping"
mkdir -p $mapping_folder

for each_sample in $samples; do
	bowtie2 -p $cores -x $refdir -U $trimmed_folder/${each_sample}_trimmed.fq -S $mapping_folder/${each_sample}.sam
done

echo "Reads mapping is done."

###3. samtools processing and format transformation using bedtools
filtered_folder="$output_folder/filtered_bam"
mkdir -p $filtered_folder
for each_sample in $samples; do
	samtools view -bh -q 2 -@ $cores -S $mapping_folder/${each_sample}.sam > $filtered_folder/${each_sample}.bam
	samtools sort -@ $cores $filtered_folder/${each_sample}.bam > $filtered_folder/${each_sample}.sorted.bam
	#bedtools bamtobed -i $filtered_folder/${each_sample}.sorted.bam > $filtered_folder/${each_sample}.sorted.bed
done

echo "Samtools processing is done."

###4. picard dedup
rmdup_sam_folder=$output_folder/rmdup_bam
rmdup_bed_folder=$output_folder/rmdup_bed
mkdir -p $rmdup_sam_folder
mkdir -p $rmdup_bed_folder

for each_sample in $samples; do
    java -jar $picard_path MarkDuplicates INPUT=$filtered_folder/${each_sample}.sorted.bam OUTPUT=$rmdup_sam_folder/${each_sample}.rmdup.bam REMOVE_DUPLICATES=true \
    ASSUME_SORTED=True METRICS_FILE=$rmdup_sam_folder/${each_sample}.picard_metrics_file.txt VALIDATION_STRINGENCY=LENIENT
    bedtools bamtobed -i $rmdup_sam_folder/${each_sample}.rmdup.bam > $rmdup_bed_folder/${each_sample}.sorted.bed
done

echo "Duplicate removal is done."

###5. call peaks
peak_folder="$output_folder/peak_calling"
mkdir -p $peak_folder

for each_sample in $samples; do
    mkdir -p $peak_folder/$each_sample
    macs2 callpeak -t $rmdup_bed_folder/${each_sample}.sorted.bed \
    -n $each_sample --nomodel --keep-dup all --shift 0 -q 0.05 \
    -B -f BED -g hs --outdir $peak_folder/$each_sample --call-summits
done
