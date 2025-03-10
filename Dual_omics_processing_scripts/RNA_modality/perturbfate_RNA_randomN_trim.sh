input_folder=$1
sample=$2
output_folder=$3

echo Trimming sample: $sample
trim_galore $input_folder/$sample.R2.fastq.gz -a CACAAGTATGCAGCGCG --three_prime_clip_R1 14 -o $output_folder
#module unload python/2.7.3
echo Trimming $sample done.
