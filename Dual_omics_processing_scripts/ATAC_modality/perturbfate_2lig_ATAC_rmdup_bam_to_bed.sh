input_folder=$1
sample=$2
output_folder=$3

###these bed files are used for snapATAC2 tsv generation, so don't filter chromosomes
samtools sort $input_folder/${sample}.bam | bedtools bamtobed -split -i - | sort -k1,1 -k2,2n - > $output_folder/${sample}.bed

echo "The bed file of ${sample} is generated."
