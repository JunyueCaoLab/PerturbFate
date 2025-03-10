
input_folder=$1
sample=$2
output_folder=$3
barcode_file=$4
cutoff=$5

mismatch=1
python="/file/path/to/bin"
python_script="/file/path/to/perturbfate_scripts/sam_split.py"

$python $python_script $input_folder/$sample.sam $barcode_file $output_folder $cutoff
echo splitting sample done: $sample
