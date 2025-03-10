import os
import copy
import sys
from multiprocessing import Pool
from functools import partial

def from_bed2tsv_SE(input_prefix, bed_dir, output_dir, chrom_len_dict):
    
    ##create the output dir
    if os.path.exists(output_dir) == False:
        os.makedirs(output_dir)
    
    ##generate a dictionary containing all read pairs
    read_pairs_dict = {}
    
    ##open the output file
    out_tsv = open(output_dir + "/" + input_prefix + ".tsv", "w")
    
    bed = open(bed_dir + "/" + input_prefix + ".bed", "r")
    
    print("Start reading the bed file: " + bed_dir + "/" + input_prefix + ".bed")
    
    for each_line in bed:
        elements = each_line.strip().split("\t")
        read_name = elements[3]
        
        ###when the input is paired-end, remove the read1/2 info to make TriSci and SciCAR data compatible
        if "/" in read_name:
            read_name = read_name.split("/")[0]
            
        cell_bc = read_name.split(",")[0]
        cell_name = input_prefix + "." + cell_bc
        
        chrom = elements[0]
        strand = elements[5]
        start = int(elements[1])
        end = int(elements[2])
        
        ###only consider chromosomes with record
        if chrom in chrom_len_dict.keys():
            chrom_len = chrom_len_dict[chrom]

            ###perform correction
            if strand == "+":
                insertion = max((start + 4), 0)
                if insertion < end:
                    out_tsv.write("\t".join([chrom, str(insertion), str(end), cell_name, "1", strand]) + "\n")    
            elif strand == "-":
                insertion = min((end - 5), chrom_len)
                if insertion > start:
                    out_tsv.write("\t".join([chrom, str(start), str(insertion), cell_name, "1", strand]) + "\n")  
            
    bed.close()
    out_tsv.close()
        
    return(input_prefix + ": the bed file has been transformed to tsv.")


def bed2tsv_parallelized_SE(sample_list_dir, bed_dir, output_dir, chrom_len_dir, core):
    
    ##load the sample id
    sample_list_file = open(sample_list_dir, "r")
    sample_list = [x.strip() for x in sample_list_file]
    sample_list_file.close()
    #print(sample_list)
    
    ##load the chrom len and format into a dictionary
    chrom_len_file = open(chrom_len_dir, "r")
    chrom_elm = [x.strip().split("\t") for x in chrom_len_file.readlines()]
    chrom_names = [a[0] for a in chrom_elm]
    chrom_len = [int(a[1]) for a in chrom_elm]
    chrom_len_dict = dict(zip(chrom_names, chrom_len))
    chrom_len_file.close()
    
    # parallelize the function
    p = Pool(processes = int(core))
    func = partial(from_bed2tsv_SE, bed_dir = bed_dir, chrom_len_dict = chrom_len_dict, output_dir = output_dir)
    result = p.map(func, sample_list)
    p.close()
    p.join()
    
    return("All tsv files are generated.")

if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    chrom_len_path = sys.argv[4]
    core = int(sys.argv[5])
    bed2tsv_parallelized_SE(sample_list_dir=sampleID,\
                            bed_dir=input_folder,\
                            output_dir=output_folder,\
                            chrom_len_dir=chrom_len_path,\
                            core=core)

