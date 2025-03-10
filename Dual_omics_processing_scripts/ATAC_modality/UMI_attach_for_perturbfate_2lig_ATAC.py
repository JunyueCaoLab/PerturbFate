import subprocess
import sys
from Levenshtein import distance
import gzip
import scipy.spatial.distance
from multiprocessing import Pool
from functools import partial
import pickle
import copy


def UMI_attach_read2_barcode_list(sample, input_folder, output_folder, ligation_1st_barcode_list, ligation_2nd_barcode_list, mismatch_rate = 1):

    # sample = 'sciATAC3_EXP10_01'
    Read1 = input_folder + "/" + sample + ".R1.fastq.gz"
    Read2 = input_folder + "/" + sample + ".R3.fastq.gz"
    Read3 = input_folder + "/" + sample + ".R2.fastq.gz"

    output_file1 = output_folder + "/" + sample + ".R1.fastq.gz"
    output_file2 = output_folder + "/" + sample + ".R2.fastq.gz"

    f1 = gzip.open(Read1)
    f2 = gzip.open(Read2)
    I5_reads = gzip.open(Read3)

    f3 = gzip.open(output_file1, 'wb')
    f4 = gzip.open(output_file2, 'wb')

    line1 = f1.readline()
    line2 = f2.readline()
    line3 = I5_reads.readline()
    total_line = 0
    filtered_line = 0
    ligation_matched_line  = 0
    
    ###generate the error list of ME sequence, use 2 mismatches
    ME = "GAGATGTGTA"
    ME_list = list(ME)
    ME_error_list = set()
    for first_error in range(len(ME)):
        for each_base in ["A", "G", "C", "T", "N"]:
            ME_list_first_error = copy.deepcopy(ME_list)
            ME_list_first_error[first_error] = each_base
            ME_list_first_error_seq = "".join(ME_list_first_error)
            ME_error_list.add(ME_list_first_error_seq)
            if first_error+1 < len(ME):
                for second_error in range(first_error+1,len(ME)):
                    for each_base_again in ["A", "G", "C", "T", "N"]:
                        ME_list_second_error = copy.deepcopy(ME_list_first_error)
                        ME_list_second_error[second_error] = each_base_again
                        ME_list_second_error_seq = "".join(ME_list_second_error)
                        ME_error_list.add(ME_list_second_error_seq)                    
    ME_error_list = list(ME_error_list)

    while (line1):
        total_line += 1
        line1_header = line1
        line2_header = line2
        line1 = f1.readline()
        line2 = f2.readline()
        line3 = I5_reads.readline()
            #print("read1: ", line1)
            # first check if the ligation barcode match with the barcode
        tmp_lig1 = line1[0:10]
        tmp_lig2 = line3[0:10]
        # print(tmp_lig.decode())
        if tmp_lig2.decode() in ligation_2nd_barcode_list and tmp_lig1.decode() in ligation_1st_barcode_list:
            ligation_bc2_match = ligation_2nd_barcode_list[tmp_lig2.decode()]
            ligation_bc1_match = ligation_1st_barcode_list[tmp_lig1.decode()]
            ligation_matched_line += 1
            
            ##match ME sequence (including the additional G) with 2 mismatches allowed
            ME_position = line1[40:50].decode()
            
            if ME_position in ME_error_list:
                filtered_line += 1
                first_line = '@' + ligation_bc1_match + ligation_bc2_match + ',' +line2_header[1:].decode()
                
                f3.write(first_line.encode())
                f4.write(first_line.encode())
                
                ###only retain the trimmed R1
                trimmed_R1_Seq = line1[60:].decode()
                second_line = line2
                f3.write(trimmed_R1_Seq.encode())
                f4.write(second_line)

                line1_third = f1.readline()
                third_line = f2.readline()
                f3.write(line1_third)
                f4.write(third_line)

                line1_fourth = f1.readline()
                four_line = f2.readline()
                trimmed_line1_fourth = line1_fourth[60:].decode()
                f3.write(trimmed_line1_fourth.encode())
                f4.write(four_line)

                line1 = f1.readline()
                line2 = f2.readline()
            else:
                line1 = f1.readline()
                line1 = f1.readline()
                line1 = f1.readline()
                line2 = f2.readline()
                line2 = f2.readline()
                line2 = f2.readline()
        else:      
            line1 = f1.readline()
            line1 = f1.readline()
            line1 = f1.readline()
            line2 = f2.readline()
            line2 = f2.readline()
            line2 = f2.readline()
        line3 = I5_reads.readline() 
        line3 = I5_reads.readline()
        line3 = I5_reads.readline()
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    I5_reads.close()
    print("sample name: %s, total line: %f, barcodes matched line: %f, \
          barcodes and Tn5 features matched line: %f, filter rate: %f" 
            %(sample, total_line,ligation_matched_line, filtered_line, float(filtered_line) / float(total_line)))


def attach_UMI_files(input_folder, sampleID, output_folder, ligation_barcode1_file, ligation_barcode2_file, core):
    
    init_message = '''
    --------------------------start attaching UMI-----------------------------
    input folder: %s
    sample ID: %s
    output_folder: %s
    ligation barcode1 file: %s
    ligation barcode2 file: %s
    ___________________________________________________________________________
    ''' %(input_folder, sampleID, output_folder, ligation_barcode1_file, ligation_barcode2_file)
    
    print(init_message)
    
    
    print("Load ligation barcodes dictionary...")

    barcodes = open(ligation_barcode1_file, "rb")
    ligation_1st_barcode_list = pickle.load(barcodes)
    barcodes.close()

    barcodes = open(ligation_barcode2_file, "rb")
    ligation_2nd_barcode_list = pickle.load(barcodes)
    barcodes.close()

    #for each sample in the sample list, use the read1 file, read2 file, output file
    # and barcode_list to run UMI_attach_read2_barcode_list
    sample_file = open(sampleID)
    sample_list = []
    for line in sample_file:
        sample = line.strip()
        sample_list.append(sample)
    sample_file.close()
    
    # parallele for the functions
    p = Pool(processes = int(core))
    #print("Processing core number: ", core_number)
    func = partial(UMI_attach_read2_barcode_list, input_folder = input_folder, output_folder=output_folder, ligation_1st_barcode_list=ligation_1st_barcode_list, ligation_2nd_barcode_list=ligation_2nd_barcode_list, mismatch_rate = 1)
    #sciRNAseq_count(sample, input_folder, exons, genes, gene_end)
    result = p.map(func, sample_list)
    p.close()
    p.join()
    
    #print the completion message
    com_message = '''~~~~~~~~~~~~~~~UMI attachment done~~~~~~~~~~~~~~~~~~'''
    print(com_message)
    
if __name__ == "__main__":
    input_folder = sys.argv[1]
    sampleID = sys.argv[2]
    output_folder = sys.argv[3]
    ligation_barcode1_file = sys.argv[4]
    ligation_barcode2_file = sys.argv[5]
    core=sys.argv[6]
    attach_UMI_files(input_folder, sampleID, output_folder, ligation_barcode1_file, ligation_barcode2_file, core)
